use std::{cmp::min, fmt::Display, ops::Mul};
use crate::glfem::color_maps;
use cgmath::num_traits::Float;
use glyphon::Color;
use wgpu::{util::DeviceExt, ColorTargetState, MultisampleState, RenderPass};

use super::color_maps::{MAGMA_DATA, PLASMA_DATA, TURBO_DATA, VIRIDIS_DATA};

pub struct MeshRenderer {
    triangle_pipeline: wgpu::RenderPipeline,
    line_pipeline: wgpu::RenderPipeline,
    index_buffer: wgpu::Buffer,
    num_indices: u32,
    vertices: Vec<Vertex>,
    lines: Vec<(u32, u32)>,
    mesh_visible: bool,
    lines_visible: bool,
}

#[derive(Debug, Clone, Copy)]
pub enum ColorMap {
    JET,
    MAGMA,
    TURBO,
    VIRIDIS,
    PLASMA,
    COOLWARM,
    GREYSCALE,
}

impl Display for ColorMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ColorMap::JET => write!(f, "Jet"),
            ColorMap::MAGMA => write!(f, "Magma"),
            ColorMap::TURBO => write!(f, "Turbo"),
            ColorMap::VIRIDIS => write!(f, "Viridis"),
            ColorMap::PLASMA => write!(f, "Plasma"),
            ColorMap::COOLWARM => write!(f, "Coolwarm"),
            ColorMap::GREYSCALE => write!(f, "Greyscale"),
        }
    }
}

impl ColorMap {
    pub fn color(&self, value: f32) -> [f32; 3] {
        match self {
            ColorMap::JET => Self::jet(value),
            ColorMap::MAGMA => Self::map(value, MAGMA_DATA),
            ColorMap::TURBO => Self::map(value, TURBO_DATA),
            ColorMap::VIRIDIS => Self::map(value, VIRIDIS_DATA),
            ColorMap::PLASMA => Self::map(value, PLASMA_DATA),
            ColorMap::COOLWARM => Self::coolwarm(value),
            ColorMap::GREYSCALE => Self::greyscale(value),
        }
    }

    pub fn list() -> Vec<ColorMap> {
        vec![
            ColorMap::JET,
            ColorMap::MAGMA,
            ColorMap::TURBO,
            ColorMap::VIRIDIS,
            ColorMap::PLASMA,
            ColorMap::COOLWARM,
            ColorMap::GREYSCALE,
        ]
    }

    pub fn default() -> Self {
        ColorMap::JET
    }

    fn map(value: f32, data: [[f32; 3]; 256]) -> [f32; 3] {
        let r = data[(value * 255.0) as usize][0];
        let g = data[(value * 255.0) as usize][1];
        let b = data[(value * 255.0) as usize][2];
        [r, g, b]
    }

    fn coolwarm(v: f32) -> [f32; 3] {
        if v <= 0.5 {
            let t = v * 2.0;
            [t, t, 1.0]
        } else {
            let t = (v - 0.5) * 2.0;
            [1.0, 1.0 - t, 1.0 - t]
        }
    }

    fn greyscale(v: f32) -> [f32; 3] {
        [v, v, v]
    }

    fn jet(v: f32) -> [f32; 3] {
        let mut v = v;
        
        let c = 
        if v < 0.25 {
            [0.0, 4.0 * v, 1.0]
        } else if v < 0.5 {
            [0.0, 1.0, 1.0 - 4.0 * (v - 0.25)]
        } else if v < 0.75 {
            [4.0 * (v - 0.5), 1.0, 0.0]
        } else {
            [1.0, 1.0 - 4.0 * (v - 0.75), 0.0]
        };

        c.map(|x| x.clamp(0.0, 1.0))
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Vertex {
    position: [f32; 3],
    value: f32,
}
#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct DrawableVertex {
    position: [f32; 3],
    color: [f32; 3],
}

impl Vertex {

    pub fn new(position: [f32; 3], value: f32) -> Self {
        Self { position, value }
    }
   
}

impl DrawableVertex {
    const NUMBER_OF_COLORS: u16 = 50;

    const ATTRIBS: [wgpu::VertexAttribute; 2] =
        wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32x3];

    pub fn desc() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Self>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &Self::ATTRIBS,
        }
    }

    pub fn from(position: [f32; 3], value: f32, color_map: ColorMap) -> Self {
        let color = color_map.color(value);
        Self {
            position: position,
            color,
        }
    }

    pub fn list_from(vertices: &Vec<Vertex>, color_map: ColorMap) -> Vec<Self> {
        vertices.iter().map(|v| Self::from(v.position, v.value, color_map)).collect()
    }

    pub fn gazou_color(mut value: f32) -> [f32; 3] {
        if value > 1.0 {
            value = 1.0;
        }
        value = value * (Self::NUMBER_OF_COLORS as f32);
        value = value - 0.00000001;
        value = value / (Self::NUMBER_OF_COLORS as f32 - 1.0);
    
        value = (1.0 - value).clamp(0.0, 1.0);
        
        let r = 3.5 * (1.0 - value) * (1.0 - value);
        let g = (1.0 - value) * (value) * 3.5;
        let b = value * value;
        [r, g, b]
    }
}

impl MeshRenderer {
    pub fn new(
        device: &wgpu::Device,
        config: &wgpu::SurfaceConfiguration,
        vertices: Vec<Vertex>,
        indices: Vec<u16>,
        lines: Vec<(u32, u32)>,
        multisample: MultisampleState,
    ) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Mesh Shader"),
            source: wgpu::ShaderSource::Wgsl(
                include_str!("../../assets/shaders/shader.wgsl").into(),
            ),
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Pipeline Layout"),
            bind_group_layouts: &[],
            push_constant_ranges: &[],
        });

        let triangle_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Render Pipeline"),
            layout: Some(&pipeline_layout),
            cache: None,
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_main"),
                compilation_options: Default::default(),
                buffers: &[DrawableVertex::desc()],
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_main"),
                compilation_options: Default::default(),
                targets: &[Some(wgpu::ColorTargetState {
                    format: config.format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None, //Some(wgpu::Face::Back),
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
            },
            depth_stencil: None,
            multisample,
            multiview: None,
        });

        let line_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Line Pipeline"),
            layout: Some(&pipeline_layout),
            cache: None,
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_main"),
                compilation_options: Default::default(),
                buffers: &[DrawableVertex::desc()],
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_main"),
                compilation_options: Default::default(),
                targets: &[Some(wgpu::ColorTargetState {
                    format: config.format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None, //Some(wgpu::Face::Back),
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
            },
            depth_stencil: None,
            multisample,
            multiview: None,
        });


        let index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Index Buffer"),
            contents: bytemuck::cast_slice(&indices),
            usage: wgpu::BufferUsages::INDEX,
        });

        // Normalize vertices to [-1, 1]

        let mut normed_vertices = Vec::new();
        let x_max = vertices.iter().map(|v| v.position[0]).fold(0.0, |a, b| a.max(b));
        let y_max = vertices.iter().map(|v| v.position[1]).fold(0.0, |a, b| a.max(b));
        let x_min = vertices.iter().map(|v| v.position[0]).fold(0.0, |a, b| a.min(b));
        let y_min = vertices.iter().map(|v| v.position[1]).fold(0.0, |a, b| a.min(b));

        for vertex in vertices {
            let mut normed_vertex = Vertex::new([0.0, 0.0, 0.0], vertex.value);
            normed_vertex.position[0] = 2.0*(vertex.position[0] - x_min) / (x_max - x_min) - 1.0;
            normed_vertex.position[1] = 2.0*(vertex.position[1] - y_min) / (y_max - y_min) - 1.0;
            normed_vertices.push(normed_vertex);
        }

        Self {
            triangle_pipeline,
            line_pipeline,
            index_buffer,
            num_indices: indices.len() as u32,
            vertices: normed_vertices,
            lines,
            mesh_visible: true,
            lines_visible: true,
        }
    }

    const SCALE: f32 = 0.5;

    pub fn transform_vertices(&self, vertices: &Vec<Vertex>, aspect_ratio: f32, user_scale: f32, translation: (f32, f32)) -> Vec<Vertex> {
        let aspect_ratio_x = if aspect_ratio > 1.0 { 1.0 } else { aspect_ratio };
        let aspect_ratio_y = if aspect_ratio < 1.0 { 1.0 } else { aspect_ratio };
        let scale = Self::SCALE * user_scale;
        vertices.iter().map(|v| Vertex::new([v.position[0] * aspect_ratio_x * scale + translation.0, v.position[1] * aspect_ratio_y * scale + translation.1, v.position[2]], v.value)).collect()
    }

    pub fn draw_lines(&self, render_pass: &mut RenderPass, device: &wgpu::Device, aspect_ratio: f32, scale: f32, translation: (f32, f32)) {
            
        let mut line_vertices = Vec::new();
            
            for (i, j) in &self.lines {
                let v1 = &self.vertices[*i as usize];
                let v2 = &self.vertices[*j as usize];
                let dx = v2.position[0] - v1.position[0];
                let dy = v2.position[1] - v1.position[1];
                let length = (dx*dx + dy*dy).sqrt();
                let nx = dy / length;
                let ny = -dx / length;
                let delta = 0.002/2.0;
                let v1p = Vertex::new([v1.position[0] - nx * delta, v1.position[1] - ny * delta, v1.position[2]], 0.0);
                let v2p = Vertex::new([v2.position[0] - nx * delta, v2.position[1] - ny * delta, v2.position[2]], 0.0);
                let v3p = Vertex::new([v1.position[0] + nx * delta, v1.position[1] + ny * delta, v1.position[2]], 0.0);
                let v4p = Vertex::new([v2.position[0] + nx * delta, v2.position[1] + ny * delta, v2.position[2]], 0.0);
                line_vertices.push(v1p);
                line_vertices.push(v2p);
                line_vertices.push(v3p);
                line_vertices.push(v2p);
                line_vertices.push(v3p);
                line_vertices.push(v4p);
            }

            let line_vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Line Vertex Buffer"),
                contents: bytemuck::cast_slice(&DrawableVertex::list_from(&self.transform_vertices(&line_vertices, aspect_ratio, scale, translation), ColorMap::GREYSCALE)),
                usage: wgpu::BufferUsages::VERTEX,
            });

            let line_index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Line Index Buffer"),
                contents: bytemuck::cast_slice(&(0..(line_vertices.len() as u32)).collect::<Vec<u32>>()),
                usage: wgpu::BufferUsages::INDEX,
            });

            render_pass.set_pipeline(&self.line_pipeline);
            render_pass.set_vertex_buffer(0, line_vertex_buffer.slice(..));
            render_pass.set_index_buffer(line_index_buffer.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..(line_vertices.len() as u32), 0, 0..1);
    }

    fn color_vertices(&self, color_map: ColorMap) {
        
    }

    pub fn toggle_mesh(&mut self) {
        self.mesh_visible = !self.mesh_visible;
    }

    pub fn toggle_lines(&mut self) {
        self.lines_visible = !self.lines_visible;
    }

    pub fn render(
        &self,
        render_pass: &mut RenderPass,
        device: &wgpu::Device,
        surface_config: &wgpu::SurfaceConfiguration,
        scale: f32,
        translation: (f32, f32),
        color_map: Option<ColorMap>,
    ) {

        let h = surface_config.height as f32;
        let w = surface_config.width as f32;
        let aspect_ratio = h/w;

        let v_max = self.vertices.iter().map(|v| v.value).fold(0.0, |a, b| a.max(b));
        let v_min = self.vertices.iter().map(|v| v.value).fold(0.0, |a, b| a.min(b));

        let color_map = color_map.unwrap_or(ColorMap::default());
        
        let vertices = DrawableVertex::list_from(&self.transform_vertices(&self.vertices, aspect_ratio, scale, translation), color_map);
        if self.mesh_visible {
            let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Vertex Buffer"),
                contents: bytemuck::cast_slice(&vertices),
                usage: wgpu::BufferUsages::VERTEX,
            });

            render_pass.set_pipeline(&self.triangle_pipeline);
            render_pass.set_vertex_buffer(0, vertex_buffer.slice(..));
            render_pass.set_index_buffer(self.index_buffer.slice(..), wgpu::IndexFormat::Uint16);
            render_pass.draw_indexed(0..self.num_indices, 0, 0..1);

        }

        if self.lines_visible {
            self.draw_lines(render_pass, device, aspect_ratio, scale, translation);
        }       
    }
}
