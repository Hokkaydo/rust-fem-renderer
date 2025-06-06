use glyphon::Color;
use glyphon::{
    fontdb::Database, Cache, FontSystem, Resolution, SwashCache, TextAtlas, TextRenderer, Viewport,
};
use winit::event::MouseButton;
use std::ptr::null;
use std::sync::Arc;
use wgpu::{
    CommandEncoderDescriptor, CompositeAlphaMode, DeviceDescriptor, Instance, InstanceDescriptor, MultisampleState, Operations, PresentMode, RenderPassColorAttachment, RenderPassDescriptor, RequestAdapterOptions, SurfaceConfiguration, TextureFormat, TextureUsages, TextureViewDescriptor
};
use winit::{dpi::LogicalSize, event::WindowEvent, event_loop::EventLoop, window::Window};

use crate::fem::FemProblem;
use crate::glfem::mesh_renderer::{ColorMap, Vertex};
use crate::glfem::{mesh_renderer, text};

fn generate_colored_grid(problem: &FemProblem) -> (Vec<Vertex>, Vec<u16>) {
    let geo = problem.geometry.clone();
    let nodes = &geo.nodes;
    let nnodes = nodes.numbers.len();


    let mut vertices = vec![Vertex::new([0.0, 0.0, 0.0], 0.0); nnodes];
    let indices = Vec::from(geo.elements.elements.to_vec().iter().map(|&x| x as u16).collect::<Vec<u16>>());

    let u = problem.solution.chunks_exact(2).map(|c| (c[0] * c[0] + c[1] * c[1]).sqrt()).take(nnodes).collect::<Vec<f64>>();
    let u_min = u.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let u_max = u.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    let u_scaled = u.iter().map(|&u| ((u-u_min)/(u_max-u_min)) as f32).collect::<Vec<f32>>();
    
    for idx in 0..geo.nodes.numbers.len()/2 {
        let node = idx as usize;
        let x = geo.nodes.x[node] as f32 + 2e4*problem.solution[2*node] as f32;
        let y = geo.nodes.y[node] as f32 + 2e4*problem.solution[2*node+1] as f32;
        let z = 0.0;
        vertices[node] = Vertex::new([x, y, z], u_scaled[node]);
    }

    (vertices, indices)
}

pub fn run(start_function: fn() -> FemProblem) {
    let event_loop = EventLoop::new().unwrap();
    event_loop
        .run_app(&mut Application {
            window_state: None,
            start_function,
        })
        .unwrap();
}

const FONT_PATH: &str = "/usr/share/fonts/TTF/JetBrainsMono-Bold.ttf";
//"/usr/share/fonts/TTF/JetBrainsMono-Bold.ttf";//"./assets/fonts/font.ttf";

struct DragState {
    translation: (f32, f32),
    last_cursos_pos: (f32, f32),
    is_dragging: bool,
}
impl DragState {
    fn new() -> Self {
        Self {
            translation: (0.0, 0.0),
            last_cursos_pos: (0.0, 0.0),
            is_dragging: false,
        }
    }
}
struct WindowState {
    device: wgpu::Device,
    queue: wgpu::Queue,
    surface: wgpu::Surface<'static>,
    surface_config: SurfaceConfiguration,
    scale: f32,
    drag_state: DragState,
    font_system: FontSystem,
    swash_cache: SwashCache,
    viewport: glyphon::Viewport,
    atlas: glyphon::TextAtlas,
    text_renderer: TextRenderer,
    texts: text::Texts,
    mesh_renderer: mesh_renderer::MeshRenderer,
    color_map: usize,
    multisample: MultisampleState,
    // Make sure that the winit window is last in the struct so that
    // it is dropped after the wgpu surface is dropped, otherwise the
    // program may crash when closed. This is probably a bug in wgpu.
    window: Arc<Window>,
}

impl WindowState {
    async fn new(window: Arc<Window>, vertices: Vec<Vertex>, indices: Vec<u16>, lines: Vec<(u32, u32)>, problem: FemProblem) -> Self {
        let physical_size = window.inner_size();

        let instance = Instance::new(&InstanceDescriptor::default());
        let adapter = instance
            .request_adapter(&RequestAdapterOptions::default())
            .await
            .unwrap();
        let (device, queue) = adapter
            .request_device(&DeviceDescriptor::default(), None)
            .await
            .unwrap();

        let surface = instance
            .create_surface(window.clone())
            .expect("Create surface");
        let swapchain_format = TextureFormat::Bgra8UnormSrgb;
        let surface_config = SurfaceConfiguration {
            usage: TextureUsages::RENDER_ATTACHMENT,
            format: swapchain_format,
            width: physical_size.width,
            height: physical_size.height,
            present_mode: PresentMode::Fifo,
            alpha_mode: CompositeAlphaMode::Opaque,
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &surface_config);


        let font_data = std::fs::read(FONT_PATH).expect("Failed to load font");

        let mut font_db = Database::default();
        font_db.load_font_data(font_data);

        let mut font_system = FontSystem::new_with_locale_and_db("fr/fr".to_string(), font_db);

        let swash_cache = SwashCache::new();
        let cache = Cache::new(&device);
        let viewport = Viewport::new(&device, &cache);
        let mut atlas = TextAtlas::new(&device, &queue, &cache, swapchain_format);

        let multisample = MultisampleState {
            count: 4, // 4x MSAA
            mask: !0,
            alpha_to_coverage_enabled: false,
        };

        let text_renderer =
            TextRenderer::new(&mut atlas, &device, multisample.clone(), None);

        let mut texts = text::Texts::new();

        let text_content = format!("Nombre d'éléments: {}", problem.geometry.elements.elem_numbers.len());
        let text = text::Text::new("elem-number", &text_content, &mut font_system, &window, 0.05, 0.9, Color::rgb(255, 0, 0), true);

        texts.add_text(text);

        let text_content = "Press H to display help";
        let text = text::Text::new("help-hint", text_content, &mut font_system, &window, 0.05, 0.05, Color::rgb(255, 0, 0), true);

        texts.add_text(text);

        let text_content = "\n      L to toggle the grid\n      M to toggle the mesh\n      Q to quit\n      C to change the colormap\n";
        let text = text::Text::new("help", text_content, &mut font_system, &window, 0.05, 0.0525, Color::rgb(255, 0, 0), false);

        texts.add_text(text);

        let text_content = "<Placeholder>";
        let text = text::Text::new("colormap", text_content, &mut font_system, &window, 0.05, 0.85, Color::rgb(255, 0, 0), true);

        texts.add_text(text);

        let mesh_renderer = mesh_renderer::MeshRenderer::new(&device, &surface_config, vertices, indices, lines, multisample);

        Self {
            device,
            queue,
            surface,
            surface_config,
            font_system,
            swash_cache,
            viewport,
            atlas,
            text_renderer,
            texts,
            mesh_renderer,
            window,
            scale: 1.0,
            drag_state: DragState::new(),
            color_map: 0,
            multisample,
        }
    }
}

struct Application {
    window_state: Option<WindowState>,
    start_function: fn() -> FemProblem,
}

impl winit::application::ApplicationHandler for Application {
    fn resumed(&mut self, event_loop: &winit::event_loop::ActiveEventLoop) {
        if self.window_state.is_some() {
            return;
        }

        let problem = (self.start_function)();
        let (vertices, indices) = generate_colored_grid(&problem);
        // let (vertices, indices) = (vec![], vec![]);

        // Vertices lines
        let mut lines: Vec<(u32, u32)> = Vec::new();

        for i in 0..problem.geometry.elements.elem_numbers.len() as usize {
            let mut nodes = Vec::new();
            for j in 0..problem.geometry.elements.n_local_nodes as usize{
                nodes.push(problem.geometry.elements.elements[problem.geometry.elements.n_local_nodes as usize * i + j]);
            }
            for j in 0..problem.geometry.elements.n_local_nodes as usize {
                let node1 = nodes[j];
                let node2 = nodes[(j+1) % problem.geometry.elements.n_local_nodes as usize];
                if node1 < node2 {
                    lines.push((node2, node1));
                }
            }
        }

        for i in 0..problem.geometry.edges.elem_numbers.len() as usize {
            let mut nodes = Vec::new();
            for j in 0..problem.geometry.edges.n_local_nodes as usize{
                nodes.push(problem.geometry.edges.elements[problem.geometry.edges.n_local_nodes as usize * i + j]);
            }
            for j in 0..problem.geometry.edges.n_local_nodes as usize {
                let node1 = nodes[j];
                let node2 = nodes[(j+1) % problem.geometry.edges.n_local_nodes as usize];
                if node1 < node2 {
                    lines.push((node2, node1));
                }
            }
        }

        // Set up window
        let (width, height) = (800, 600);
        let window_attributes = Window::default_attributes()
            .with_inner_size(LogicalSize::new(width as f64, height as f64))
            .with_title("Rusted Fem");
        let window = Arc::new(event_loop.create_window(window_attributes).unwrap());

        self.window_state = Some(pollster::block_on(WindowState::new(
            window,
            vertices,
            indices,
            lines,
            problem,
        )));
    }

    fn window_event(
        &mut self,
        event_loop: &winit::event_loop::ActiveEventLoop,
        _window_id: winit::window::WindowId,
        event: WindowEvent,
    ) {
        let Some(state) = &mut self.window_state else {
            return;
        };

        let WindowState {
            window,
            device,
            queue,
            surface,
            surface_config,
            viewport,
            atlas,
            texts,
            text_renderer,
            font_system,
            swash_cache,
            mesh_renderer,
            scale,
            drag_state,
            multisample,
            ..
        } = state;

        match event {
            WindowEvent::Resized(size) => {
                surface_config.width = size.width;
                surface_config.height = size.height;
                surface.configure(&device, &surface_config);
                window.request_redraw();
            }
            WindowEvent::MouseWheel { device_id: _, delta, phase: _ } => {
                match delta {
                    winit::event::MouseScrollDelta::LineDelta(_, y) => {
                        *scale += y / 10.0;
                    }
                    winit::event::MouseScrollDelta::PixelDelta(pos) => {
                        *scale += pos.y as f32/ 10.0;
                    }
                }
                window.request_redraw();
            }
            WindowEvent::MouseInput { device_id: _, state, button } => {
                if button != MouseButton::Left {
                    return;
                }
                drag_state.is_dragging = state == winit::event::ElementState::Pressed;
            }
            WindowEvent::KeyboardInput { device_id: _, event, is_synthetic } => {
                if is_synthetic {
                    return;
                }
                
                if event.state == winit::event::ElementState::Pressed {
                    match event.logical_key.to_text().unwrap_or("".into()) {
                        "h" | "H" => {
                            if let Some(text) = texts.get_text("help") {
                                text.toggle_visibility();
                            }
                        }
                        "l" | "L" => {
                            mesh_renderer.toggle_lines();
                        }
                        "m" | "M" => {
                            mesh_renderer.toggle_mesh();
                        }
                        "q" | "Q" => {
                            event_loop.exit();
                        }
                        "c" | "C" => {
                            state.color_map = (state.color_map + 1) % ColorMap::list().len();
                        }
                        _ => {}
                    }
                    window.request_redraw();
                }
                
            }
            // on move, translate the coordinates
            WindowEvent::CursorMoved { device_id: _, position } => {
                
                let (x, y) = (position.x as f32, position.y as f32);
                let (width, height) = (surface_config.width as f32, surface_config.height as f32);
                // Convert to normalized device coordinates
                let new_cursor_pos = (
                    2.0 * x / width - 1.0,
                    1.0 - 2.0 * y / height,
                );

                if drag_state.is_dragging {
                    // Compute displacement from last position
                    let dx = new_cursor_pos.0 - drag_state.last_cursos_pos.0;
                    let dy = new_cursor_pos.1 - drag_state.last_cursos_pos.1;

                    // Accumulate translation
                    drag_state.translation.0 += dx;
                    drag_state.translation.1 += dy;
                }

                // Update last cursor position
                drag_state.last_cursos_pos = new_cursor_pos;
                window.request_redraw();
            }
            WindowEvent::RedrawRequested => {
                viewport.update(
                    &queue,
                    Resolution {
                        width: surface_config.width,
                        height: surface_config.height,
                    },
                );

                if let Some(text) = texts.get_text("colormap") {
                    text.set_text(&format!("ColorMap: {}", ColorMap::list()[state.color_map]), font_system);
                }

                texts.fonts_resize(font_system, surface_config.height as f32/30.0);

                let text_areas = texts.render(&window);

                text_renderer
                    .prepare(
                        device,
                        queue,
                        font_system,
                        atlas,
                        viewport,
                        text_areas,
                        swash_cache,
                    )
                    .unwrap();

                let msaa_texture = device.create_texture(&wgpu::TextureDescriptor {
                    label: Some("MSAA Texture"),
                    size: wgpu::Extent3d {
                        width: surface_config.width,
                        height: surface_config.height,
                        depth_or_array_layers: 1,
                    },
                    view_formats: &[],
                    mip_level_count: 1,
                    sample_count: multisample.count,
                    dimension: wgpu::TextureDimension::D2,
                    format: surface_config.format,
                    usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
                });
                
                let msaa_view = msaa_texture.create_view(&wgpu::TextureViewDescriptor::default());
                let frame = surface.get_current_texture().unwrap();
                let frame_view = frame.texture.create_view(&TextureViewDescriptor::default());
                let (msaa_view, view) = if multisample.count > 1 {
                    (&msaa_view, Some(&frame_view))
                } else {
                    (&frame_view, None)
                };
                    
                frame.texture.create_view(&TextureViewDescriptor::default());

                let mut encoder = device.create_command_encoder(&CommandEncoderDescriptor::default());

                let mut render_pass = encoder.begin_render_pass(&RenderPassDescriptor {
                    label: None,
                    color_attachments: &[Some(RenderPassColorAttachment {
                        view: &msaa_view,
                        resolve_target: view,
                        ops: Operations {
                            load: wgpu::LoadOp::Clear(wgpu::Color {
                                r: 255.0/255.0,
                                g: 255.0/255.0,
                                b: 181.0/255.0,
                                a: 1.0,
                            }),
                            store: wgpu::StoreOp::Discard,
                        },
                    })],
                    depth_stencil_attachment: None,
                    occlusion_query_set: None,
                    timestamp_writes: None,
                });

                let color_map = ColorMap::list()[state.color_map];

                mesh_renderer.render(
                    &mut render_pass,
                    &device,
                    &surface_config,
                    *scale,
                    drag_state.translation,
                    Some(color_map),
                );
                
                text_renderer
                    .render(&atlas, &state.viewport, &mut render_pass)
                    .unwrap();
            
                drop(render_pass);

                let command_buffer = encoder.finish();
                queue.submit(Some(command_buffer));
                frame.present();
                atlas.trim();
            }
            WindowEvent::CloseRequested => event_loop.exit(),
            _ => {}
        }
    }
}