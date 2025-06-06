#![allow(non_snake_case)]

use crate::Error;
use anyhow::{Error, Ok};
use std::ops::{Index, IndexMut};
use std::{f64::NAN, rc::Rc};

mod fem_io;

pub struct ProjectParameters {
    pub E_ground: f64,
    pub nu_ground: f64,
    pub rho_ground: f64,
    pub A_ground: f64,
    pub B_ground: f64,
    pub C_ground: f64,
    pub E_tyre: f64,
    pub nu_tyre: f64,
    pub rho_tyre: f64,
    pub A_tyre: f64,
    pub B_tyre: f64,
    pub C_tyre: f64,
    pub E_wheel: f64,
    pub nu_wheel: f64,
    pub rho_wheel: f64,
    pub A_wheel: f64,
    pub B_wheel: f64,
    pub C_wheel: f64,
}

pub struct FemNodes {
    pub numbers: Vec<u32>,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

pub struct FemMesh {
    pub n_local_nodes: u8,
    pub nodes: Rc<FemNodes>,
    pub elements: Vec<u32>,
    pub elem_numbers: Vec<u32>,
}

pub struct FemDomain {
    pub mesh: Rc<FemMesh>,
    pub elements: Vec<u32>,
    pub name: String,
}

#[derive(PartialEq)]
pub enum FemElementType {
    FemTriangle,
    FemQuad,
    FemLine,
}

#[derive(PartialEq, Clone, Copy)]
pub enum FemBoundaryType {
    DirichletX,
    DirichletY,
    DirichletXY,
    DirichletN,
    DirichletT,
    DirichletNT,
    NeumannX,
    NeumannY,
    NeumannN,
    NeumannT,
    Undefined = -1,
}

#[derive(PartialEq)]
pub enum FemElasticCase {
    PlanarStress,
    PlanarStrain,
    Axisymmetric,
    ProjectHybrid,
}

pub enum FemRenumType {
    No,
    XNum,
    YNum,
    RCMK,
}

pub enum FemSolverType {
    Full,
    Band,
    Front,
}

pub struct FemGeo {
    pub lx_ground: f64,
    pub ly_ground: f64,
    pub radius_tyre: f64,
    pub radius_wheel: f64,
    pub radius_hole: f64,
    pub tyre_initial_position: f64,

    pub geo_size: fn(f64, f64) -> f64,
    pub h: f64,
    pub element_type: FemElementType,
    pub nodes: Rc<FemNodes>,
    pub elements: Rc<FemMesh>,
    pub edges: Rc<FemMesh>,
    pub domains: Vec<Rc<FemDomain>>,
}

pub struct FemDiscrete {
    pub n: u32,
    pub element_type: FemElementType,
    pub x2: fn() -> ([f64; 4], [f64; 4]),
    pub phi2: fn(xsi: f64, eta: f64) -> [f64; 4],
    pub dphi2dx: fn(xsi: f64, eta: f64) -> ([f64; 4], [f64; 4]),
    pub x: fn() -> [f64; 2],
    pub phi: fn(xsi: f64) -> [f64; 2],
    pub dphidx: fn(xsi: f64) -> [f64; 2],
}

pub struct FemIntegration {
    pub n: u32,
    pub xsi: Vec<f64>,
    pub eta: Vec<f64>,
    pub w: Vec<f64>,
}

pub struct FemFullSystem {
    pub B: Vec<f64>,
    pub A: FemFullSystemStiffnessMatrix,
    pub size: usize,
}

pub struct FemFullSystemStiffnessMatrix {
    pub a: Vec<*mut f64>,
    pub size: usize,
}

pub struct FemBoundaryCondition {
    pub domain: Rc<FemDomain>,
    pub boundary_type: FemBoundaryType,
    pub value1: f64,
    pub value2: f64,
}

#[derive(Clone)]
pub struct FemConstrainedNode {
    pub boundary_type: FemBoundaryType,
    pub nx: f64,
    pub ny: f64,
    pub value1: f64,
    pub value2: f64,
}

pub struct FemProblem {
    pub E: f64,
    pub nu: f64,
    pub rho: f64,
    pub gx: f64,
    pub gy: f64,
    pub A: f64,
    pub B: f64,
    pub C: f64,
    pub material_parameters: ProjectParameters,
    pub planar_strain_stress: FemElasticCase,
    pub conditions: Vec<FemBoundaryCondition>,
    pub solution: Vec<f64>,
    pub residual: Vec<f64>,
    pub geometry: Rc<FemGeo>,
    pub space: FemDiscrete,
    pub rule: FemIntegration,
    pub space_edge: FemDiscrete,
    pub rule_edge: FemIntegration,
    pub system: FemFullSystem,
    pub constrained_nodes: Vec<FemConstrainedNode>,
}

impl FemProblem {
    pub fn elasticity_add_boundary_condition(
        &mut self,
        domain_name: &str,
        boundary_type: FemBoundaryType,
        value1: f64,
        mut value2: f64,
    ) {
        let domain = self
            .geometry
            .get_domain(domain_name)
            .expect("Domain not found");
        value2 = if (boundary_type != FemBoundaryType::DirichletXY)
            && (boundary_type != FemBoundaryType::DirichletNT)
        {
            NAN
        } else {
            value2
        };

        let dirichlet_types = vec![
            FemBoundaryType::DirichletX,
            FemBoundaryType::DirichletY,
            FemBoundaryType::DirichletXY,
            FemBoundaryType::DirichletNT,
        ];
        let is_dirichlet = dirichlet_types.contains(&boundary_type);

        if is_dirichlet {
            // Ensure that there is only one Dirichlet boundary condition per domain
            for condition in &self.conditions {
                if condition.domain.name != domain_name {
                    continue;
                }
                if dirichlet_types.contains(&condition.boundary_type) {
                    println!(
                        "There is already a Dirichlet boundary condition in domain {}",
                        domain_name
                    );
                    Error!("Only one Dirichlet boundary condition is allowed per domain");
                }
            }
            let constrained_node = FemConstrainedNode {
                boundary_type,
                nx: NAN,
                ny: NAN,
                value1,
                value2,
            };
            let elem = &domain.elements;
            let mesh = &domain.mesh;
            let n_elem = elem.len();
            if boundary_type == FemBoundaryType::DirichletX
                || boundary_type == FemBoundaryType::DirichletY
                || boundary_type == FemBoundaryType::DirichletXY
            {
                for i_elem in 0..n_elem {
                    for i in 0..2 {
                        let node = mesh.elements[2 * elem[i_elem] as usize + i] as usize;
                        self.constrained_nodes[node] = constrained_node.clone();
                    }
                }
            } else {
                // need to compute normals
                let mut normals_x = vec![0.0; self.geometry.nodes.numbers.len()];
                let mut normals_y = vec![0.0; self.geometry.nodes.numbers.len()];

                let mesh_elems = &mesh.elements;
                let nodes = &self.geometry.nodes;

                for i_elem in 0..n_elem as usize {
                    let node0 = mesh_elems[2 * elem[i_elem] as usize] as usize;
                    let node1 = mesh_elems[2 * elem[i_elem] as usize + 1] as usize;
                    let tx = nodes.x[node1] - nodes.x[node0];
                    let ty = nodes.y[node1] - nodes.y[node0];
                    let nx = ty;
                    let ny = -tx;

                    normals_x[node0] += nx;
                    normals_y[node0] += ny;
                    normals_x[node1] += nx;
                    normals_y[node1] += ny;
                }

                for i_elem in 0..n_elem as usize {
                    for i in 0..2 as usize {
                        let node = mesh_elems[2 * elem[i_elem] as usize + i] as usize;
                        let nx = normals_x[node];
                        let ny = normals_y[node];
                        let n = (nx * nx + ny * ny).sqrt();
                        self.constrained_nodes[node] = constrained_node.clone();
                        self.constrained_nodes[node].nx = nx / n;
                        self.constrained_nodes[node].ny = ny / n;
                    }
                }
            }
        }

        let boundary_condition = FemBoundaryCondition {
            domain,
            boundary_type,
            value1,
            value2,
        };
        self.conditions.push(boundary_condition);
    }
}

impl FemGeo {
    pub fn get_size_default(&self, _: f64, _: f64) -> f64 {
        return self.h;
    }

    pub fn set_size_callback(&mut self, f: fn(f64, f64) -> f64) {
        self.geo_size = f;
    }

    pub fn get_domain(&self, name: &str) -> Option<Rc<FemDomain>> {
        for domain in &self.domains {
            if domain.name == name {
                return Some(domain.clone());
            }
        }
        None
    }
}

impl FemDomain {
    pub fn setName(&mut self, name: &str) {
        self.name = name.to_string();
    }
}

impl FemIntegration {
    const GAUSS_QUAD4_XSI: [f64; 4] = [
        -0.577350269189626,
        -0.577350269189626,
        0.577350269189626,
        0.577350269189626,
    ];
    const GAUSS_QUAD4_ETA: [f64; 4] = [
        0.577350269189626,
        -0.577350269189626,
        -0.577350269189626,
        0.577350269189626,
    ];
    const GAUSS_QUAD4_WEIGHT: [f64; 4] = [
        1.000000000000000,
        1.000000000000000,
        1.000000000000000,
        1.00000000000000,
    ];
    const GAUSS_TRI3_XSI: [f64; 3] = [0.166666666666667, 0.666666666666667, 0.166666666666667];
    const GAUSS_TRI3_ETA: [f64; 3] = [0.166666666666667, 0.166666666666667, 0.666666666666667];
    const GAUSS_TRI3_WEIGHT: [f64; 3] = [0.166666666666667, 0.166666666666667, 0.166666666666667];
    const GAUSS_EDGE2_XSI: [f64; 2] = [0.577350269189626, -0.577350269189626];
    const GAUSS_EDGE2_WEIGHT: [f64; 2] = [1.000000000000000, 1.000000000000000];

    pub fn new(n: u32, element_type: FemElementType) -> Result<Self, Error> {
        match (n, element_type) {
            (4, FemElementType::FemQuad) => {
                let xsi = Vec::from(Self::GAUSS_QUAD4_XSI);
                let eta = Vec::from(Self::GAUSS_QUAD4_ETA);
                let w = Vec::from(Self::GAUSS_QUAD4_WEIGHT);
                Ok(Self { n, xsi, eta, w })
            }
            (3, FemElementType::FemTriangle) => {
                let xsi = Vec::from(Self::GAUSS_TRI3_XSI);
                let eta = Vec::from(Self::GAUSS_TRI3_ETA);
                let w = Vec::from(Self::GAUSS_TRI3_WEIGHT);
                Ok(Self { n, xsi, eta, w })
            }
            (2, FemElementType::FemLine) => {
                let xsi = Vec::from(Self::GAUSS_EDGE2_XSI);
                let eta = vec![0.0];
                let w = Vec::from(Self::GAUSS_EDGE2_WEIGHT);
                Ok(Self { n, xsi, eta, w })
            }
            _ => Err(Error::msg("Invalid number of integration points")),
        }
    }
}

impl FemDiscrete {
    const DUMMY_2: [f64; 2] = [0.0, 0.0];
    const DUMMY_4: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
    const DUMMY_42: ([f64; 4], [f64; 4]) = (Self::DUMMY_4, Self::DUMMY_4);

    pub fn new(n: i32, element_type: FemElementType) -> Result<Self, Error> {
        match (n, element_type) {
            (2, FemElementType::FemLine) => Ok(Self {
                n: 2,
                element_type: FemElementType::FemLine,
                x2: || Self::DUMMY_42,
                phi2: |_, _| Self::DUMMY_4,
                dphi2dx: |_, _| Self::DUMMY_42,
                x: || [1.0, -1.0],
                phi: |x| [0.5 * (1.0 - x), 0.5 * (1.0 + x)],
                dphidx: |_| [-0.5, 0.5],
            }),
            (3, FemElementType::FemTriangle) => Ok(Self {
                n: 3,
                element_type: FemElementType::FemTriangle,
                x2: || ([0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]),
                phi2: |xsi, eta| [1.0 - xsi - eta, xsi, eta, 0.0],
                dphi2dx: |_, _| ([-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0]),
                x: || Self::DUMMY_2,
                phi: |_| Self::DUMMY_2,
                dphidx: |_| Self::DUMMY_2,
            }),
            (4, FemElementType::FemQuad) => Ok(Self {
                n: 4,
                element_type: FemElementType::FemQuad,
                x2: || ([1.0, -1.0, -1.0, 1.0], [1.0, 1.0, -1.0, -1.0]),
                phi2: |xsi, eta| {
                    [
                        0.25 * (1.0 + xsi) * (1.0 + eta),
                        0.25 * (1.0 - xsi) * (1.0 + eta),
                        0.25 * (1.0 - xsi) * (1.0 - eta),
                        0.25 * (1.0 + xsi) * (1.0 - eta),
                    ]
                },
                dphi2dx: |xsi, eta| {
                    (
                        [
                            0.25 * (1.0 + eta),
                            -0.25 * (1.0 + eta),
                            -0.25 * (1.0 - eta),
                            0.25 * (1.0 - eta),
                        ],
                        [
                            0.25 * (1.0 + xsi),
                            0.25 * (1.0 - xsi),
                            -0.25 * (1.0 - xsi),
                            -0.25 * (1.0 + xsi),
                        ],
                    )
                },
                x: || Self::DUMMY_2,
                phi: |_| Self::DUMMY_2,
                dphidx: |_| Self::DUMMY_2,
            }),
            _ => Err(Error::msg("Invalid number of discrete points")),
        }
    }
}

impl FemFullSystem {
    pub fn new(size: usize) -> Self {
        let mut b = vec![0.0; size * (size + 1)];
        let mut a: Vec<*mut f64> = Vec::with_capacity(size);
        unsafe { a.push(b.as_mut_ptr().add(size)); }
        for i in 1..size {
            let ptr = unsafe { a[i - 1].add(size) };
            a.push(ptr);
        }

        Self {
            A: FemFullSystemStiffnessMatrix { a, size },
            B: b,
            size,
        }
    }

    pub fn clear(&mut self) {
        for i in 0..self.size {
            for j in 0..self.size {
                self.A[i][j] = 0.0;
            }
            self.B[i] = 0.0;
        }
    }

    pub fn gaussian_elimination(&mut self) -> &Vec<f64> {
        let A = &mut self.A;
        let B = &mut self.B;
        let size = self.size;

        /* Gaussian elimination */
        for i in 0..size {
            let pivot = A[i][i];
            if pivot.abs() < 1.0e-16 {
                println!("Pivot index: {}", i);
                println!("Pivot value: {}", pivot);
                Error!("Pivot is too small");
            }
            for j in i + 1..size {
                let factor = A[j][i] / pivot;
                for k in i..size {
                    A[j][k] -= factor * A[i][k];
                }
                B[j] -= factor * B[i];
            }
        }

        /* Backward substitution */
        for i in (0..size).rev() {
            let mut sum = 0.0;
            for j in i + 1..size {
                sum += A[i][j] * B[j];
            }
            B[i] = (B[i] - sum) / A[i][i];
        }

        B
    }

    pub fn band_elimination(&mut self, band: usize) -> &Vec<f64> {
        let A = &mut self.A;
        let B = &mut self.B;
        let size = self.size;

        /* Incomplete Cholesky factorization */

        for i in 0..size {
            let pivot = A[i][i];
            if pivot.abs() < 1.0e-16 {
                println!("Pivot index: {}", i);
                println!("Pivot value: {}", pivot);
                Error!("Pivot is too small");
            }
            let jend = size.min(i + band + 1); // copilot min i+band+1
            for j in i + 1..jend {
                let factor = A[j][i] / pivot;
                for k in j..jend {
                    A[j][k] -= factor * A[i][k];
                }
                B[j] -= factor * B[i];
            }
        }

        for i in (0..size).rev() {
            let mut sum = 0.0;
            for j in i + 1..size.min(i + band) {
                sum += A[i][j] * B[j];
            }
            B[i] = (B[i] - sum) / A[i][i];
        }

        B
    }

    pub fn constrain(&mut self, node: usize, value: f64) {
        let A = &mut self.A;
        let B = &mut self.B;
        let size = self.size;

        for i in 0..size {
            B[i] -= A[i][node] * value;
            A[i][node] = 0.0;
        }

        for i in 0..size {
            A[node][i] = 0.0;
        }

        A[node][node] = 1.0;
        B[node] = value;
    }
}

impl Index<usize> for FemFullSystemStiffnessMatrix {
    type Output = [f64];

    fn index(&self, i: usize) -> &Self::Output {
        let ptr = self.a[i] as *const f64;
        unsafe { std::slice::from_raw_parts(ptr, self.size + 1) }
    }
}

impl IndexMut<usize> for FemFullSystemStiffnessMatrix {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        let ptr = self.a[i];
        unsafe { std::slice::from_raw_parts_mut(ptr, self.size + 1) }
    }
}

impl ProjectParameters {
    pub fn default() -> Self {
        Self {
            A_ground: 0.0,
            A_tyre: 0.0,
            A_wheel: 0.0,
            E_ground: 0.0,
            E_tyre: 0.0,
            E_wheel: 0.0,
            nu_ground: 0.0,
            nu_tyre: 0.0,
            nu_wheel: 0.0,
            rho_ground: 0.0,
            rho_tyre: 0.0,
            rho_wheel: 0.0,
            B_ground: 0.0,
            B_tyre: 0.0,
            B_wheel: 0.0,
            C_ground: 0.0,
            C_tyre: 0.0,
            C_wheel: 0.0,
        }
    }
}

impl FemBoundaryCondition {
    pub fn is_neumann(&self) -> bool {
        match self.boundary_type {
            FemBoundaryType::NeumannX
            | FemBoundaryType::NeumannY
            | FemBoundaryType::NeumannN
            | FemBoundaryType::NeumannT => true,
            _ => false,
        }
    }
}