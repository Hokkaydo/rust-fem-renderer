use crate::Error;
use anyhow::{Error, Ok};
use std::{
    fs::File,
    io::{self, BufRead, BufReader, Write},
    rc::Rc,
};

use crate::fem::{
    FemBoundaryType, FemConstrainedNode, FemDiscrete, FemDomain, FemElasticCase, FemElementType,
    FemFullSystem, FemGeo, FemIntegration, FemMesh, FemNodes, FemProblem, FemRenumType,
    ProjectParameters,
};

impl FemProblem {
    pub fn new(
        geometry: Rc<FemGeo>,
        filename: &str,
        renum_type: FemRenumType,
    ) -> Result<Self, Error> {
        let file = File::open(filename)?;
        let reader = BufReader::new(file);

        // Initialize node numbering

        let n_nodes = geometry.nodes.x.len();

        // Renumber mesh nodes
        // fem_mesh_renumber(&geometry.elements, renum_type);

        // Initialize the problem structure
        let mut problem = Self {
            conditions: Vec::new(),
            solution: vec![0.0; 2 * n_nodes],
            residual: vec![0.0; 2 * n_nodes],
            constrained_nodes: vec![
                FemConstrainedNode {
                    boundary_type: FemBoundaryType::Undefined,
                    nx: f64::NAN,
                    ny: f64::NAN,
                    value1: f64::NAN,
                    value2: f64::NAN,
                };
                n_nodes
            ],
            material_parameters: ProjectParameters::default(),
            geometry: geometry.clone(),
            space: match geometry.elements.n_local_nodes {
                3 => FemDiscrete::new(3, FemElementType::FemTriangle)?,
                4 => FemDiscrete::new(4, FemElementType::FemQuad)?,
                _ => {
                    Error!("Unsupported element type");
                    panic!("Unsupported element type")
                }
            },
            rule: match geometry.elements.n_local_nodes {
                3 => FemIntegration::new(3, FemElementType::FemTriangle)?,
                4 => FemIntegration::new(4, FemElementType::FemQuad)?,
                _ => {
                    Error!("Unsupported element type");
                    panic!("Unsupported element type")
                }
            },
            space_edge: FemDiscrete::new(2, FemElementType::FemLine)?,
            rule_edge: FemIntegration::new(2, FemElementType::FemLine)?,
            system: FemFullSystem::new(2 * n_nodes),
            planar_strain_stress: FemElasticCase::PlanarStress,
            E: 0.0,
            nu: 0.0,
            rho: 0.0,
            gx: 0.0,
            gy: 0.0,
            A: 0.0,
            B: 0.0,
            C: 0.0,
        };

        // Parse file
        for line in reader.lines() {
            let line = line?;
            if line.starts_with("Type of problem") {
                if line.contains("Planar stresses") {
                    problem.planar_strain_stress = FemElasticCase::PlanarStress;
                } else if line.contains("Planar strains") {
                    problem.planar_strain_stress = FemElasticCase::PlanarStrain;
                } else if line.contains("Axi-symetric problem") {
                    problem.planar_strain_stress = FemElasticCase::Axisymmetric;
                } else if line.contains("Hybrid problem") {
                    problem.planar_strain_stress = FemElasticCase::ProjectHybrid;
                }
            } else if line.starts_with("Young modulus") {
                problem.E = Self::parse_value(&line)?;
            } else if line.starts_with("Poisson ratio") {
                problem.nu = Self::parse_value(&line)?;
            } else if line.starts_with("Mass density") {
                problem.rho = Self::parse_value(&line)?;
            } else if line.starts_with("Ground Young modul.") {
                problem.material_parameters.E_ground = Self::parse_value(&line)?;
            } else if line.starts_with("Ground Poisson rat.") {
                problem.material_parameters.nu_ground = Self::parse_value(&line)?;
            } else if line.starts_with("Ground Mass density") {
                problem.material_parameters.rho_ground = Self::parse_value(&line)?;
            } else if line.starts_with("Tyre Young modulus") {
                problem.material_parameters.E_tyre = Self::parse_value(&line)?;
            } else if line.starts_with("Tyre Poisson ratio") {
                problem.material_parameters.nu_tyre = Self::parse_value(&line)?;
            } else if line.starts_with("Tyre Mass density") {
                problem.material_parameters.rho_tyre = Self::parse_value(&line)?;
            } else if line.starts_with("Wheel Young modulus") {
                problem.material_parameters.E_wheel = Self::parse_value(&line)?;
            } else if line.starts_with("Wheel Poisson ratio") {
                problem.material_parameters.nu_wheel = Self::parse_value(&line)?;
            } else if line.starts_with("Wheel Mass density") {
                problem.material_parameters.rho_wheel = Self::parse_value(&line)?;
            } else if line.starts_with("Gravity-X") {
                problem.gx = Self::parse_value(&line)?;
            } else if line.starts_with("Gravity-Y") {
                problem.gy = Self::parse_value(&line)?;
            } else if line.starts_with("Boundary condition") {
                let tokens: Vec<&str> = line.split_whitespace().collect();
                if tokens.len() >= 5 {
                    let boundary_type = match tokens[3] {
                        "Dirichlet-X" => FemBoundaryType::DirichletX,
                        "Dirichlet-Y" => FemBoundaryType::DirichletY,
                        "Neumann-X" => FemBoundaryType::NeumannX,
                        "Neumann-Y" => FemBoundaryType::NeumannY,
                        _ => FemBoundaryType::Undefined,
                    };
                    let value1: f64 = tokens[5].replace(",", "").parse().unwrap_or(0.0);
                    let value2: f64 = tokens[6].replace(":", "").parse().unwrap_or(0.0);
                    let domain = tokens.last().unwrap();
                    Self::elasticity_add_boundary_condition(
                        &mut problem,
                        &domain,
                        boundary_type,
                        value1,
                        value2,
                    );
                }
            }
        }

        // Compute material properties

        let E = problem.E;
        let nu = problem.nu;
        let E_ground = problem.material_parameters.E_ground;
        let nu_ground = problem.material_parameters.nu_ground;
        let E_tyre = problem.material_parameters.E_tyre;
        let nu_tyre = problem.material_parameters.nu_tyre;
        let E_wheel = problem.material_parameters.E_wheel;
        let nu_wheel = problem.material_parameters.nu_wheel;

        match problem.planar_strain_stress {
            FemElasticCase::PlanarStress => {
                problem.A = E / (1.0 - nu * nu);
                problem.B = E * nu / (1.0 - nu * nu);
                problem.C = E / (2.0 * (1.0 + nu));
            }
            FemElasticCase::PlanarStrain | FemElasticCase::Axisymmetric => {
                problem.A = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
                problem.B = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
                problem.C = E / (2.0 * (1.0 + nu));
            }
            FemElasticCase::ProjectHybrid => {
                let material = &mut problem.material_parameters;
                material.A_ground =
                    E_ground * (1.0 - nu_ground) / ((1.0 + nu_ground) * (1.0 - 2.0 * nu_ground));
                material.B_ground =
                    E_ground * nu_ground / ((1.0 + nu_ground) * (1.0 - 2.0 * nu_ground));
                material.C_ground = E_ground / (2.0 * (1.0 + nu_ground));

                material.A_tyre = E_tyre / (1.0 - nu_tyre * nu_tyre);
                material.B_tyre = E_tyre * nu_tyre / (1.0 - nu_tyre * nu_tyre);
                material.C_tyre = E_tyre / (2.0 * (1.0 + nu_tyre));

                material.A_wheel = E_wheel / (1.0 - nu_wheel * nu_wheel);
                material.B_wheel = E_wheel * nu_wheel / (1.0 - nu_wheel * nu_wheel);
                material.C_wheel = E_wheel / (2.0 * (1.0 + nu_wheel));
            }
        }

        Ok(problem)
    }

    fn parse_value(line: &str) -> io::Result<f64> {
        line.split(':')
            .nth(1)
            .and_then(|s| s.trim().parse().ok())
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Invalid numerical value"))
    }

    pub fn print(&self) {
        print!("\n\n ======================================================================================= \n\n");
        println!(" Linear Elasticity Problem");
        if self.planar_strain_stress == FemElasticCase::ProjectHybrid {
            println!(
                "  Ground Young Modulus    E_ground    = {:14.7e}[N/m2]",
                self.material_parameters.E_ground
            );
            println!(
                "  Ground Poisson's Ratio  nu_ground   = {:14.7e}",
                self.material_parameters.nu_ground
            );
            println!(
                "  Ground Density          rho_ground  = {:14.7e}[kg/m3]",
                self.material_parameters.rho_ground
            );
            println!(
                "  Tyre Young Modulus      E_tyre      = {:14.7e}[N/m2]",
                self.material_parameters.E_tyre
            );
            println!(
                "  Tyre Poisson's Ratio    nu_tyre     = {:14.7e}",
                self.material_parameters.nu_tyre
            );
            println!(
                "  Tyre Density            rho_tyre    = {:14.7e}[kg/m3]",
                self.material_parameters.rho_tyre
            );
            println!(
                "  Wheel Young Modulus     E_wheel     = {:14.7e}[N/m2]",
                self.material_parameters.E_wheel
            );
            println!(
                "  Wheel Poisson's Ratio   nu_wheel    = {:14.7e}",
                self.material_parameters.nu_wheel
            );
            println!(
                "  Wheel Density           rho_wheel   = {:14.7e}[kg/m3]",
                self.material_parameters.rho_wheel
            );
        } else {
            println!("  Young Modulus       E   = {:14.7e}[N/m2]", self.E);
            println!("  Poisson's Ratio     nu  = {:14.7e}", self.nu);
            println!("  Density             rho = {:14.7e}[kg/m3]", self.rho);
        }
        println!("  Gravity-X           gx = {:14.7e}[m/s2]", self.gx);
        println!("  Gravity-Y           gy = {:14.7e}[m/s2]", self.gy);

        if self.planar_strain_stress == FemElasticCase::PlanarStrain {
            println!("  Planar strains formulation");
        }
        if self.planar_strain_stress == FemElasticCase::PlanarStress {
            println!("  Planar stress formulation");
        }
        if self.planar_strain_stress == FemElasticCase::Axisymmetric {
            println!("  Axisymmetric formulation");
        }
        if self.planar_strain_stress == FemElasticCase::ProjectHybrid {
            println!("  Project Hybrid formulation");
        }
        println!("  Boundary conditions");
        for condition in &self.conditions {
            print!("  {:20} :", condition.domain.name);
            if condition.boundary_type == FemBoundaryType::DirichletX {
                println!(
                    " imposing {:9.e} as the horizontal displacement",
                    condition.value1
                );
            }
            if condition.boundary_type == FemBoundaryType::DirichletY {
                println!(
                    " imposing {:9.e} as the vertical displacement",
                    condition.value1
                );
            }
            if condition.boundary_type == FemBoundaryType::DirichletXY {
                println!(
                    " imposing {:9.e} as the horizontal displacement",
                    condition.value1
                );
                println!(
                    " imposing {:9.e} as the vertical displacement",
                    condition.value2
                );
            }
            if condition.boundary_type == FemBoundaryType::DirichletN {
                println!(
                    " imposing {:9.e} as the normal displacement",
                    condition.value1
                );
            }
            if condition.boundary_type == FemBoundaryType::DirichletT {
                println!(
                    " imposing {:9.e} as the tangential displacement",
                    condition.value1
                );
            }
            if condition.boundary_type == FemBoundaryType::DirichletNT {
                println!(
                    " imposing {:9.e} as the normal displacement",
                    condition.value1
                );
                println!(
                    " imposing {:9.e} as the tangential displacement",
                    condition.value2
                );
            }
            if condition.boundary_type == FemBoundaryType::NeumannX {
                println!(
                    " imposing {:9.e} as the horizontal traction",
                    condition.value1
                );
            }
            if condition.boundary_type == FemBoundaryType::NeumannY {
                println!(
                    " imposing {:9.e} as the vertical traction",
                    condition.value1
                );
            }
            if condition.boundary_type == FemBoundaryType::NeumannN {
                println!(" imposing {:9.e} as the normal traction", condition.value1);
            }
            if condition.boundary_type == FemBoundaryType::NeumannT {
                println!(
                    " imposing {:9.e} as the tangential traction",
                    condition.value1
                );
            }
        }
        print!(" ======================================================================================= \n\n");
    }

    pub fn write(&self, filename: &str) -> Result<(), Error> {
        let mut file = File::create(filename)?;

        write!(file, "Type of problem    :  ")?;
        match self.planar_strain_stress {
            FemElasticCase::PlanarStress => writeln!(file, "Planar stress"),
            FemElasticCase::PlanarStrain => writeln!(file, "Planar strains"),
            FemElasticCase::Axisymmetric => writeln!(file, "Axi-symmetric problem"),
            FemElasticCase::ProjectHybrid => writeln!(file, "Hybrid problem"),
        }?;

        if self.planar_strain_stress == FemElasticCase::ProjectHybrid {
            writeln!(
                file,
                "Ground Young modul.: {:14.7e}",
                self.material_parameters.E_ground
            )?;
            writeln!(
                file,
                "Ground Poisson rat.: {:14.7e}",
                self.material_parameters.nu_ground
            )?;
            writeln!(
                file,
                "Ground Mass density: {:14.7e}",
                self.material_parameters.rho_ground
            )?;
            writeln!(
                file,
                "Tyre Young modulus : {:14.7e}",
                self.material_parameters.E_tyre
            )?;
            writeln!(
                file,
                "Tyre Poisson ratio : {:14.7e}",
                self.material_parameters.nu_tyre
            )?;
            writeln!(
                file,
                "Tyre Mass density  : {:14.7e}",
                self.material_parameters.rho_tyre
            )?;
            writeln!(
                file,
                "Wheel Young modulus: {:14.7e}",
                self.material_parameters.E_wheel
            )?;
            writeln!(
                file,
                "Wheel Poisson ratio: {:14.7e}",
                self.material_parameters.nu_wheel
            )?;
            writeln!(
                file,
                "Wheel Mass density : {:14.7e}",
                self.material_parameters.rho_wheel
            )?;
        } else {
            writeln!(file, "Young modulus      : {:14.7e}", self.E)?;
            writeln!(file, "Poisson ratio      : {:14.7e}", self.nu)?;
            writeln!(file, "Mass density       : {:14.7e}", self.rho)?;
        }

        writeln!(file, "Gravity-X          : {:14.7e}", self.gx)?;
        writeln!(file, "Gravity-Y          : {:14.7e}", self.gy)?;

        for condition in &self.conditions {
            write!(file, "Boundary condition : ")?;
            match condition.boundary_type {
                FemBoundaryType::DirichletX => write!(
                    file,
                    " Dirichlet-X        = {:14.7e}, {:14.7e}",
                    condition.value1, 0.0
                ),
                FemBoundaryType::DirichletY => write!(
                    file,
                    " Dirichlet-Y        = {:14.7e}, {:14.7e}",
                    condition.value1, 0.0
                ),
                FemBoundaryType::DirichletXY => write!(
                    file,
                    " Dirichlet-XY       = {:14.7e}, {:14.7e}",
                    condition.value1, condition.value2
                ),
                FemBoundaryType::DirichletN => write!(
                    file,
                    " Dirichlet-N        = {:14.7e}, {:14.7e}",
                    condition.value1, 0.0
                ),
                FemBoundaryType::DirichletT => write!(
                    file,
                    " Dirichlet-T        = {:14.7e}, {:14.7e}",
                    condition.value1, 0.0
                ),
                FemBoundaryType::DirichletNT => write!(
                    file,
                    " Dirichlet-NT       = {:14.7e}, {:14.7e}",
                    condition.value1, condition.value2
                ),
                FemBoundaryType::NeumannX => write!(
                    file,
                    " Neumann-X        = {:14.7e}, {:14.7e}",
                    condition.value1, 0.0
                ),
                FemBoundaryType::NeumannY => write!(
                    file,
                    " Neumann-Y        = {:14.7e}, {:14.7e}",
                    condition.value1, 0.0
                ),
                FemBoundaryType::NeumannN => write!(
                    file,
                    " Neumann-N        = {:14.7e}, {:14.7e}",
                    condition.value1, 0.0
                ),
                FemBoundaryType::NeumannT => write!(
                    file,
                    " Neumann-T        = {:14.7e}, {:14.7e}",
                    condition.value1, 0.0
                ),
                FemBoundaryType::Undefined => {
                    write!(file, " Undefined        = {:14.7e}, {:14.7e}", 0.0, 0.0)
                }
            }?;
            writeln!(file, ": {}", condition.domain.name)?;
        }
        Ok(())
    }

    pub fn write_solution(&self, filename: &str, n_fields_opt: Option<usize>) -> Result<(), Error> {
        let mut file = File::create(filename)?;
        let n_nodes = self.geometry.nodes.x.len();
        let n_fields = n_fields_opt.unwrap_or(2); // 2 fields by default, only value used in original fem.c
        writeln!(file, "Size {},{}", n_nodes, n_fields)?;
        for i in 0..n_nodes as usize {
            for j in 0..n_fields - 1 {
                write!(file, "{:.18e},", self.solution[n_fields * i + j])?;
            }
            write!(file, "{:18e}", self.solution[n_fields * i + n_fields - 1])?;
            writeln!(file)?;
        }
        Ok(())
    }

    pub fn read_solution(&self, filename: &str) -> Result<Vec<f64>, Error> {
        let file = File::open(filename)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        let size_line = lines.next().ok_or(io::Error::new(
            io::ErrorKind::Other,
            "Failed to read size line",
        ))??;
        let mut parts = size_line.split(' ').last().unwrap().split(',');
        let n_nodes: usize = parts.next().unwrap().parse().unwrap();
        let n_fields: usize = parts.next().unwrap().parse().unwrap();

        let mut solution = vec![0.0; n_nodes * n_fields];
        for i in 0..n_nodes {
            let line = lines.next().ok_or(io::Error::new(
                io::ErrorKind::Other,
                "Failed to read solution line",
            ))??;
            let mut parts = line.split(",").map(|x| x.parse::<f64>().unwrap());
            for j in 0..n_fields {
                solution[n_fields * i + j] = parts.next().expect("Failed to read solution value");
            }
        }
        println!(
            "Read solution with {} nodes and {} fields",
            n_nodes, n_fields
        );
        Ok(solution)
    }
}

impl FemGeo {
    pub fn new(filename: &str) -> Result<Self, Error> {
        let mut geometry = Self {
            lx_ground: 0.0,
            ly_ground: 0.0,
            radius_tyre: 0.0,
            radius_wheel: 0.0,
            radius_hole: 0.0,
            tyre_initial_position: 0.0,
            geo_size: |_, _| 0.0,
            h: 0.0,
            element_type: FemElementType::FemTriangle,
            nodes: Rc::new(FemNodes {
                numbers: Vec::new(),
                x: Vec::new(),
                y: Vec::new(),
            }),
            elements: Rc::new(FemMesh {
                n_local_nodes: 0,
                nodes: Rc::new(FemNodes {
                    numbers: Vec::new(),
                    x: Vec::new(),
                    y: Vec::new(),
                }),
                elements: Vec::new(),
                elem_numbers: Vec::new(),
            }),
            edges: Rc::new(FemMesh {
                n_local_nodes: 0,
                elem_numbers: Vec::new(),
                nodes: Rc::new(FemNodes {
                    numbers: Vec::new(),
                    x: Vec::new(),
                    y: Vec::new(),
                }),
                elements: Vec::new(),
            }),
            domains: Vec::new(),
        };
        let file = File::open(filename)?;
        let reader = io::BufReader::new(file);
        let mut lines = reader.lines();

        // Read number of nodes
        let nodes_line = lines
            .next()
            .ok_or(io::Error::new(io::ErrorKind::Other, "Failed to read line"))??;
        let n_nodes: usize = nodes_line
            .split_whitespace()
            .last()
            .unwrap()
            .parse()
            .unwrap();
        let mut the_nodes = FemNodes {
            numbers: (0..2*n_nodes as u32).collect(),
            x: vec![0.0; n_nodes],
            y: vec![0.0; n_nodes],
        };

        for i in 0..n_nodes {
            let node_line = lines.next().ok_or(io::Error::new(
                io::ErrorKind::Other,
                "Failed to read node line",
            ))??;
            let mut parts = node_line
                .split_whitespace()
                .skip(2)
                .map(|x| x.parse::<f64>().unwrap());
            the_nodes.x[i] = parts.next().expect("Failed to read node's x");
            the_nodes.y[i] = parts.next().expect("Failed to read node's y");
        }

        geometry.nodes = Rc::new(the_nodes);

        // Read number of edges
        let edges_line = lines.next().ok_or(io::Error::new(
            io::ErrorKind::Other,
            "Failed to read edges line",
        ))??;
        let n_elem: usize = edges_line
            .split_whitespace()
            .last()
            .unwrap()
            .parse()
            .unwrap();
        let mut the_edges = FemMesh {
            n_local_nodes: 2,
            nodes: geometry.nodes.clone(),
            elements: vec![0; 2 * n_elem],
            elem_numbers: (0..n_elem as u32).collect(),
        };

        for i in 0..n_elem {
            let edge_line = lines.next().ok_or(io::Error::new(
                io::ErrorKind::Other,
                "Failed to read edge line",
            ))??;
            let mut parts = edge_line
                .split_whitespace()
                .skip(2)
                .map(|x| x.parse::<u32>().unwrap());
            the_edges.elements[2 * i] = parts.next().expect("Failed to read edge's first node");
            the_edges.elements[2 * i + 1] =
                parts.next().expect("Failed to read edge's second node");
        }
        geometry.edges = Rc::new(the_edges);

        // Read number of elements
        let element_type_line = lines.next().ok_or(io::Error::new(
            io::ErrorKind::Other,
            "Failed to read element type line",
        ))??;
        let mut parts = element_type_line.split_whitespace().skip(2);
        let _element_type = parts.next().unwrap();
        let n_elem: usize = parts.next().unwrap().parse().unwrap();

        let mut elements = FemMesh {
            n_local_nodes: 0,
            nodes: geometry.nodes.clone(),
            elements: vec![0; 0],
            elem_numbers: (0..n_elem as u32).collect(),
        };
        // Read element type (triangles or quads)
        if _element_type.to_lowercase() == "triangles" {
            elements.n_local_nodes = 3;
            elements.elements = vec![0; 3 * n_elem];
            for i in 0..n_elem {
                let element_line = lines.next().ok_or(io::Error::new(
                    io::ErrorKind::Other,
                    "Failed to read triangle element line",
                ))??;

                let mut parts = element_line
                    .split_whitespace()
                    .skip(2)
                    .map(|x| x.parse::<u32>().unwrap());
                elements.elements[3 * i] =
                    parts.next().expect("Failed to read element's first node");
                elements.elements[3 * i + 1] =
                    parts.next().expect("Failed to read element's second node");
                elements.elements[3 * i + 2] =
                    parts.next().expect("Failed to read element's third node");
            }
        }

        if _element_type.to_lowercase() == "quads" {
            elements.n_local_nodes = 4;
            elements.elements = vec![0; 4 * n_elem];
            for i in 0..n_elem {
                let element_line = lines.next().ok_or(io::Error::new(
                    io::ErrorKind::Other,
                    "Failed to read quad element line",
                ))??;
                let mut parts = element_line
                    .split_whitespace()
                    .skip(2)
                    .map(|x| x.parse::<u32>().unwrap());
                elements.elements[4 * i] =
                    parts.next().expect("Failed to read element's first node");
                elements.elements[4 * i + 1] =
                    parts.next().expect("Failed to read element's second node");
                elements.elements[4 * i + 2] =
                    parts.next().expect("Failed to read element's third node");
                elements.elements[4 * i + 3] =
                    parts.next().expect("Failed to read element's fourth node");
            }
        }

        geometry.elements = Rc::new(elements);

        // Read domains
        let domains_line = lines.next().ok_or(io::Error::new(
            io::ErrorKind::Other,
            "Failed to read domains line",
        ))??;
        let n_domains: usize = domains_line
            .split_whitespace()
            .last()
            .expect("Failed to read number of domains")
            .parse()
            .expect("Failed to parse number of domains");

        let mut _trash = lines.next().unwrap()?; // Skip first domain id
       
        for _i_domain in 0..n_domains {
            let mut domain = FemDomain {
                name: String::new(),
                elements: Vec::new(),
                mesh: geometry.edges.clone(),
            };
            let name_line = lines.next().ok_or(io::Error::new(
                io::ErrorKind::Other,
                "Failed to read domain name line",
            ))??;
            
            domain.name = name_line.split_whitespace().last().unwrap().to_string();
            lines.next(); // Skip number of elements
            domain.elements = Vec::new();

            while let Some(line) = lines.next() {
                let elem_line: String = line?;
                if elem_line.contains("Domain") {
                    // next domain's id found -> skip it
                    break;
                }
                let parts = elem_line.split_whitespace();
                for part in parts {
                    domain.elements.push(part.parse().unwrap());
                }
            }

            geometry.domains.push(Rc::new(domain));
        }

        Ok(geometry)
    }

    pub fn print(&self) {
        let nodes = &self.nodes;
        println!("Number of nodes: {}", nodes.x.len());
        for i in 0..nodes.x.len() {
            println!("{:6} : {:14.7e} {:14.7e}", i, nodes.x[i], nodes.y[i]);
        }

        let edges = &self.edges;
        println!("Number of edges: {}", edges.elem_numbers.len());
        for i in 0..edges.elem_numbers.len() {
            println!(
                "{:6} : {:6} {:6}",
                i,
                edges.elements[2 * i],
                edges.elements[2 * i + 1]
            );
        }

        let elements = &self.elements;
        if elements.n_local_nodes == 3 {
            println!("Number of triangles: {}", elements.elements.len() / 3);
            for i in 0..elements.elements.len() / 3 {
                println!(
                    "{:6} : {:6} {:6} {:6}",
                    i,
                    elements.elements[3 * i],
                    elements.elements[3 * i + 1],
                    elements.elements[3 * i + 2]
                );
            }
        }
        if elements.n_local_nodes == 4 {
            println!("Number of quads: {}", elements.elements.len() / 4);
            for i in 0..elements.elements.len() / 4 {
                println!(
                    "{:6} : {:6} {:6} {:6} {:6}",
                    i,
                    elements.elements[4 * i],
                    elements.elements[4 * i + 1],
                    elements.elements[4 * i + 2],
                    elements.elements[4 * i + 3]
                );
            }
        }

        let domains = &self.domains;
        println!("Number of domains: {}", domains.len());
        for i in 0..domains.len() {
            println!("  Domain: {}", i);
            println!("  Name: {}", domains[i].name);
            println!("  Number of elements: {}", domains[i].elements.len());
            for j in 0..domains[i].elements.len() {
                print!("{:6} ", domains[i].elements[j]);
                if ((j + 1) != domains[i].elements.len()) && ((j + 1) % 10 == 0) {
                    println!();
                }
            }
            println!();
        }
    }

    pub fn write(&self, filename: &str) -> Result<(), Error> {
        let mut file = File::create(filename)?;
        let nodes = &self.nodes;

        writeln!(file, "Number of nodes {}", nodes.x.len())?;
        for i in 0..nodes.x.len() {
            writeln!(file, "{:6} : {:14.7e} {:14.7e}", i, nodes.x[i], nodes.y[i])?;
        }

        let edges = &self.edges;
        writeln!(file, "Number of edges {}", edges.elem_numbers.len())?;
        for i in 0..edges.elem_numbers.len() {
            writeln!(
                file,
                "{:6} : {:6} {:6}",
                i,
                edges.elements[2 * i],
                edges.elements[2 * i + 1]
            )?;
        }

        let elements = &self.elements;
        if elements.n_local_nodes == 3 {
            writeln!(file, "Number of triangles {}", elements.elements.len() / 3)?;
            for i in 0..elements.elements.len() / 3 {
                writeln!(
                    file,
                    "{:6} : {:6} {:6} {:6}",
                    i,
                    elements.elements[3 * i],
                    elements.elements[3 * i + 1],
                    elements.elements[3 * i + 2]
                )?;
            }
        }
        if elements.n_local_nodes == 4 {
            writeln!(file, "Number of quads {}", elements.elements.len() / 4)?;
            for i in 0..elements.elements.len() / 4 {
                writeln!(
                    file,
                    "{:6} : {:6} {:6} {:6} {:6}",
                    i,
                    elements.elements[4 * i],
                    elements.elements[4 * i + 1],
                    elements.elements[4 * i + 2],
                    elements.elements[4 * i + 3]
                )?;
            }
        }

        let domains = &self.domains;
        writeln!(file, "Number of domains {}", domains.len())?;
        for i in 0..domains.len() {
            writeln!(file, "  Domain :      {}", i)?;
            writeln!(file, "  Name : {}", domains[i].name)?;
            writeln!(
                file,
                "  Number of elements :     {}",
                domains[i].elements.len()
            )?;
            for j in 0..domains[i].elements.len() {
                write!(file, "{:6}", domains[i].elements[j])?;
                if ((j + 1) != domains[i].elements.len()) && ((j + 1) % 10 == 0) {
                    writeln!(file)?;
                }
            }
            writeln!(file)?;
        }

        Ok(())
    }
}

impl FemDiscrete {
    pub fn print(&self) {
        let (xsi, eta) = (self.x2)();

        for i in 0..self.n {
            let phi = (self.phi2)(xsi[i as usize], eta[i as usize]);
            let (dphidxsi, dphideta) = (self.dphi2dx)(xsi[i as usize], eta[i as usize]);
            for j in 0..self.n {
                println!(
                    "{:6} :
                    (xsi = {:14.7e}, eta = {:14.7e}),
                     phi({}) = {:14.7e}, \t dphidxsi({}) = {:14.7e}, dphideta({}) = {:14.7e}",
                    i,
                    xsi[i as usize],
                    eta[i as usize],
                    j,
                    phi[j as usize],
                    j,
                    dphidxsi[j as usize],
                    j,
                    dphideta[j as usize]
                );
            }
            println!();
        }
    }
}

impl FemFullSystem {
    pub fn print(&self) {
        println!("Size: {}", self.size);
        println!("B:");
        for i in 0..self.size {
            print!("{:14.7e} ", self.B[i]);
            if (i + 1) % 10 == 0 {
                println!();
            }
        }
        println!();
        println!("A:");
        for i in 0..self.size {
            for j in 0..self.size {
                print!("{:14.7e} ", self.A[i][j]);
                if (j + 1) % 10 == 0 {
                    println!();
                }
            }
            println!();
        }
    }

    pub fn write_matrix(&self, filename: &str) -> Result<(), Error> {
        let mut file = File::create(filename)?;
        writeln!(file, "Size: {}", self.size)?;
        for i in 0..self.size {
            for j in 0..self.size {
                write!(file, "{:1.18e},", self.A[i][j])?;
            }
            writeln!(file)?;
        }
        Ok(())
    }

    pub fn write_vector(&self, filename: &str) -> Result<(), Error> {
        let mut file = File::create(filename)?;
        writeln!(file, "Size: {}", self.size)?;
        for i in 0..self.size {
            writeln!(file, "{:1.18e}", self.B[i])?;
        }
        
        Ok(())
    }
}
