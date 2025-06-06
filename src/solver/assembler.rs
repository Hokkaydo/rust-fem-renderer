use anyhow::{Error, Ok};

use crate::{fem::{FemBoundaryType, FemElasticCase, FemElementType, FemProblem}, program, Warning};



impl FemProblem {

    fn get_parameters_from_position(&self, x_mean: f64, y_mean: f64) -> (f64, f64, f64, f64) {
        let mut a = 0.0;
        let mut b = 0.0;
        let mut c = 0.0;
        let mut rho = 0.0;
        let the_problem = self;
        let a_ground = the_problem.material_parameters.A_ground;
        let b_ground = the_problem.material_parameters.B_ground;
        let c_ground = the_problem.material_parameters.C_ground;
        let a_tyre = the_problem.material_parameters.A_tyre;
        let b_tyre = the_problem.material_parameters.B_tyre;
        let c_tyre = the_problem.material_parameters.C_tyre;
        let a_wheel = the_problem.material_parameters.A_wheel;
        let b_wheel = the_problem.material_parameters.B_wheel;
        let c_wheel = the_problem.material_parameters.C_wheel;
        let rho_ground = the_problem.material_parameters.rho_ground;
        let rho_tyre = the_problem.material_parameters.rho_tyre;
        let rho_wheel = the_problem.material_parameters.rho_wheel;
        let lx_ground = the_problem.geometry.lx_ground;
        let ly_ground = the_problem.geometry.ly_ground;
        let radius_tyre = the_problem.geometry.radius_tyre;
        let radius_wheel = the_problem.geometry.radius_wheel;
        let tyre_initial_position = the_problem.geometry.tyre_initial_position;
        let d = ((x_mean - lx_ground / 2.0) * (x_mean - lx_ground / 2.0) + (y_mean - tyre_initial_position) * (y_mean - tyre_initial_position)).sqrt();
        if d <= radius_tyre && d > radius_wheel {
            a = a_tyre;
            b = b_tyre;
            c = c_tyre;
            rho = rho_tyre;
        } else if d <= radius_wheel {
            a = a_wheel;
            b = b_wheel;
            c = c_wheel;
            rho = rho_wheel;
        } else if y_mean <= ly_ground {
            a = a_ground;
            b = b_ground;
            c = c_ground;
            rho = rho_ground;
        } else {
            println!("Point outside the geometry, x_mean = {:9.2e}, y_mean = {:9.2e}", x_mean, y_mean);
        }
        (a, b, c, rho)
    }

    pub fn assemble_elements(&mut self) {
        let gx: f64 = self.gx;
        let gy: f64 = self.gy;

        let mesh = &self.geometry.elements;
        let nodes = &self.geometry.nodes;
        let geometry = &self.geometry;

        let mut a = self.A;
        let mut b = self.B;
        let mut c = self.C;
        let mut rho = self.rho;

        let n_local = mesh.n_local_nodes as usize;
        
        let mut map: [usize; 4] = [0; 4];
        let mut x = [0.0; 4];
        let mut y = [0.0; 4];
        let mut map_x: [usize; 4] = [0; 4];
        let mut map_y : [usize; 4]= [0; 4];

        for i_elem in mesh.elem_numbers.iter() {
            let i_elem = *i_elem as usize;
            for i in 0..n_local as usize {
                map[i] = mesh.elements[i_elem * n_local + i] as usize;
                x[i] = nodes.x[map[i] as usize];
                y[i] = nodes.y[map[i] as usize];
                map_x[i] = mesh.nodes.numbers[2*map[i]] as usize;
                map_y[i] = mesh.nodes.numbers[2*map[i] + 1] as usize;
            }
            let (x_mean, y_mean) = geometry.element_type.mean(x.to_vec(), y.to_vec()).expect("Should be able to get mean");
            
            if self.planar_strain_stress == FemElasticCase::ProjectHybrid {
                (a, b, c, rho) = self.get_parameters_from_position(x_mean, y_mean);
            }

            for i_integ in 0..self.rule.n as usize {
                let xsi = self.rule.xsi[i_integ];
                let eta = self.rule.eta[i_integ];
                let w = self.rule.w[i_integ];

                let phi = (self.space.phi2)(xsi, eta);
                let (dphidxsi, dphideta) = (self.space.dphi2dx)(xsi, eta);
                
                let mut dxdxsi = 0.0;
                let mut dydxsi = 0.0;
                let mut dxdeta = 0.0;
                let mut dydeta = 0.0;

                for i in 0..self.space.n as usize {
                    dxdxsi += dphidxsi[i] * x[i];
                    dydxsi += dphidxsi[i] * y[i];
                    dxdeta += dphideta[i] * x[i];
                    dydeta += dphideta[i] * y[i];
                }   

                let mut jac = dxdxsi * dydeta - dydxsi * dxdeta;
                if jac <= 0.0 {
                    Warning!("Negative Jacobian determinant");
                }
                jac = jac.abs();
                
                let mut dphidx = vec![0.0; self.space.n as usize];
                let mut dphidy = vec![0.0; self.space.n as usize];
                for i in 0..self.space.n as usize {
                    dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                    dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
                }

                if self.planar_strain_stress == FemElasticCase::Axisymmetric {
                    let mut x_interp = 0.0;
                    for i in 0..self.space.n as usize {
                        x_interp += phi[i] * x[i];
                    }
                    for i in 0..self.space.n as usize {
                        for j in 0..self.space.n as usize {
                            self.system.A[map_x[i]][map_x[j]] += (dphidx[i] * a * dphidx[j] * x_interp + dphidy[i] * c * dphidy[j] * x_interp + dphidx[i] * b * phi[j] + phi[i] * (b * dphidx[j] + a * phi[j] / x_interp)) * jac * w;
						    self.system.A[map_x[i]][map_y[j]] += (dphidx[i] * b * dphidy[j] * x_interp + dphidy[i] * c * dphidx[j] * x_interp + phi[i] * b * dphidy[j]) * jac * w;
						    self.system.A[map_y[i]][map_x[j]] += (dphidy[i] * b * dphidx[j] * x_interp + dphidx[i] * c * dphidy[j] * x_interp + dphidy[i] * b * phi[j]) * jac * w;
						    self.system.A[map_y[i]][map_y[j]] += (dphidy[i] * a * dphidy[j] * x_interp + dphidx[i] * c * dphidx[j] * x_interp) * jac * w;
                        }
                    }
                    for i in 0..self.space.n as usize {
                        self.system.B[map_x[i]] += phi[i] * gx * rho * x_interp * jac * w;
                        self.system.B[map_y[i]] += phi[i] * gy * rho * x_interp * jac * w;
                    }
                } else {
                    for i in 0..self.space.n as usize {
                        for j in 0..self.space.n as usize {
                            self.system.A[map_x[i]][map_x[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * w;
                            self.system.A[map_x[i]][map_y[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * w;
                            self.system.A[map_y[i]][map_x[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * w;
                            self.system.A[map_y[i]][map_y[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * w;
                        }
                    }
                    for i in 0..self.space.n as usize {
                        self.system.B[map_x[i]] += phi[i] * gx * rho * jac * w;
                        self.system.B[map_y[i]] += phi[i] * gy * rho * jac * w;
                    }
                }

            }
        }
    }

    pub fn apply_neumann(&mut self) {

        let n_local = 2;
        
        for i_bnd in 0..self.conditions.len() as usize {
            let condition = &self.conditions[i_bnd];
            let boundary_type = condition.boundary_type;
            if condition.is_neumann() {
                continue;
            }

            for i_elem in condition.domain.elements.iter() {
                let i_elem = *i_elem as usize;
                let mut map: [usize; 4] = [0; 4];
                let mut x = [0.0; 4];
                let mut y = [0.0; 4];

                for i in 0..n_local {
                    map[i] = self.geometry.edges.elements[i_elem * n_local + i] as usize;
                    x[i] = self.geometry.nodes.x[map[i] as usize];
                    y[i] = self.geometry.nodes.y[map[i] as usize];
                }

                let mut tx = x[1] - x[0];
                let mut ty = y[1] - y[0];
                let l = (tx * tx + ty * ty).sqrt();
                let jac = l / 2.0;
                let nx = ty/l;
                let ny = -tx/l;
                tx /= l;
                ty /= l;

                let mut f_x = 0.0;
                let mut f_y = 0.0;

                if boundary_type == FemBoundaryType::NeumannX {
                    f_x = condition.value1;
                } else if boundary_type == FemBoundaryType::NeumannY {
                    f_y = condition.value1;
                } else if boundary_type == FemBoundaryType::NeumannN {
                    f_x = condition.value1 * nx;
                    f_y = condition.value1 * ny;
                } else if boundary_type == FemBoundaryType::NeumannT {
                    f_x = condition.value1 * ty;
                    f_y = condition.value1 * tx;
                }

                for i_integ in 0..self.rule.n as usize {
                    let xsi = self.rule.xsi[i_integ];
                    let w = self.rule.w[i_integ];

                    let phi = (self.space_edge.phi)(xsi);

                    if self.planar_strain_stress == FemElasticCase::Axisymmetric {
                        let mut x_interp = 0.0;
                        for i in 0..self.space_edge.n as usize {
                            x_interp += phi[i] * x[i];
                        }
                        for i in 0..self.space_edge.n as usize {
                            self.system.B[2 * map[i] + 0] += phi[i] * jac * w * f_x * x_interp;
                            self.system.B[2 * map[i] + 1] += phi[i] * jac * w * f_y * x_interp;
                        }
                    } else {
                        for i in 0..self.space_edge.n as usize {
                            self.system.B[2 * map[i] + 0] += phi[i] * jac * w * f_x;
                            self.system.B[2 * map[i] + 1] += phi[i] * jac * w * f_y;
                        }
                    }
                    
                }
            }

        }

    }

    pub fn apply_dirichlet(&mut self) {

        let numbers = &self.geometry.nodes.numbers;
        for node in 0..self.geometry.nodes.x.len() as usize {
            let constrained_node = &self.constrained_nodes[node];
            if constrained_node.boundary_type == FemBoundaryType::Undefined {
                continue;
            }
            let boundary_type = constrained_node.boundary_type;

            if boundary_type == FemBoundaryType::DirichletX {
                self.system.constrain(numbers[2*node] as usize, constrained_node.value1);
            }
            if boundary_type == FemBoundaryType::DirichletY {
                self.system.constrain(numbers[2*node + 1] as usize, constrained_node.value1);
            }
            if boundary_type == FemBoundaryType::DirichletXY {
                self.system.constrain(numbers[2*node + 1] as usize, constrained_node.value1);
                self.system.constrain(numbers[2*node + 1] as usize, constrained_node.value2);
            }
            if boundary_type == FemBoundaryType::DirichletN {
                let value = constrained_node.value1;
                let nx = constrained_node.nx;
                let ny = constrained_node.ny;
                let tx = ny;
                let ty = -nx;

                let mut c_x = vec![0.0; self.system.size];
                let mut c_y = vec![0.0; self.system.size];
                let mut l_x = vec![0.0; self.system.size];
                let mut l_y = vec![0.0; self.system.size];

                for i in 0..self.system.size as usize {
                    c_x[i] = self.system.A[i][numbers[2*node] as usize];
                    c_y[i] = self.system.A[i][numbers[2*node + 1] as usize];
                    l_x[i] = self.system.A[numbers[2*node] as usize][i];
                    l_y[i] = self.system.A[numbers[2*node + 1] as usize][i];
                }

                let a_tt = tx * (tx * self.system.A[numbers[2*node] as usize][numbers[2*node] as usize] + ty * self.system.A[numbers[2*node + 1] as usize][numbers[2*node] as usize]) + ty * (tx * self.system.A[numbers[2*node] as usize][numbers[2*node + 1] as usize] + ty * self.system.A[numbers[2*node + 1] as usize][numbers[2*node + 1] as usize]);
                let a_tn = nx * (tx * self.system.A[numbers[2*node] as usize][numbers[2*node] as usize] + ty * self.system.A[numbers[2*node + 1] as usize][numbers[2*node] as usize]) + ny * (tx * self.system.A[numbers[2*node] as usize][numbers[2*node + 1] as usize] + ty * self.system.A[numbers[2*node + 1] as usize][numbers[2*node + 1] as usize]);
                let b_t = tx * self.system.B[numbers[2*node] as usize] + ty * self.system.B[numbers[2*node + 1] as usize];
                
                self.system.A[numbers[2*node] as usize][numbers[2*node] as usize] = nx * nx + a_tt * tx * tx;
                self.system.A[numbers[2*node] as usize][numbers[2*node + 1] as usize] = nx * ny + a_tt * tx * ty;
                self.system.A[numbers[2*node + 1] as usize][numbers[2*node] as usize] = ny * nx + a_tt * ty * tx;
                self.system.A[numbers[2*node + 1] as usize][numbers[2*node + 1] as usize] = ny * ny + a_tt * ty * ty;
                
                self.system.B[numbers[2*node] as usize] = nx * value + tx * (b_t - a_tn * value);
                self.system.B[numbers[2*node + 1] as usize] = ny * value + ty * (b_t - a_tn * value);
    
                for i in 0..self.system.size as usize {
                    if i != numbers[2*node] as usize && i != numbers[2*node + 1] as usize {
                        self.system.A[i][numbers[2*node] as usize] = tx * (tx * c_x[i] + ty * c_y[i]);
                        self.system.A[i][numbers[2*node + 1] as usize] = ty * (tx * c_x[i] + ty * c_y[i]);
                        self.system.A[numbers[2*node] as usize][i] = tx * (tx * l_x[i] + ty * l_y[i]);
                        self.system.A[numbers[2*node + 1] as usize][i] = ty * (tx * l_x[i] + ty * l_y[i]);
                        self.system.B[i] -= value * (nx * c_x[i] + ny * c_y[i]);
                    }
                }
            }

            if boundary_type == FemBoundaryType::DirichletT {

                let value = constrained_node.value1;
                let nx = constrained_node.nx;
                let ny = constrained_node.ny;
                let tx = ny;
                let ty = -nx;

                let mut c_x = vec![0.0; self.system.size];
                let mut c_y = vec![0.0; self.system.size];
                let mut l_x = vec![0.0; self.system.size];
                let mut l_y = vec![0.0; self.system.size];

                for i in 0..self.system.size as usize {
                    c_x[i] = self.system.A[i][numbers[2*node] as usize];
                    c_y[i] = self.system.A[i][numbers[2*node + 1] as usize];
                    l_x[i] = self.system.A[numbers[2*node] as usize][i];
                    l_y[i] = self.system.A[numbers[2*node + 1] as usize][i];
                }

                let a_nn = nx * (nx * self.system.A[numbers[2*node] as usize][numbers[2*node] as usize] + ny * self.system.A[numbers[2*node] as usize][numbers[2*node + 1] as usize]) + ny * (nx * self.system.A[numbers[2*node + 1] as usize][numbers[2*node] as usize] + ny * self.system.A[numbers[2*node + 1] as usize][numbers[2*node + 1] as usize]);
                let a_nt = nx * (tx * self.system.A[numbers[2*node] as usize][numbers[2*node] as usize] + ty * self.system.A[numbers[2*node + 1] as usize][numbers[2*node] as usize]) + ny * (tx * self.system.A[numbers[2*node] as usize][numbers[2*node + 1] as usize] + ty * self.system.A[numbers[2*node + 1] as usize][numbers[2*node + 1] as usize]);
                let b_n = nx * self.system.B[numbers[2*node] as usize] + ny * self.system.B[numbers[2*node + 1] as usize];
                
                self.system.A[numbers[2*node] as usize][numbers[2*node] as usize] = tx * tx + a_nn * nx * nx;
                self.system.A[numbers[2*node] as usize][numbers[2*node + 1] as usize] = tx * ty + a_nn * nx * ny;
                self.system.A[numbers[2*node + 1] as usize][numbers[2*node] as usize] = ty * tx + a_nn * ny * nx;
                self.system.A[numbers[2*node + 1] as usize][numbers[2*node + 1] as usize] = ty * ty + a_nn * ny * ny;

                self.system.B[numbers[2*node] as usize] = tx * value + nx * (b_n - a_nt * value);
                self.system.B[numbers[2*node + 1] as usize] = ty * value + ny * (b_n - a_nt * value);

                for i in 0..self.system.size as usize {
                    if i != numbers[2*node] as usize && i != numbers[2*node + 1] as usize {
                        self.system.A[i][numbers[2*node] as usize] = nx * (nx * c_x[i] + ny * c_y[i]);
                        self.system.A[i][numbers[2*node + 1] as usize] = ny * (nx * c_x[i] + ny * c_y[i]);
                        self.system.A[numbers[2*node] as usize][i] = nx * (nx * l_x[i] + ny * l_y[i]);
                        self.system.A[numbers[2*node + 1] as usize][i] = ny * (nx * l_x[i] + ny * l_y[i]);
                        self.system.B[i] -= value * (tx * c_x[i] + ty * c_y[i]);
                    }
                }
            }

            if boundary_type == FemBoundaryType::DirichletNT {
                let value_n = constrained_node.value1;
                let value_t = constrained_node.value2;
                let nx = constrained_node.nx;
                let ny = constrained_node.ny;
                let tx = ny;
                let ty = -nx;
                self.system.constrain(numbers[2*node + 0] as usize, value_n * nx + value_t * tx);
                self.system.constrain(numbers[2*node + 1] as usize, value_n * ny + value_t * ty);
            }
        }
    }

    pub fn assemble(&mut self) {
        self.system.clear();
        self.assemble_elements();
        self.apply_neumann();
        self.apply_dirichlet();
    }
}


impl FemElementType {
    pub fn mean(&self, x: Vec<f64>, y: Vec<f64>) -> Result<(f64, f64), Error> {
        return Ok(match self {
            case @ Self::FemTriangle => {
                let x_mean = (x[0] + x[1] + x[2]) / 3.0;
                let y_mean = (y[0] + y[1] + y[2]) / 3.0;
                (x_mean, y_mean)
            }
            case @ Self::FemQuad => {
                let x_mean = (x[0] + x[1] + x[2] + x[3]) / 4.0;
                let y_mean = (y[0] + y[1] + y[2] + y[3]) / 4.0;
                (x_mean, y_mean)
            }
            case @ _ => {
                return Err(anyhow::anyhow!("Unknown element type"));
            }
        });
       
    }
}