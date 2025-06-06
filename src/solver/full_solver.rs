use crate::{fem::FemProblem, Error};

impl FemProblem {
    pub fn full_solver(&mut self) {

        self.system.write_matrix("A_full.txt");
        self.system.write_vector("B_full.txt");
        for i in 0..self.system.size {
            if self.system.A[i][i].abs() < 1e-10 {
                Error!("Matrix is singular");
            }
            for j in i + 1..self.system.size {
                let ratio = self.system.A[j][i] / self.system.A[i][i];
                for k in i..self.system.size {
                    self.system.A[j][k] -= ratio * self.system.A[i][k];
                }
                self.system.B[j] -= ratio * self.system.B[i];
            }
        }   

        for i in (0..self.system.size).rev() {
            self.solution[i] = self.system.B[i];
            for j in i + 1..self.system.size {
                self.solution[i] -= self.system.A[i][j] * self.solution[j];
            }
            self.solution[i] /= self.system.A[i][i];
        }     
        // self.system.write_matrix("A_full.txt");
        // self.system.write_vector("B_full.txt");
    }
}