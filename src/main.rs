#![allow(unused)]
// TODO: remove
mod fem;
mod glfem;
mod program;
mod solver;
use fem::{FemGeo, FemProblem, FemRenumType};
use std::rc::Rc;
use std::env;

fn job() -> FemProblem {
    let args: Vec<String> = env::args().collect();
    let suffix = if args.len() > 1 {
        args[1].as_str()
    } else {
        ""
    };
    let geo = FemGeo::new(&format!("./assets/mesh/mesh{suffix}.txt")).expect("exist");
    let mut problem = FemProblem::new(
        Rc::new(geo),
        &format!("./assets/problem/problem{suffix}.txt"),
        FemRenumType::No,
    )
    .expect("exist");

    // problem.solution = problem.system.gaussian_elimination().to_vec();
    // problem.write_solution(&format!("./assets/solution/solution{suffix}.txt"), Some(2)).expect("exist");
    problem.solution = problem.read_solution(&format!("./assets/solution/UV{suffix}.txt")).expect("exist");
    // problem.assemble();
    // problem.full_solver();
    // problem.write_solution(&format!("./assets/solution/UV_test.txt"), Some(2)).expect("exist");
    problem
}

fn main() {
    program::run(job);
}
