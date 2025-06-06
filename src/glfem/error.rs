#[macro_export]
macro_rules! Error {
    ($a:expr) => {
        crate::glfem::error::fem_error($a, line!(), file!())
    };
}

#[macro_export]
macro_rules! ErrorScan {
    ($a:expr) => {
        crate::glfem::error::fem_error_scan($a, line!(), file!())
    };
}

#[macro_export]
macro_rules! Warning {
    ($a:expr) => {
        crate::glfem::error::fem_warning($a, line!(), file!())
    };
}

pub fn fem_error(message: &str, line: u32, file: &str) {
    println!("\n-------------------------------------------------------------------------------- ");
    println!(
        "  Error in {}:{} at line {} : \n {}",
        file, line, line, message
    );
    println!("--------------------------------------------------------------------- Yek Yek !!\n");
}

pub fn fem_error_scan(test: i32, line: u32, file: &str) {
    if test >= 0 {
        return;
    }
    println!("\n-------------------------------------------------------------------------------- ");
    println!(
        "  Error in scan or read in {}:{} at line {} : ",
        file, line, line
    );
    println!("--------------------------------------------------------------------- Yek Yek !!\n");
}

pub fn fem_warning(message: &str, line: u32, file: &str) {
    println!("\n-------------------------------------------------------------------------------- ");
    println!(
        "  Warning in {}:{} at line {} : \n {}",
        file, line, line, message
    );
    println!("--------------------------------------------------------------------- Yek Yek !!\n");
}
