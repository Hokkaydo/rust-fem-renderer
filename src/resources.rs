use std::path::{Path, PathBuf};
use std::fs;
use std::io::{self, Read};
use std::ffi;

#[derive(Debug)]
pub enum Error {
    Io(io::Error),
    FileContainsNil,
    FailedToGetExePath
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Error {
        Error::Io(e)
    }
}

pub struct Resources {
    root_path: PathBuf
}

impl Resources {
    pub fn from_relative_exe_path(rel_path: &Path) -> Result<Resources, Error> {
        let exe_file_name = std::env::current_exe().map_err(|_| Error::FailedToGetExePath)?;
        let exe_path = exe_file_name.parent().ok_or(Error::FailedToGetExePath)?;
        
        Ok(Resources {
            root_path: exe_path.join(rel_path)
        })
    }

    pub fn resolve(&self, resource_name: &str) -> Result<ffi::CString, Error> {
        let mut file = fs::File::open(resource_name_to_path(&self.root_path, resource_name))?;
        let mut buffer: Vec<u8> = Vec::with_capacity(file.metadata()?.len() as usize + 1);
        file.read_to_end(&mut buffer)?;

        if buffer.iter().find(|x| **x == 0).is_some() {
            return Err(Error::FileContainsNil);
        }
        Ok(unsafe { ffi::CString::from_vec_unchecked(buffer) })
    }
}

fn resource_name_to_path(root_dir: &Path, location: &str) -> PathBuf {
    let mut path: PathBuf = root_dir.into();
    for part in location.split("/") {
        path = path.join(part);
    }
    path
}