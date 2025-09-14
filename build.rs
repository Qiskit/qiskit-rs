use git2;
use std::path::Path;
use std::process::Command;
use std::env;

enum InstallMethod {
    CLONE,
    PATH(String)
}

// There are two installation methods:
// - Default (no path specified): Automatically clones and builds the qiskit c api from source
// - Path (Manually specified path): Uses qiskit c api binary or source from a path
//     Use envvar QISKIT_CEXT_PATH to specify the path to the qiskit binary
// By default, qiskit-rs uses Source (Automatic)
fn check_installation_method() -> InstallMethod {
    match env::var("QISKIT_CEXT_PATH") {
        Ok(val) => InstallMethod::PATH(val),
        Err(e) => match e {
            env::VarError::NotPresent => InstallMethod::CLONE,
            env::VarError::NotUnicode(env) => panic!("Envvar QISKIT_CEXT_PATH is not unicode: {env:?}")
        }
    }
}


fn clone_qiskit(source_path: &Path) {
    let url = "https://github.com/Qiskit/qiskit.git";
    let _ = match git2::Repository::clone(url, source_path) {
        Ok(_) => println!("Repository successfully cloned"),
        Err(e) => match e.code() {
            git2::ErrorCode::Exists => {println!("Repository already exists")}
            _ => panic!("Git clone failed: {e:?}")
        }
    };
}


fn build_qiskit(source_path: &Path) {  
    let _ = Command::new("make")
        .current_dir(source_path)
        .arg("c")
        .status()
        .expect("Dynamically linked library generation failed");
}


fn build_qiskit_from_source() {
    let out_dir = std::env::var("OUT_DIR").unwrap();
    let source_path = Path::new(&out_dir).join("qiskit_c_lib");
    let source_path = source_path.as_path();
   
    clone_qiskit(&source_path);

    match source_path.try_exists() {
        Ok(b) => {
            match b {
                true => {},
                false => panic!("Qiskit source path does not exist")
            }
        },
        Err(e) => panic!("{e:?}")
    }
   
    let repo_dir_str: &str = source_path
        .to_str()
        .unwrap();

    build_qiskit(&source_path);

    println!("cargo:rustc-env=LD_LIBRARY_PATH={}/dist/c/lib", repo_dir_str);
    println!("cargo:rustc-link-search={}/dist/c/lib", repo_dir_str);
    println!("cargo:rustc-link-lib=qiskit");

    let bindings = bindgen::Builder::default()
        .header(format!("{}/dist/c/include/qiskit.h", repo_dir_str))
        .header(format!("{}/dist/c/include/qiskit/complex.h", repo_dir_str))
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("Unable to generate bindings");

    bindings
        .write_to_file("src/qiskit_ffi.rs")
        .expect("Couldn't write bindings!");
}

fn build_qiskit_from_path(qiskit_path_str: String) {
    let qiskit_path = Path::new(&qiskit_path_str);

    match qiskit_path.try_exists() {
        Ok(b) => {
            match b {
                true => {},
                false => panic!("Qiskit path does not exist")
            }
        },
        Err(e) => panic!("{e:?}")
    }

    println!("cargo:rustc-env=LD_LIBRARY_PATH={}/dist/c/lib", qiskit_path_str);
    println!("cargo:rustc-link-search={}/dist/c/lib", qiskit_path_str);
    println!("cargo:rustc-link-lib=qiskit");

    let bindings = bindgen::Builder::default()
        .header(format!("{}/dist/c/include/qiskit.h", qiskit_path_str))
        .header(format!("{}/dist/c/include/qiskit/complex.h", qiskit_path_str))
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("Unable to generate bindings");

    bindings
        .write_to_file("src/qiskit_ffi.rs")
        .expect("Couldn't write bindings!");
}


fn main() { 
    println!("cargo:rerun-if-changed=build.rs");
    
    // println!("cargo::rerun-if-env-changed=QISKIT_CEXT_PATH");

    let install_method = check_installation_method();

    match install_method {
        InstallMethod::CLONE => {
            build_qiskit_from_source();
        },
        InstallMethod::PATH(path) => {
            build_qiskit_from_path(path);
        }
    };  
}
