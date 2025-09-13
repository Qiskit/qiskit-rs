use git2;
use std::path::Path;
use std::process::Command;
use std::env;

enum InstallMethod {
    SOURCE(String),
    BIN(String)
}

fn check_installation_method() -> InstallMethod {
    let qiskit_rs_source_install: Option<String> = match env::var("QISKIT_RS_SOURCE_INSTALL") {
        Ok(val) => Some(val),
        Err(e) => match e {
            env::VarError::NotPresent => None,
            env::VarError::NotUnicode(env) => panic!("Envvar QISKIT_RS_SOURCE_INSTALL is not unicode: {env:?}")
        }
    };

    let qiskit_rs_binary_install: Option<String> = match env::var("QISKIT_CEXT_DIST") {
        Ok(val) => Some(val),
        Err(e) => match e {
            env::VarError::NotPresent => None,
            env::VarError::NotUnicode(env) => panic!("Envvar QISKIT_RS_BIN_INSTALL is not unicode: {env:?}")
        }
    };

    assert!(!(qiskit_rs_source_install.is_some() && qiskit_rs_binary_install.is_some()), "Only one installation method should be specified");

    if qiskit_rs_source_install.is_some() {
        return InstallMethod::SOURCE(qiskit_rs_source_install.unwrap());
    }
    InstallMethod::BIN(qiskit_rs_binary_install.unwrap())
}


fn build_qiskit_from_source() {
    let curr_dir_str = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let curr_dir: &Path = Path::new(&curr_dir_str);

    let qiskit_c_lib = curr_dir.join("qiskit_c_lib");
    let qiskit_source_dir: &Path = qiskit_c_lib.as_path();

    let repo_dir_str: &str = qiskit_source_dir
        .to_str()
        .expect("Qiskit source directory could not be found");

    println!("Cloning qiskit from source into {}", repo_dir_str);

    let url = "https://github.com/Qiskit/qiskit.git";
    let _ = match git2::Repository::clone(url, &qiskit_source_dir) {
        Ok(_) => println!("Repository successfully cloned"),
        Err(e) => match e.code() {
            git2::ErrorCode::Exists => {println!("Repository already exists")}
            _ => panic!("Git clone failed: {e:?}")
        }
    };

    println!("Generating dynamically linked qiskit libraries");

    let _ = Command::new("make")
        .current_dir(qiskit_source_dir)
        .arg("c")
        .status()
        .expect("Dynamically linked library generation failed");

    println!("Dynamically linked libraries generated at {}", repo_dir_str);

    println!("cargo:rustc-env=LD_LIBRARY_PATH={}/dist/c/lib", repo_dir_str);
    println!("cargo:rustc-link-search={}/dist/c/lib", repo_dir_str);
    println!("cargo:rustc-link-lib=qiskit");
}

fn build_qiskit_from_dist(dist_path_str: String) {
    let dist_path = Path::new(&dist_path_str);

    match dist_path.try_exists() {
        Ok(b) => {
            match b {
                true => {},
                false => panic!("Qiskit dist path does not exist")
            }
        },
        Err(e) => panic!("{e:?}")
    }

    println!("cargo:rustc-env=LD_LIBRARY_PATH={}/dist/c/lib", dist_path_str);
    println!("cargo:rustc-link-search={}/dist/c/lib", dist_path_str);
    println!("cargo:rustc-link-lib=qiskit");
}


fn main() { 
    println!("cargo:rerun-if-changed=build.rs");

    let install_method = check_installation_method();

    match install_method {
        InstallMethod::SOURCE(_) => {
            build_qiskit_from_source();
        },
        InstallMethod::BIN(dist_path) => {
            // panic!("Binary build not implemented")
            build_qiskit_from_dist(dist_path);
        }
    }; 

    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("Unable to generate bindings");

    bindings
        .write_to_file("src/qiskit_ffi.rs")
        .expect("Couldn't write bindings!");
}
