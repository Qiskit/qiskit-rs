use git2::Repository;
use std::path::Path;
use std::process::Command;
use std::io::{self, Write};

fn build_qiskit() {
    
}

fn main() {
    let curr_dir_str = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let curr_dir: &Path = Path::new(&curr_dir_str);

    let qiskit_c_lib = curr_dir.join("qiskit_c_lib");
    let repo_dir: &Path = qiskit_c_lib.as_path();
    let repo_dir_str: &str = repo_dir.to_str().unwrap();

    println!("Cloning qiskit from source into {}", repo_dir_str);

    let url = "https://github.com/Qiskit/qiskit.git";
    let _ = Repository::clone(url, &repo_dir);

    println!("Generating dynamically linked qiskit libraries");

    let output = Command::new("make")
        .current_dir(repo_dir)
        .arg("c")
        .output()
        .unwrap();

    io::stdout().write_all(&output.stdout).unwrap();
    io::stderr().write_all(&output.stderr).unwrap();

    // let output = Command::new("cargo")
    //     .current_dir(repo_dir)
    //     .arg("rustc")
    //     .arg("--release")
    //     .arg("--crate-type")
    //     .arg("cdylib")
    //     .arg("-p")
    //     .arg("qiskit-cext")
    //     .args(["--target-dir", &format!("{}/different_target", env!("CARGO_MANIFEST_DIR"))])
    //     .output()
    //     .unwrap();
    // 
    // let _ = Command::new("mkdir")
    //     .current_dir(repo_dir)
    //     .arg("-p")
    //     .arg("dist/c/lib")
    //     .output()
    //     .unwrap();

    // let _ = Command::new("mkdir")
    //     .current_dir(repo_dir)
    //     .arg("-p")
    //     .arg("dist/c/include/qiskit")
    //     .output()
    //     .unwrap();

    // let _ = Command::new("cp")
    //     .current_dir(repo_dir)
    //     .arg("target/release/*")
    //     .arg("dist/c/lib")
    //     .output()
    //     .unwrap();

    // let _ = Command::new("cp")
    //     .current_dir(repo_dir)
    //     .arg("target/qiskit.h")
    //     .arg("dist/c/lib/qiskit.h")
    //     .output()
    //     .unwrap();

    // let _ = Command::new("cp")
    //     .current_dir(repo_dir)
    //     .arg("crates/cext/include/complex.h")
    //     .arg("dist/c/include/qiskit/complex.h")
    //     .output()
    //     .unwrap();

    println!("Dynamically linked libraries generated at {}", repo_dir_str); 

    println!("cargo:rustc-env=LD_LIBRARY_PATH={}/dist/c/lib", repo_dir_str);
    println!("cargo:rustc-link-search={}/dist/c/lib", repo_dir_str);
    println!("cargo:rustc-link-lib=qiskit");

    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("Unable to generate bindings");

    bindings
        .write_to_file("src/qiskit_ffi.rs")
        .expect("Couldn't write bindings!");
}
