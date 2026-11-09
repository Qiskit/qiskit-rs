use std::path::Path;

use qiskit_sys_bindgen::*;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    let install_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    println!(
        "cargo:warning=Installation Path: {:?}",
        install_path.to_str()
    );
    generate_and_install_rust_pyo3_ffi(install_path).unwrap();
}
