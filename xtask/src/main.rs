use std::{env, path::Path, process::Command};

use anyhow::anyhow;

fn generate_pyo3() -> anyhow::Result<()> {
    let manifest_dir = env::var("CARGO_MANIFEST_DIR")?;
    // where we find the Qiskit submodule
    let qiskit_path = Path::new(&manifest_dir).join("../qiskit");
    // where we write the qiskit-pyo3-ffi crate to
    let pyo3_path = Path::new(&manifest_dir).join("../qiskit-pyo3-ffi");

    let cargo = env::var("CARGO").unwrap_or_else(|_| "cargo".to_string());
    let status = Command::new(cargo)
        .current_dir(&qiskit_path)
        .arg("run")
        .arg("--manifest-path")
        .arg(qiskit_path.join("Cargo.toml"))
        .args(["-p", "qiskit-bindgen-cli", "--", "generate-pyo3"])
        .arg("--cext-path")
        .arg(qiskit_path.join("crates/cext"))
        .arg("--output-path")
        .arg(&pyo3_path)
        .status()?;
    if !status.success() {
        return Err(anyhow!("generate-pyo3 failed"));
    }
    Ok(())
}

fn main() {
    let task = env::args().nth(1);
    let exit = match task.as_deref() {
        Some("generate-pyo3") => generate_pyo3(),
        Some(name) => Err(anyhow!("Invalid task: {name}")),
        None => Err(anyhow!("Usage: cargo xtask <task>")),
    };
    exit.unwrap_or_else(|e| eprintln!("{:?}", e));
}
