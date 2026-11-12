use std::{env, fs, path::Path, process::Command};

use anyhow::anyhow;

/// Generate the `qiskit-pyo3-ffi` using the repo's Qiskit/qiskit submodule.
fn generate_pyo3() -> anyhow::Result<()> {
    let manifest_dir = env::var("CARGO_MANIFEST_DIR")?;
    // where we find the Qiskit submodule
    let qiskit_path = Path::new(&manifest_dir).join("../qiskit");
    // where we write the final qiskit-pyo3-ffi crate to
    let pyo3_path = Path::new(&manifest_dir).join("../qiskit-pyo3-ffi");
    let pyo3_src_path = pyo3_path.join("src");
    // where we generate the crate into, and from where we'll copy it into the final location
    let build_path = pyo3_path.join("build");

    // Clear any existing `qiskit-pyo3-ffi/src` and `build` directories
    for to_remove in [&pyo3_src_path, &build_path] {
        if to_remove.exists() {
            if to_remove.is_dir() {
                fs::remove_dir_all(&to_remove)?;
            } else {
                return Err(anyhow!("Aborted due to unexpected file: {:?}", to_remove));
            }
        }
    }

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
        .arg(&build_path)
        .status()?;
    if !status.success() {
        return Err(anyhow!("generate-pyo3 failed"));
    }

    // Copy the src into the right location
    fs::rename(build_path.join("src"), pyo3_src_path)?;

    // Delete the build directory
    fs::remove_dir_all(&build_path)?;

    Ok(())
}

fn main() {
    let task = env::args().nth(1);
    let exit = match task.as_deref() {
        Some("generate-pyo3") => generate_pyo3(),
        Some(name) => Err(anyhow!("Invalid task: {name}")),
        None => Err(anyhow!("Usage: cargo xtask <task>")),
    };
    if let Err(e) = exit {
        eprintln!("{:?}", e);
        std::process::exit(1);
    }
}
