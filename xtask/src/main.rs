use std::{env, fs, path::Path, process::Command};

use anyhow::anyhow;
use toml_edit::DocumentMut;

static PYO3_FFI_VERSION: &str = "0.1.0";

/// Generate the `qiskit-pyo3-ffi` using the repo's Qiskit/qiskit submodule.
/// The generated crate's Cargo.toml is modified with the [PYO3_FFI_VERSION] version number.
fn generate_pyo3() -> anyhow::Result<()> {
    let manifest_dir = env::var("CARGO_MANIFEST_DIR")?;
    // where we find the Qiskit submodule
    let qiskit_path = Path::new(&manifest_dir).join("../qiskit");
    // where we write the qiskit-pyo3-ffi crate to
    let pyo3_path = Path::new(&manifest_dir).join("../qiskit-pyo3-ffi");

    // Clear the existing directory
    if pyo3_path.exists() {
        if pyo3_path.is_dir() {
            fs::remove_dir_all(&pyo3_path)?;
        } else {
            return Err(anyhow!(
                "Aborted due to unexpected file at: {:?}",
                pyo3_path
            ));
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
        .arg(&pyo3_path)
        .status()?;
    if !status.success() {
        return Err(anyhow!("generate-pyo3 failed"));
    }

    patch_manifest(&pyo3_path.join("Cargo.toml"), PYO3_FFI_VERSION)
}

/// Patch the Cargo.toml with the set version.
///
/// This also sets the repository to qiskit-rs.
///
/// # Args
///
/// * `cargo_path` - The full path to the `Cargo.toml` to patch, including the filename.
/// * `version` - The version string to use, e.g. "0.1.0". No verification is done on this.
fn patch_manifest(cargo_path: &Path, version: &str) -> anyhow::Result<()> {
    let mut doc: DocumentMut = fs::read_to_string(cargo_path)?.parse()?;
    doc["package"]["version"] = toml_edit::value(version);
    doc["package"]["repository"] = toml_edit::value("https://github.com/Qiskit/qiskit-rs");
    fs::write(cargo_path, doc.to_string())?;
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
