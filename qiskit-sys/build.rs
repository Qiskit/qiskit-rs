// This code is part of Qiskit Rust bindings.
//
// (C) Copyright IBM 2025
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

enum InstallMethod {
    Clone,
    Path(String),
}

#[derive(Debug)]
struct CargoCallbacks;

impl bindgen::callbacks::ParseCallbacks for CargoCallbacks {
    fn process_comment(&self, comment: &str) -> Option<String> {
        Some(format!("````ignore\n{}\n````", comment))
    }
}

// There are two installation methods:
// - Clone (no path specified): Automatically clones and builds the qiskit c api from source
//     Set envvar export QISKIT_CEXT_INSTALL_METHOD="clone" to use the clone method. WARNING, cloning and building from
//     source is very slow.
// - Path (Manually specified path): Uses qiskit c api binary or source from a path
//     export QISKIT_CEXT_INSTALL_METHOD="path"
//     export QISKIT_CEXT_PATH="path/to/qiskit-cext-dir"
fn check_installation_method() -> InstallMethod {
    let qiskit_cext_path = env::var("QISKIT_CEXT_PATH");
    match env::var("QISKIT_CEXT_INSTALL_METHOD") {
        Ok(val) => match val.as_str() {
            "path" => InstallMethod::Path(qiskit_cext_path.expect("QISKIT_CEXT_PATH is unset")),
            "clone" => InstallMethod::Clone,
            _ => panic!(
                "\"{}\" is not a valid input to QISKIT_CEXT_INSTALL_METHOD, please specify one of the following options: (\"path\", \"clone\")",
                val
            ),
        },
        Err(e) => match e {
            env::VarError::NotPresent => InstallMethod::Clone,
            env::VarError::NotUnicode(_) => {
                panic!("Envvar QISKIT_CEXT_INSTALL_METHOD is not unicode")
            }
        },
    }
}

fn build_qiskit(source_path: &Path, new_path: &Path) {
    Command::new("make")
        .current_dir(source_path)
        .arg("c")
        .status()
        .expect("Dynamically linked library generation failed");

    let new_to_dist = new_path.join("dist");

    // We must copy the generated dynamic library and header files
    // to our output folder so it is automatically included in the
    // rpath and cargo is able to correctly link to the dynamic
    // library and the header files.
    Command::new("cp")
        .current_dir(source_path.join("dist"))
        .args([
            "-r",
            ".",
            new_to_dist
                .to_str()
                .expect("Path should be convertible to string"),
        ])
        .status()
        .expect("Source path for c files should exist");
}

fn build_qiskit_from_source() {
    let out_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let source_path = Path::new(&out_dir).join("qiskit");
    let source_path = source_path.as_path();

    match source_path.try_exists() {
        Ok(b) => match b {
            true => {}
            false => panic!("Qiskit source path does not exist"),
        },
        Err(e) => panic!("{e:?}"),
    }

    // Path to which we will copy the generated c files.
    let new_path = std::env::var("OUT_DIR").expect("OUT_DIR env variable should have been set.");
    let new_path = Path::new(&new_path);

    build_qiskit(source_path, new_path);
    let repo_dir_str = new_path.to_str().unwrap();
    generate_bindings(repo_dir_str);
}

fn generate_bindings(qiskit_path_str: &str) {
    let qiskit_path = Path::new(&qiskit_path_str);

    match qiskit_path.try_exists() {
        Ok(b) => match b {
            true => {}
            false => panic!("Qiskit path does not exist"),
        },
        Err(e) => panic!("{e:?}"),
    }

    println!("cargo:rustc-link-search={}/dist/c/lib", qiskit_path_str);
    println!("cargo:rustc-link-lib=qiskit");

    let bindings = bindgen::Builder::default()
        .header(format!("{}/dist/c/include/qiskit.h", qiskit_path_str))
        // Link the include folder directly to clang (this includes complex.h and others).
        .clang_arg(format!("-I{}/dist/c/include/", qiskit_path_str))
        .parse_callbacks(Box::new(CargoCallbacks))
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}

/// Ensures that the version of Qiskit being checked out matches the
/// package version of ``qiskit-sys``.
fn check_or_update_submodule() -> Result<(), git2::Error> {
    let main_repo = git2::Repository::open_from_env()?;

    println!("{:?}", main_repo.commondir());

    let mut submodule = main_repo.find_submodule("qiskit")?;

    // Check if the repository is initialized and update it if not before checking out
    let repo = if let Ok(repo) = submodule.open() {
        repo
    } else {
        submodule.update(true, None).expect("WHAT?");
        submodule.open()?
    };

    let current_commit = repo.head()?;

    let tag_string: String = env!("CARGO_PKG_VERSION").into();

    // If current version is called "dev" do not update version in use.
    if !tag_string.contains("dev") {
        let (tag, Some(reference)) = repo.revparse_ext(&tag_string)? else {
            return Err(git2::Error::from_str(&format!(
                "Could not find tag referring to {}",
                tag_string
            )));
        };

        let ref_as_commit = reference.peel_to_commit()?;
        let curr_as_commit = current_commit.peel_to_commit()?;
        if ref_as_commit.id() != curr_as_commit.id() {
            println!(
                "cargo::warning=Current commit '{:?}' does not reference the package associated tag '{}'. The submodule will now checkout the correct tag.",
                curr_as_commit.id(),
                tag.as_tag().unwrap().name().unwrap()
            );
            repo.reset(ref_as_commit.as_object(), git2::ResetType::Hard, None)?;
        }
    } else {
        println!(
            "cargo::warning= The current version of qiskit-sys {:?} is in developer mode and may be prone to fail if the incorrect version of Qiskit is linked.",
            tag_string
        );
    }

    Ok(())
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo::rerun-if-env-changed=QISKIT_CEXT_INSTALL_METHOD");
    println!("cargo::rerun-if-env-changed=QISKIT_CEXT_PATH");
    let install_method = check_installation_method();

    match install_method {
        InstallMethod::Clone => {
            println!("cargo::warning=Cloning and building from source is very slow");
            check_or_update_submodule().expect("Submodule doesn't exist");
            build_qiskit_from_source();
        }
        InstallMethod::Path(path) => {
            generate_bindings(&path);
        }
    };
}
