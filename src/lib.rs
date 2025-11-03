//! # Qiskit bindings for Rust
//!
//! This library exposes the C API for Qiskit in Rust.
//!
//!
//! ## Installation
//!
//! The Qiskit C API needs to be installed to use qiskit-rs. There are two
//! supported installation methods that can be specified:
//!
//! #### Path
//!
//! Uses qiskit c api binary or pre-installed source directory from a path.
//!
//! ```bash
//! export QISKIT_CEXT_INSTALL_METHOD="path"
//! export QISKIT_CEXT_PATH="<path/to/qiskit-cext-dir>"
//! ```
//!
//! #### Clone
//!
//! Automatically clones and builds the qiskit c api from source.
//!       
//! <div class="warning">Cloning and building from source is very slow!</div>
//!
//! ```bash
//! export QISKIT_CEXT_INSTALL_METHOD="clone"
//! ```
//!
//! ## Creating a Circuit
//!
//! Create a simple bell state circuit in qiskit-rs
//!
//! ```
//! use qiskit_rs::QuantumCircuit;
//!
//! // Initialize a circuit with 2 quantum registers and 2 classical registers
//! let mut qc = QuantumCircuit::new(2, 2);
//!
//! // Apply circuit operations
//! qc.h(0);
//! qc.cx(0, 1);
//!
//! // Apply measurements
//! qc.measure(0, 0);
//! qc.measure(1, 1);
//! ```

#![warn(missing_docs)]
pub mod qiskit;

pub use qiskit::{ClassicalRegister, QiskitError, QuantumCircuit, QuantumRegister};
