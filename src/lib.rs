//! Qiskit bindings for Rust
//!
//! This library exposes the C API for Qiskit in Rust.
//!
//!
//! ## Building
//!
//! The Qiskit C API needs to be installed to use qiskit-rs. There are two 
//! supported installation methods that can be specified:
//! 
//! ##### Path (Manually specified path): Uses qiskit c api binary or source from a path
//! ```bash
//! export QISKIT_CEXT_INSTALL_METHOD="path"
//! export QISKIT_CEXT_PATH="<path/to/qiskit-cext-dir>"
//! ```
//!
//! ##### Clone (no path specified): Automatically clones and builds the qiskit c api from source
//!       
//! WARNING, cloning and building from source is very slow.
//! ```bash
//! export QISKIT_CEXT_INSTALL_METHOD="clone"
//! ``` 
//!
//! 
//!
//! ## Creating a Circuit
//!
//! Create a simple bell state circuit in qiskit-rs
//!
//! ```
//! use qiskit_rs::QuantumCircuit;
//!
//! fn main() {
//!     // Initialize a circuit with 2 quantum registers and 2 classical registers
//!     let mut qc = QuantumCircuit::new(2, 2);
//!
//!     // Apply circuit operations
//!     qc.h(0);
//!     qc.cx(0, 1);
//!
//!     // Apply measurements
//!     qc.measure(0, 0);
//!     qc.measure(1, 1);
//! }
//! ```


#![allow(non_snake_case)]
#![allow(improper_ctypes)]
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]
// #[cfg(not(doctest))]
mod qiskit_ffi;
pub mod qiskit;

pub use qiskit::{
    QuantumCircuit,
    QuantumRegister,
    ClassicalRegister,
    QiskitError,
};
