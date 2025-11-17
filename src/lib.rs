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

//! # qiskit-rs
//!
//! This library exposes the [C API for Qiskit](https://docs.quantum.ibm.com/api/qiskit-c) in Rust.
//!
//! ## Installation
//!
//! ```bash
//! cargo add --git https://github.com/Qiskit/qiskit-rs qiskit
//! ```
//!
//! ## Example
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
//!
//! ## Advanced Installation
//!
//! //! The Qiskit C API needs to be installed to use qiskit-rs. By default,
//! qiskit-rs uses the "clone" method.
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
//! #### Path
//!
//! Uses qiskit c api binary or pre-installed source directory from a path.
//!
//! ```bash
//! export QISKIT_CEXT_INSTALL_METHOD="path"
//! export QISKIT_CEXT_PATH="<path/to/qiskit-cext-dir>"
//! ```
//!
//! ## License
//!
//! [Apache License 2.0](https://github.com/Qiskit/qiskit/blob/main/LICENSE.txt)

#![warn(missing_docs)]
/// The main qiskit-rs module
pub mod qiskit;

pub use qiskit::{ClassicalRegister, QiskitError, QuantumCircuit, QuantumRegister, Observable, Complex64};
