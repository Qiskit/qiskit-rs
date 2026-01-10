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

use qiskit_rs::circuit::{ClassicalRegister, QuantumCircuit, QuantumRegister};

#[test]
fn test_initialize_registers() {
    QuantumRegister::new(2, "qreg");
    ClassicalRegister::new(2, "creg");
}

#[test]
fn test_add_quantum_register() {
    let qreg = QuantumRegister::new(2, "qreg");
    let mut qc = QuantumCircuit::new(2, 0);
    qc.add_quantum_register(qreg);
}

#[test]
fn test_add_classical_register() {
    let creg = ClassicalRegister::new(2, "creg");
    let mut qc = QuantumCircuit::new(2, 0);
    qc.add_classical_register(creg);
}
