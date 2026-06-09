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

use qiskit_rs::{QiskitError, QuantumCircuit};

#[test]
fn test_single_qubit_gates() {
    let gate_funcs = [
        QuantumCircuit::h,
        QuantumCircuit::id,
        QuantumCircuit::s,
        QuantumCircuit::sdg,
        QuantumCircuit::sx,
        QuantumCircuit::sxdg,
        QuantumCircuit::t,
        QuantumCircuit::tdg,
        QuantumCircuit::x,
        QuantumCircuit::y,
        QuantumCircuit::z,
    ];

    for gate in gate_funcs {
        let mut qc = QuantumCircuit::new(1, 0);
        let ret = gate(&mut qc, 0);
        assert_eq!(ret, QiskitError::Success);
    }
}

#[test]
fn test_two_qubit_gates() {
    let gate_funcs = [
        QuantumCircuit::dcx,
        QuantumCircuit::ecr,
        QuantumCircuit::iswap,
        QuantumCircuit::cx,
    ];

    for gate in gate_funcs {
        let mut qc = QuantumCircuit::new(2, 0);
        let ret = gate(&mut qc, 0, 1);
        assert_eq!(ret, QiskitError::Success);
    }
}

#[test]
fn test_single_param_gates() {
    let gate_funcs = [
        QuantumCircuit::p,
        QuantumCircuit::rx,
        QuantumCircuit::ry,
        QuantumCircuit::rz,
    ];

    for gate in gate_funcs {
        let mut qc = QuantumCircuit::new(1, 0);
        let ret = gate(&mut qc, 0.0, 0);
        assert_eq!(ret, QiskitError::Success);
    }
}
