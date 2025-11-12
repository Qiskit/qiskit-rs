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

use qiskit_rs::QuantumCircuit;

#[test]
fn test_ghz() {
    let num_qubits = 10;
    let mut qc = QuantumCircuit::new(num_qubits, num_qubits);
    qc.h(0);
    for i in 0..(num_qubits - 1) {
        qc.cx(i, i + 1);
    }
    for i in 0..num_qubits {
        qc.measure(i, i);
    }
    assert_eq!(qc.num_qubits(), num_qubits);
    assert_eq!(qc.num_clbits(), num_qubits);
}
