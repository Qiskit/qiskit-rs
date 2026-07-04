// This code is part of Qiskit Rust bindings.
//
// (C) Copyright IBM 2026
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

use qiskit_rs::{QuantumCircuit, Target, transpile, transpiler::TranspilerOptions};

#[test]
fn test_transpile() {
    let qc = QuantumCircuit::new(2, 2);
    let options = TranspilerOptions {
        optimization_level: 0,
        seed: 0,
        approximation_degree: 0.0,
    };
    let target = Target::new(2);
    let (isa_circuit, layout) = transpile(qc, target, options);
    assert_eq!(isa_circuit.num_qubits(), 2);
}
