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

use std::f64::consts::FRAC_PI_2;

use qiskit_rs::{QiskitError, QuantumCircuit};
use qiskit_sys::qk_param_as_real;

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

#[test]
fn test_circuit_instructions() {
    let mut qc = QuantumCircuit::new(100, 100);
    qc.rz(FRAC_PI_2, 0);
    qc.sx(0);
    qc.rz(FRAC_PI_2, 0);
    for target in 1..100u32 {
        qc.cx(0, target);
        qc.measure(target, target);
    }
    let res = qc.instructions();
    let mut target: u32 = 1;
    for (idx, inst) in res.enumerate() {
        if idx == 0 || idx == 2 {
            assert_eq!(inst.name, "rz");
            assert_eq!(&[0,], inst.qubits);
            assert_eq!(inst.clbits, &[]);
            // TODO: Use a proper QkParam rust alternative.
            assert_eq!(
                &[FRAC_PI_2,],
                &[unsafe { qk_param_as_real(inst.params[0]) }]
            );
        } else if idx == 1 {
            assert_eq!(inst.name, "sx");
            assert_eq!(&[0,], inst.qubits);
            assert_eq!(inst.clbits, &[]);
            assert_eq!(inst.params, &[]);
        } else {
            let expected_name = if (idx - 3) % 2 == 0 { "cx" } else { "measure" };
            assert_eq!(expected_name, inst.name);
            assert_eq!(inst.params, &[]);
            if expected_name == "measure" {
                assert_eq!(inst.qubits, &[target]);
                assert_eq!(inst.clbits, &[target]);
                target += 1;
            } else {
                assert_eq!(inst.qubits, &[0, target]);
                assert_eq!(inst.clbits, &[]);
            }
        }
    }
}

#[test]
fn test_too_few_qubits_0() {
    let mut qc = QuantumCircuit::new(0, 0);
    assert_eq!(qc.id(0), QiskitError::MismatchedQubits);
    assert_eq!(qc.x(0), QiskitError::MismatchedQubits);
    assert_eq!(qc.y(0), QiskitError::MismatchedQubits);
    assert_eq!(qc.z(0), QiskitError::MismatchedQubits);
    assert_eq!(qc.h(0), QiskitError::MismatchedQubits);
    assert_eq!(qc.s(0), QiskitError::MismatchedQubits);
    assert_eq!(qc.sx(0), QiskitError::MismatchedQubits);
    assert_eq!(qc.sdg(0), QiskitError::MismatchedQubits);
    assert_eq!(qc.sxdg(0), QiskitError::MismatchedQubits);
    assert_eq!(qc.t(0), QiskitError::MismatchedQubits);
    assert_eq!(qc.tdg(0), QiskitError::MismatchedQubits);
    assert_eq!(qc.u(1., 1., 1., 0), QiskitError::MismatchedQubits);
    assert_eq!(qc.rx(1., 0), QiskitError::MismatchedQubits);
    assert_eq!(qc.ry(1., 0), QiskitError::MismatchedQubits);
    assert_eq!(qc.rz(1., 0), QiskitError::MismatchedQubits);
    assert_eq!(qc.p(1., 0), QiskitError::MismatchedQubits);
    assert_eq!(qc.r(1., 1., 0), QiskitError::MismatchedQubits);
    assert_eq!(qc.dcx(0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.ecr(0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.iswap(0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.rxx(1., 0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.ryy(1., 0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.rzz(1., 0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.rzx(1., 0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.cx(0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.rccx(0, 1, 2), QiskitError::MismatchedQubits);
    assert_eq!(qc.rcccx(0, 1, 2, 3), QiskitError::MismatchedQubits);
}

#[test]
fn test_too_few_qubits_1() {
    let mut qc = QuantumCircuit::new(1, 1);
    assert_eq!(qc.id(0), QiskitError::Success);
    assert_eq!(qc.x(0), QiskitError::Success);
    assert_eq!(qc.y(0), QiskitError::Success);
    assert_eq!(qc.z(0), QiskitError::Success);
    assert_eq!(qc.h(0), QiskitError::Success);
    assert_eq!(qc.s(0), QiskitError::Success);
    assert_eq!(qc.sx(0), QiskitError::Success);
    assert_eq!(qc.sdg(0), QiskitError::Success);
    assert_eq!(qc.sxdg(0), QiskitError::Success);
    assert_eq!(qc.t(0), QiskitError::Success);
    assert_eq!(qc.tdg(0), QiskitError::Success);
    assert_eq!(qc.u(1., 1., 1., 0), QiskitError::Success);
    assert_eq!(qc.rx(1., 0), QiskitError::Success);
    assert_eq!(qc.ry(1., 0), QiskitError::Success);
    assert_eq!(qc.rz(1., 0), QiskitError::Success);
    assert_eq!(qc.p(1., 0), QiskitError::Success);
    assert_eq!(qc.r(1., 1., 0), QiskitError::Success);
    assert_eq!(qc.dcx(0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.ecr(0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.iswap(0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.rxx(1., 0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.ryy(1., 0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.rzz(1., 0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.rzx(1., 0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.cx(0, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.rccx(0, 1, 2), QiskitError::MismatchedQubits);
    assert_eq!(qc.rcccx(0, 1, 2, 3), QiskitError::MismatchedQubits);
}

#[test]
fn test_too_few_qubits_2() {
    let mut qc = QuantumCircuit::new(2, 2);
    assert_eq!(qc.dcx(0, 1), QiskitError::Success);
    assert_eq!(qc.ecr(0, 1), QiskitError::Success);
    assert_eq!(qc.iswap(0, 1), QiskitError::Success);
    assert_eq!(qc.rxx(1., 0, 1), QiskitError::Success);
    assert_eq!(qc.ryy(1., 0, 1), QiskitError::Success);
    assert_eq!(qc.rzz(1., 0, 1), QiskitError::Success);
    assert_eq!(qc.rzx(1., 0, 1), QiskitError::Success);
    assert_eq!(qc.cx(0, 1), QiskitError::Success);
    assert_eq!(qc.rccx(0, 1, 2), QiskitError::MismatchedQubits);
    assert_eq!(qc.rcccx(0, 1, 2, 3), QiskitError::MismatchedQubits);
}

#[test]
fn test_too_few_qubits_3() {
    let mut qc = QuantumCircuit::new(3, 3);
    assert_eq!(qc.rccx(0, 1, 2), QiskitError::Success);
    assert_eq!(qc.rcccx(0, 1, 2, 3), QiskitError::MismatchedQubits);
}

#[test]
fn test_too_few_qubits_4() {
    let mut qc = QuantumCircuit::new(4, 4);
    assert_eq!(qc.rcccx(0, 1, 2, 3), QiskitError::Success);
}

#[test]
fn test_invalid_qubit_index() {
    let mut qc = QuantumCircuit::new(5, 0);
    assert_eq!(qc.id(u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.x(u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.y(u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.z(u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.h(u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.s(u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.sx(u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.sdg(u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.sxdg(u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.t(u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.tdg(u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.u(1., 1., 1., u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.rx(1., u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.ry(1., u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.rz(1., u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.p(1., u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.r(1., 1., u32::MAX), QiskitError::MismatchedQubits);
    assert_eq!(qc.dcx(u32::MAX, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.ecr(u32::MAX, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.iswap(u32::MAX, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.rxx(1., u32::MAX, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.ryy(1., u32::MAX, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.rzz(1., u32::MAX, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.rzx(1., u32::MAX, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.cx(u32::MAX, 1), QiskitError::MismatchedQubits);
    assert_eq!(qc.rccx(u32::MAX, 1, 2), QiskitError::MismatchedQubits);
    assert_eq!(qc.rcccx(u32::MAX, 1, 2, 3), QiskitError::MismatchedQubits);
}
