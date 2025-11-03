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
