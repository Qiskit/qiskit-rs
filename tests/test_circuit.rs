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
