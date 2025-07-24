use qiskit_rs::{
    QuantumCircuit,
    QuantumRegister,
    ClassicalRegister,
};

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
