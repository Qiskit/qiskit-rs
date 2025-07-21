use crate::qiskit_ffi;
use std::ffi::CString;

pub struct QuantumCircuit {
    circuit: *mut qiskit_ffi::QkCircuit,
}

impl QuantumCircuit {
    pub fn new(num_qubits: u32, num_clbits: u32) -> QuantumCircuit {
        let qc: *mut qiskit_ffi::QkCircuit = unsafe { qiskit_ffi::qk_circuit_new(num_qubits, num_clbits) };
        QuantumCircuit{
            circuit: qc,
        }
    }
    pub fn num_qubits(&mut self) -> u32 {
        unsafe { qiskit_ffi::qk_circuit_num_qubits(self.circuit) }
    }
    pub fn num_clbits(&mut self) -> u32 {
        unsafe { qiskit_ffi::qk_circuit_num_clbits(self.circuit) }
    }
    pub fn h(&mut self, qubit: u32) {
        unsafe { qiskit_ffi::qk_circuit_gate(self.circuit, qiskit_ffi::QkGate_QkGate_H, [qubit].as_ptr(), std::ptr::null()) };
    }
    pub fn cx(&mut self, control_qubit: u32, target_qubit: u32) -> u32 {
        unsafe { qiskit_ffi::qk_circuit_gate(self.circuit, qiskit_ffi::QkGate_QkGate_CX, [control_qubit, target_qubit].as_ptr(), std::ptr::null()) }
    }
    pub fn measure(&mut self, qubit: u32, clbit: u32) {
        unsafe { qiskit_ffi::qk_circuit_measure(self.circuit, qubit, clbit) };
    }
    pub fn add_quantum_register(&mut self, register: QuantumRegister) {
        unsafe { qiskit_ffi::qk_circuit_add_quantum_register(self.circuit, register.register) };
    }
    pub fn add_classical_register(&mut self, register: ClassicalRegister) {
        unsafe { qiskit_ffi::qk_circuit_add_classical_register(self.circuit, register.register) };
    }
    pub fn copy(&mut self) -> QuantumCircuit {
        QuantumCircuit{
            circuit: unsafe { qiskit_ffi::qk_circuit_copy(self.circuit) },
        }
    }
}

impl Drop for QuantumCircuit {
    fn drop(&mut self) {
        unsafe { qiskit_ffi::qk_circuit_free(self.circuit) };
    }
}

pub struct QuantumRegister {
    register: *mut qiskit_ffi::QkQuantumRegister,
}

impl QuantumRegister {
    pub fn new(num_qubits: u32, name: &str) -> QuantumRegister {
        let cname = CString::new(name)
            .expect("String to CString conversion failed");
        let cname = cname.as_ptr();
        QuantumRegister {
            register: unsafe { qiskit_ffi::qk_quantum_register_new(num_qubits, cname as *const i8) }
        }
    }
}

impl Drop for QuantumRegister {
    fn drop(&mut self) {
        unsafe { qiskit_ffi::qk_quantum_register_free(self.register) };
    }
}

pub struct ClassicalRegister {
    register: *mut qiskit_ffi::QkClassicalRegister,
}

impl ClassicalRegister {
    pub fn new(num_clbits: u32, name: &str) -> ClassicalRegister {
        let cname = CString::new(name)
            .expect("String to CString conversion failed");
        let cname = cname.as_ptr();
        ClassicalRegister {
            register: unsafe { qiskit_ffi::qk_classical_register_new(num_clbits, cname as *const i8) }
        }
    }
}

impl Drop for ClassicalRegister {
    fn drop(&mut self) {
        unsafe { qiskit_ffi::qk_classical_register_free(self.register) };
    }
}
