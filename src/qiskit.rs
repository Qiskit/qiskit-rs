use crate::qiskit_ffi;
use crate::qiskit_ffi::{qk_circuit_gate,QkExitCode};
use std::ffi::CString;

#[derive(PartialEq, Eq, Debug)]
pub enum QiskitError {
    Success,
    CInputError,
    NullPointerError,
    AlignmentError,
    IndexError,
    ArithmeticError,
    MismatchedQubits,
    ExpectedUnitary,
    TargetError,
    TargetInstAlreadyExists,
    TargetQargMismatch,
    TargetInvalidQargsKey,
    TargetInvalidInstKey,
    None
}


fn qk_to_qiskit_error(err: qiskit_ffi::QkExitCode) -> QiskitError {
    match err {
        qiskit_ffi::QkExitCode_QkExitCode_Success                   => QiskitError::Success,
        qiskit_ffi::QkExitCode_QkExitCode_CInputError               => QiskitError::CInputError,
        qiskit_ffi::QkExitCode_QkExitCode_NullPointerError          => QiskitError::NullPointerError,
        qiskit_ffi::QkExitCode_QkExitCode_AlignmentError            => QiskitError::AlignmentError,
        qiskit_ffi::QkExitCode_QkExitCode_IndexError                => QiskitError::IndexError,
        qiskit_ffi::QkExitCode_QkExitCode_ArithmeticError           => QiskitError::ArithmeticError,
        qiskit_ffi::QkExitCode_QkExitCode_MismatchedQubits          => QiskitError::MismatchedQubits,
        qiskit_ffi::QkExitCode_QkExitCode_ExpectedUnitary           => QiskitError::ExpectedUnitary,
        qiskit_ffi::QkExitCode_QkExitCode_TargetError               => QiskitError::TargetError,
        qiskit_ffi::QkExitCode_QkExitCode_TargetInstAlreadyExists   => QiskitError::TargetInstAlreadyExists,
        qiskit_ffi::QkExitCode_QkExitCode_TargetQargMismatch        => QiskitError::TargetQargMismatch,
        qiskit_ffi::QkExitCode_QkExitCode_TargetInvalidQargsKey     => QiskitError::TargetInvalidQargsKey,
        qiskit_ffi::QkExitCode_QkExitCode_TargetInvalidInstKey      => QiskitError::TargetInvalidInstKey,
        _ => QiskitError::None
    }
}


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
    fn gate(&mut self, gate: qiskit_ffi::QkGate, qubits: &[u32], params: &[f64]) -> QiskitError {
        let retval: QkExitCode;
        if params.len() == 0 {
            retval = unsafe { qk_circuit_gate(self.circuit, gate, qubits.as_ptr(), std::ptr::null()) };
        } else {
            retval = unsafe { qk_circuit_gate(self.circuit, gate, qubits.as_ptr(), params.as_ptr()) };
        }
        qk_to_qiskit_error(retval)
    }
    pub fn dcx(&mut self, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_DCX, &[qubit1, qubit2], &[])
    }
    pub fn ecr(&mut self, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_ECR, &[qubit1, qubit2], &[])
    }
    pub fn h(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_H, &[qubit], &[])
    }
    pub fn id(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_I, &[qubit], &[])
    }
    pub fn iswap(&mut self, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_ISwap, &[qubit1, qubit2], &[])
    }
    pub fn p(&mut self, theta: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_Phase, &[qubit], &[theta])
    }
    pub fn r(&mut self, theta: f64, phi: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_R, &[qubit], &[theta, phi])
    }
    pub fn rcccx(&mut self, control_qubit1: u32, control_qubit2: u32, control_qubit3: u32, target_qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_RC3X, &[control_qubit1, control_qubit2, control_qubit3, target_qubit], &[])
    }
    pub fn rccx(&mut self, control_qubit1: u32, control_qubit2: u32, target_qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_RCCX, &[control_qubit1, control_qubit2, target_qubit], &[])
    }
    pub fn rx(&mut self, theta: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_RX, &[qubit], &[theta])
    } 
    pub fn rxx(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_RXX, &[qubit1, qubit2], &[theta])
    }
    pub fn ry(&mut self, theta: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_RY, &[qubit], &[theta])
    }
    pub fn ryy(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_RY, &[qubit1, qubit2], &[theta])
    }
    pub fn rz(&mut self, phi: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_RZ, &[qubit], &[phi])
    }
    pub fn rzx(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_RZX, &[qubit1, qubit2], &[theta])
    }
    pub fn rzz(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_RZZ, &[qubit1, qubit2], &[theta])
    }
    pub fn s(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_S, &[qubit], &[])
    }
    pub fn sdg(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_Sdg, &[qubit], &[])
    }
    pub fn sx(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_SX, &[qubit], &[])
    }
    pub fn sxdg(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_SXdg, &[qubit], &[])
    }
    pub fn t(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_T, &[qubit], &[])
    }
    pub fn tdg(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_Tdg, &[qubit], &[])
    }
    pub fn u(&mut self, theta: f64, phi: f64, lam: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_U, &[qubit], &[theta, phi, lam])
    }
    pub fn x(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_X, &[qubit], &[])
    }
    pub fn y(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_Y, &[qubit], &[])
    }
    pub fn z(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_Z, &[qubit], &[])
    }
    pub fn cx(&mut self, control_qubit: u32, target_qubit: u32) -> QiskitError {
        self.gate(qiskit_ffi::QkGate_QkGate_CX, &[control_qubit, target_qubit], &[])
    }
    pub fn measure(&mut self, qubit: u32, clbit: u32) -> QiskitError {
        let retval = unsafe { qiskit_ffi::qk_circuit_measure(self.circuit, qubit, clbit) };
        qk_to_qiskit_error(retval)
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
