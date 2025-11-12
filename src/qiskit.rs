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


use qiskit_sys::qk_circuit_gate;
use std::ffi::{CStr, CString};

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
    None,
}

fn qk_to_qiskit_error(err: qiskit_sys::QkExitCode) -> QiskitError {
    match err {
        qiskit_sys::QkExitCode_QkExitCode_Success => QiskitError::Success,
        qiskit_sys::QkExitCode_QkExitCode_CInputError => QiskitError::CInputError,
        qiskit_sys::QkExitCode_QkExitCode_NullPointerError => QiskitError::NullPointerError,
        qiskit_sys::QkExitCode_QkExitCode_AlignmentError => QiskitError::AlignmentError,
        qiskit_sys::QkExitCode_QkExitCode_IndexError => QiskitError::IndexError,
        qiskit_sys::QkExitCode_QkExitCode_ArithmeticError => QiskitError::ArithmeticError,
        qiskit_sys::QkExitCode_QkExitCode_MismatchedQubits => QiskitError::MismatchedQubits,
        qiskit_sys::QkExitCode_QkExitCode_ExpectedUnitary => QiskitError::ExpectedUnitary,
        qiskit_sys::QkExitCode_QkExitCode_TargetError => QiskitError::TargetError,
        qiskit_sys::QkExitCode_QkExitCode_TargetInstAlreadyExists => {
            QiskitError::TargetInstAlreadyExists
        }
        qiskit_sys::QkExitCode_QkExitCode_TargetQargMismatch => QiskitError::TargetQargMismatch,
        qiskit_sys::QkExitCode_QkExitCode_TargetInvalidQargsKey => {
            QiskitError::TargetInvalidQargsKey
        }
        qiskit_sys::QkExitCode_QkExitCode_TargetInvalidInstKey => QiskitError::TargetInvalidInstKey,
        _ => QiskitError::None,
    }
}

/// The core representation of a quantum circuit.
pub struct QuantumCircuit {
    circuit: *mut qiskit_sys::QkCircuit,
}

impl QuantumCircuit {
    /// Create a new quantum circuit.
    ///
    /// Ex: Create a quantum circuit with 10 qubits and 10 classical bits:
    ///
    /// ```
    /// use qiskit_rs::QuantumCircuit;
    ///
    /// let qc = QuantumCircuit::new(10, 10);
    /// ```
    pub fn new(num_qubits: u32, num_clbits: u32) -> QuantumCircuit {
        let qc: *mut qiskit_sys::QkCircuit =
            unsafe { qiskit_sys::qk_circuit_new(num_qubits, num_clbits) };
        QuantumCircuit { circuit: qc }
    }
    /// Return the number of qubits in a QuantumCircuit.
    pub fn num_qubits(&mut self) -> u32 {
        unsafe { qiskit_sys::qk_circuit_num_qubits(self.circuit) }
    }
    /// Return the number of classical bits in a QuantumCircuit.
    pub fn num_clbits(&mut self) -> u32 {
        unsafe { qiskit_sys::qk_circuit_num_clbits(self.circuit) }
    }

    fn gate(&mut self, gate: qiskit_sys::QkGate, qubits: &[u32], params: &[f64]) -> QiskitError {
        let retval = if params.is_empty() {
            unsafe { qk_circuit_gate(self.circuit, gate, qubits.as_ptr(), std::ptr::null()) }
        } else {
            unsafe { qk_circuit_gate(self.circuit, gate, qubits.as_ptr(), params.as_ptr()) }
        };
        qk_to_qiskit_error(retval)
    }
    pub fn dcx(&mut self, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_DCX, &[qubit1, qubit2], &[])
    }
    pub fn ecr(&mut self, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_ECR, &[qubit1, qubit2], &[])
    }
    pub fn h(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_H, &[qubit], &[])
    }
    pub fn id(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_I, &[qubit], &[])
    }
    pub fn iswap(&mut self, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_ISwap, &[qubit1, qubit2], &[])
    }
    pub fn p(&mut self, theta: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_Phase, &[qubit], &[theta])
    }
    pub fn r(&mut self, theta: f64, phi: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_R, &[qubit], &[theta, phi])
    }
    pub fn rcccx(
        &mut self,
        control_qubit1: u32,
        control_qubit2: u32,
        control_qubit3: u32,
        target_qubit: u32,
    ) -> QiskitError {
        self.gate(
            qiskit_sys::QkGate_QkGate_RC3X,
            &[control_qubit1, control_qubit2, control_qubit3, target_qubit],
            &[],
        )
    }
    pub fn rccx(
        &mut self,
        control_qubit1: u32,
        control_qubit2: u32,
        target_qubit: u32,
    ) -> QiskitError {
        self.gate(
            qiskit_sys::QkGate_QkGate_RCCX,
            &[control_qubit1, control_qubit2, target_qubit],
            &[],
        )
    }
    pub fn rx(&mut self, theta: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_RX, &[qubit], &[theta])
    }
    pub fn rxx(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_RXX, &[qubit1, qubit2], &[theta])
    }
    pub fn ry(&mut self, theta: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_RY, &[qubit], &[theta])
    }
    pub fn ryy(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_RY, &[qubit1, qubit2], &[theta])
    }
    pub fn rz(&mut self, phi: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_RZ, &[qubit], &[phi])
    }
    pub fn rzx(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_RZX, &[qubit1, qubit2], &[theta])
    }
    pub fn rzz(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_RZZ, &[qubit1, qubit2], &[theta])
    }
    pub fn s(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_S, &[qubit], &[])
    }
    pub fn sdg(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_Sdg, &[qubit], &[])
    }
    pub fn sx(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_SX, &[qubit], &[])
    }
    pub fn sxdg(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_SXdg, &[qubit], &[])
    }
    pub fn t(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_T, &[qubit], &[])
    }
    pub fn tdg(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_Tdg, &[qubit], &[])
    }
    pub fn u(&mut self, theta: f64, phi: f64, lam: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_U, &[qubit], &[theta, phi, lam])
    }
    pub fn x(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_X, &[qubit], &[])
    }
    pub fn y(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_Y, &[qubit], &[])
    }
    pub fn z(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_Z, &[qubit], &[])
    }
    pub fn cx(&mut self, control_qubit: u32, target_qubit: u32) -> QiskitError {
        self.gate(
            qiskit_sys::QkGate_QkGate_CX,
            &[control_qubit, target_qubit],
            &[],
        )
    }
    pub fn measure(&mut self, qubit: u32, clbit: u32) -> QiskitError {
        let retval = unsafe { qiskit_sys::qk_circuit_measure(self.circuit, qubit, clbit) };
        qk_to_qiskit_error(retval)
    }
    pub fn add_quantum_register(&mut self, register: QuantumRegister) {
        unsafe { qiskit_sys::qk_circuit_add_quantum_register(self.circuit, register.register) };
    }
    pub fn add_classical_register(&mut self, register: ClassicalRegister) {
        unsafe { qiskit_sys::qk_circuit_add_classical_register(self.circuit, register.register) };
    }
    pub fn copy(&mut self) -> QuantumCircuit {
        QuantumCircuit {
            circuit: unsafe { qiskit_sys::qk_circuit_copy(self.circuit) },
        }
    }

    /// Return the number of instructions in the circuit
    pub fn num_instructions(&self) -> usize {
        unsafe { qiskit_sys::qk_circuit_num_instructions(self.circuit) }
    }

    /// Return an iterator of all the instructions in the circuit
    pub fn instructions(&self) -> impl ExactSizeIterator<Item = CircuitInstruction<'_>> + '_ {
        let num_inst = self.num_instructions();
        CircuitInstructions {
            len: num_inst,
            circuit: self,
            index: 0,
        }
    }
}

impl Drop for QuantumCircuit {
    fn drop(&mut self) {
        unsafe { qiskit_sys::qk_circuit_free(self.circuit) };
    }
}

pub struct QuantumRegister {
    register: *mut qiskit_sys::QkQuantumRegister,
}

impl QuantumRegister {
    pub fn new(num_qubits: u32, name: &str) -> QuantumRegister {
        let cname = CString::new(name).expect("String to CString conversion failed");
        let cname = cname.as_ptr();
        QuantumRegister {
            register: unsafe { qiskit_sys::qk_quantum_register_new(num_qubits, cname) },
        }
    }
}

impl Drop for QuantumRegister {
    fn drop(&mut self) {
        unsafe { qiskit_sys::qk_quantum_register_free(self.register) };
    }
}

pub struct ClassicalRegister {
    register: *mut qiskit_sys::QkClassicalRegister,
}

impl ClassicalRegister {
    pub fn new(num_clbits: u32, name: &str) -> ClassicalRegister {
        let cname = CString::new(name).expect("String to CString conversion failed");
        let cname = cname.as_ptr();
        ClassicalRegister {
            register: unsafe { qiskit_sys::qk_classical_register_new(num_clbits, cname) },
        }
    }
}

impl Drop for ClassicalRegister {
    fn drop(&mut self) {
        unsafe { qiskit_sys::qk_classical_register_free(self.register) };
    }
}

/// A view of an instruction in a [`QuantumCircuit`]
///
/// This struct contains references to all the standard data
/// about an instruction in the circuit.
#[derive(Debug)]
pub struct CircuitInstruction<'a> {
    /// The name of the operation for the instruction
    pub name: &'a str,
    /// The qubits the instruction acts upon
    pub qubits: &'a [u32],
    /// The clbits the instruction acts upon
    pub clbits: &'a [u32],
    /// The parameters for the instruction
    pub params: &'a [f64],
    inst: qiskit_sys::QkCircuitInstruction,
}

impl<'a> Drop for CircuitInstruction<'a> {
    fn drop(&mut self) {
        unsafe {
            qiskit_sys::qk_circuit_instruction_clear(&mut self.inst);
        }
    }
}

pub struct CircuitInstructions<'a> {
    len: usize,
    index: usize,
    circuit: &'a QuantumCircuit,
}

impl<'a> ExactSizeIterator for CircuitInstructions<'a> {
    fn len(&self) -> usize {
        self.len - self.index
    }
}

impl<'a> Iterator for CircuitInstructions<'a> {
    type Item = CircuitInstruction<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.len {
            return None;
        }
        let out = unsafe {
            let mut inst = qiskit_sys::QkCircuitInstruction {
                name: std::ptr::null_mut(),
                qubits: std::ptr::null_mut(),
                clbits: std::ptr::null_mut(),
                params: std::ptr::null_mut(),
                num_qubits: u32::MAX,
                num_clbits: u32::MAX,
                num_params: u32::MAX,
            };

            qiskit_sys::qk_circuit_get_instruction(self.circuit.circuit, self.index, &mut inst);
            let qubits = std::slice::from_raw_parts(inst.qubits, inst.num_qubits as usize);
            let clbits = std::slice::from_raw_parts(inst.clbits, inst.num_clbits as usize);
            let params = std::slice::from_raw_parts(inst.params, inst.num_params as usize);
            let name = CStr::from_ptr(inst.name).to_str().unwrap();
            Some(CircuitInstruction {
                name,
                qubits,
                clbits,
                params,
                inst,
            })
        };
        self.index += 1;
        out
    }
}

#[cfg(test)]
mod tests {
    use super::QuantumCircuit;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn test_circuit_instructions() {
        let mut qc = QuantumCircuit::new(100, 100);
        qc.rz(FRAC_PI_2, 0);
        qc.sx(0);
        qc.rz(FRAC_PI_2, 0);
        for target in 0..100u32 {
            qc.cx(0, target);
            qc.measure(target, target);
        }
        let res = qc.instructions();
        let mut target: u32 = 0;
        for (idx, inst) in res.enumerate() {
            if idx == 0 || idx == 2 {
                assert_eq!(inst.name, "rz");
                assert_eq!(&[0,], inst.qubits);
                assert_eq!(inst.clbits, &[]);
                assert_eq!(&[FRAC_PI_2,], inst.params);
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
}
