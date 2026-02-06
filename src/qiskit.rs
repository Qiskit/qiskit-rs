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
use std::collections::BTreeSet;
use std::ffi::{CStr, CString};

#[derive(PartialEq, Eq, Debug)]
/// The error enum that enumerates the different error types possible from Qiskit.
pub enum QiskitError {
    /// Success.
    Success,

    /// Error related to data input.
    CInputError,

    /// Unexpected null pointer.
    NullPointerError,

    /// Pointer is not aligned to expected data.
    AlignmentError,

    /// Index out of bounds.
    IndexError,

    /// Error related to arithmetic operations or similar.
    ArithmeticError,

    /// Mismatching number of qubits.
    MismatchedQubits,

    /// Matrix is not unitary.
    ExpectedUnitary,

    /// Target related error.
    TargetError,

    /// Instruction already exists in the Target.
    TargetInstAlreadyExists,

    /// Properties with incorrect qargs was added.
    TargetQargMismatch,

    /// Trying to query into the target with non-existent qargs.
    TargetInvalidQargsKey,

    /// Querying an operation that doesn't exist in the Target.
    TargetInvalidInstKey,
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
        _ => panic!("Invalid option for QiskitError"),
    }
}

/// The core representation of a quantum circuit.
pub struct QuantumCircuit {
    circuit: *mut qiskit_sys::QkCircuit,
}

impl QuantumCircuit {
    /// Create a new quantum circuit.
    ///
    /// # Example
    ///
    /// Create a quantum circuit with 10 qubits and 10 classical bits:
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
    ///
    /// # Example
    ///
    /// ```
    /// use qiskit_rs::QuantumCircuit;
    ///
    /// let mut qc = QuantumCircuit::new(10, 10);
    /// let n = qc.num_qubits();
    /// ```
    pub fn num_qubits(&mut self) -> u32 {
        unsafe { qiskit_sys::qk_circuit_num_qubits(self.circuit) }
    }
    /// Return the number of classical bits in a QuantumCircuit.
    ///
    /// # Example
    ///
    /// ```
    /// use qiskit_rs::QuantumCircuit;
    ///
    /// let mut qc = QuantumCircuit::new(10, 10);
    /// let n = qc.num_clbits();
    /// ```
    pub fn num_clbits(&mut self) -> u32 {
        unsafe { qiskit_sys::qk_circuit_num_clbits(self.circuit) }
    }

    /// Apply a [QkGate] to the circuit.
    ///
    /// # Safety
    ///
    /// The caller must ensure that all qubit indices are less than num_qubits()
    /// and that no invalid combination of qubit indices is used.
    /// For example, [cx] can not use the same qubit index
    /// for both control and target qubit.
    ///
    /// The caller must also ensure that the correct number of parameters is passed
    /// for the [QkGate] used.
    unsafe fn gate(
        &mut self,
        gate: qiskit_sys::QkGate,
        qubits: &[u32],
        params: &[f64],
    ) -> QiskitError {
        let retval = if params.is_empty() {
            unsafe { qk_circuit_gate(self.circuit, gate, qubits.as_ptr(), std::ptr::null()) }
        } else {
            unsafe { qk_circuit_gate(self.circuit, gate, qubits.as_ptr(), params.as_ptr()) }
        };
        qk_to_qiskit_error(retval)
    }

    /// Check if all qubit indices are unique.
    fn unique_qubits_check(&self, qubits: &[u32]) -> bool {
        qubits.iter().collect::<BTreeSet<_>>().len() == qubits.len()
    }

    /// Check if all qubits are within the range of the circuit.
    ///
    /// Mutable as `num_qubits()` takes a mutable reference.
    fn qubits_in_range_check(&mut self, qubits: &[u32]) -> bool {
        qubits.iter().all(|q| q < &self.num_qubits())
    }

    /// Apply a double-CNOT gate.
    pub fn dcx(&mut self, qubit1: u32, qubit2: u32) -> QiskitError {
        if !self.unique_qubits_check(&[qubit1, qubit2])
            || !self.qubits_in_range_check(&[qubit1, qubit2])
        {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_DCX, &[qubit1, qubit2], &[]) }
    }

    /// Apply an echoed cross-resonance gate.
    pub fn ecr(&mut self, qubit1: u32, qubit2: u32) -> QiskitError {
        if !self.unique_qubits_check(&[qubit1, qubit2])
            || !self.qubits_in_range_check(&[qubit1, qubit2])
        {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_ECR, &[qubit1, qubit2], &[]) }
    }
    /// Apply a Hadamard gate.
    pub fn h(&mut self, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_H, &[qubit], &[]) }
    }
    /// Apply an Identity gate.
    pub fn id(&mut self, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_I, &[qubit], &[]) }
    }
    /// Apply an iSWAP gate.
    pub fn iswap(&mut self, qubit1: u32, qubit2: u32) -> QiskitError {
        if !self.unique_qubits_check(&[qubit1, qubit2])
            || !self.qubits_in_range_check(&[qubit1, qubit2])
        {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_ISwap, &[qubit1, qubit2], &[]) }
    }
    /// Apply a Phase gate, a single-qubit rotation about the Z axis.
    pub fn p(&mut self, theta: f64, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_Phase, &[qubit], &[theta]) }
    }
    /// Apply an RGate
    pub fn r(&mut self, theta: f64, phi: f64, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_R, &[qubit], &[theta, phi]) }
    }
    /// Apply a simplified 3-controlled Toffoli gate.
    pub fn rcccx(
        &mut self,
        control_qubit1: u32,
        control_qubit2: u32,
        control_qubit3: u32,
        target_qubit: u32,
    ) -> QiskitError {
        let qubits = &[control_qubit1, control_qubit2, control_qubit3, target_qubit];
        if !self.unique_qubits_check(qubits) || !self.qubits_in_range_check(qubits) {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_RC3X, qubits, &[]) }
    }
    /// Apply a simplified Toffoli gate.
    ///
    /// # Arguments
    ///
    /// * `control_qubit1`: First control qubit
    /// * `control_qubit2`: Second control qubit
    /// * `target_qubit`: Qubit to apply the gate to
    ///
    /// # Example
    ///
    /// ```
    /// use qiskit_rs::QuantumCircuit;
    /// use std::f64::consts::PI;
    ///
    /// let mut qc = QuantumCircuit::new(1, 1);
    /// qc.rx(PI / 2.0, 0);
    /// ```
    pub fn rccx(
        &mut self,
        control_qubit1: u32,
        control_qubit2: u32,
        target_qubit: u32,
    ) -> QiskitError {
        let qubits = &[control_qubit1, control_qubit2, target_qubit];
        if !self.unique_qubits_check(qubits) || !self.qubits_in_range_check(qubits) {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_RCCX, qubits, &[]) }
    }
    /// Apply a single-qubit rotation about the X axis.
    ///
    /// # Arguments
    ///
    /// * `theta`: Rotation angle
    /// * `qubit`: Qubit to apply the gate to
    ///
    /// # Example
    ///
    /// ```
    /// use qiskit_rs::QuantumCircuit;
    /// use std::f64::consts::PI;
    ///
    /// let mut qc = QuantumCircuit::new(1, 1);
    /// qc.rx(PI / 2.0, 0);
    /// ```
    pub fn rx(&mut self, theta: f64, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_RX, &[qubit], &[theta]) }
    }
    /// Apply a 2-qubit rotation about XX.
    ///
    /// # Arguments
    ///
    /// * `theta`: Rotation angle
    /// * `qubit1`: First qubit to apply the gate to
    /// * `qubit2`: Second qubit to apply the gate to
    ///
    /// # Example
    ///
    /// ```
    /// use qiskit_rs::QuantumCircuit;
    /// use std::f64::consts::PI;
    ///
    /// let mut qc = QuantumCircuit::new(2, 2);
    /// qc.rxx(PI / 2.0, 0, 1);
    /// ```
    pub fn rxx(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        if !self.unique_qubits_check(&[qubit1, qubit2])
            || !self.qubits_in_range_check(&[qubit1, qubit2])
        {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_RXX, &[qubit1, qubit2], &[theta]) }
    }

    /// Apply a single-qubit rotation about the Y axis.
    pub fn ry(&mut self, theta: f64, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_RY, &[qubit], &[theta]) }
    }
    /// Apply a 2-qubit rotation about YY.
    pub fn ryy(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        if !self.unique_qubits_check(&[qubit1, qubit2])
            || !self.qubits_in_range_check(&[qubit1, qubit2])
        {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_RYY, &[qubit1, qubit2], &[theta]) }
    }
    /// Apply a single-qubit rotation about the Z axis.
    pub fn rz(&mut self, phi: f64, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_RZ, &[qubit], &[phi]) }
    }
    /// Apply a 2-qubit rotation about ZX.
    pub fn rzx(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        if !self.unique_qubits_check(&[qubit1, qubit2])
            || !self.qubits_in_range_check(&[qubit1, qubit2])
        {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_RZX, &[qubit1, qubit2], &[theta]) }
    }
    /// Apply a 2-qubit rotation about ZX.
    pub fn rzz(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        if !self.unique_qubits_check(&[qubit1, qubit2])
            || !self.qubits_in_range_check(&[qubit1, qubit2])
        {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_RZZ, &[qubit1, qubit2], &[theta]) }
    }
    /// Apply a single qubit S gate.
    pub fn s(&mut self, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_S, &[qubit], &[]) }
    }
    /// Apply a single qubit S-adjoint gate.
    pub fn sdg(&mut self, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_Sdg, &[qubit], &[]) }
    }
    /// Apply a single-qubit Sqrt(X) gate.
    pub fn sx(&mut self, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_SX, &[qubit], &[]) }
    }
    /// Apply an inverse single-qubit Sqrt(X) gate.
    pub fn sxdg(&mut self, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_SXdg, &[qubit], &[]) }
    }
    /// Apply a single qubit T gate.
    pub fn t(&mut self, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_T, &[qubit], &[]) }
    }
    /// Apply a single qubit T-adjoint gate.
    pub fn tdg(&mut self, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_Tdg, &[qubit], &[]) }
    }
    /// Apply a generic single-qubit rotation.
    pub fn u(&mut self, theta: f64, phi: f64, lam: f64, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_U, &[qubit], &[theta, phi, lam]) }
    }
    /// Apply a single-qubit Pauli-X gate.
    pub fn x(&mut self, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_X, &[qubit], &[]) }
    }
    /// Apply a single-qubit Pauli-Y gate.
    pub fn y(&mut self, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_Y, &[qubit], &[]) }
    }
    /// Apply a single-qubit Pauli-Z gate.
    pub fn z(&mut self, qubit: u32) -> QiskitError {
        if qubit >= self.num_qubits() {
            return QiskitError::MismatchedQubits;
        }
        unsafe { self.gate(qiskit_sys::QkGate_QkGate_Z, &[qubit], &[]) }
    }
    /// Apply a controlled-X gate.
    pub fn cx(&mut self, control_qubit: u32, target_qubit: u32) -> QiskitError {
        if !self.unique_qubits_check(&[control_qubit, target_qubit])
            || !self.qubits_in_range_check(&[control_qubit, target_qubit])
        {
            return QiskitError::MismatchedQubits;
        }
        unsafe {
            self.gate(
                qiskit_sys::QkGate_QkGate_CX,
                &[control_qubit, target_qubit],
                &[],
            )
        }
    }
    /// Measure a qubit in the Z basis into a classical bit.
    pub fn measure(&mut self, qubit: u32, clbit: u32) -> QiskitError {
        let retval = unsafe { qiskit_sys::qk_circuit_measure(self.circuit, qubit, clbit) };
        qk_to_qiskit_error(retval)
    }
    /// Add a quantum register to the circuit.
    pub fn add_quantum_register(&mut self, register: QuantumRegister) {
        unsafe { qiskit_sys::qk_circuit_add_quantum_register(self.circuit, register.register) };
    }
    /// Add a classical register to the circuit.
    pub fn add_classical_register(&mut self, register: ClassicalRegister) {
        unsafe { qiskit_sys::qk_circuit_add_classical_register(self.circuit, register.register) };
    }
    /// Create a deepcopy of the circuit.
    pub fn copy(&mut self) -> QuantumCircuit {
        QuantumCircuit {
            circuit: unsafe { qiskit_sys::qk_circuit_copy(self.circuit) },
        }
    }

    /// Return the number of instructions in the circuit.
    pub fn num_instructions(&self) -> usize {
        unsafe { qiskit_sys::qk_circuit_num_instructions(self.circuit) }
    }

    /// Return an iterator of all the instructions in the circuit.
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

/// A quantum register.
pub struct QuantumRegister {
    register: *mut qiskit_sys::QkQuantumRegister,
}

impl QuantumRegister {
    /// Create a new quantum register.
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

/// A classical register.
pub struct ClassicalRegister {
    register: *mut qiskit_sys::QkClassicalRegister,
}

impl ClassicalRegister {
    /// Create a new classical register.
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

/// A list of circuit instructions for a QuantumCircuit
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
    use crate::QiskitError;
    use std::{f64::consts::FRAC_PI_2, u32};

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
}
