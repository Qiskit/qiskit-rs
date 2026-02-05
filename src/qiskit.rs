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

    fn gate(&mut self, gate: qiskit_sys::QkGate, qubits: &[u32], params: &[f64]) -> QiskitError {
        let retval = if params.is_empty() {
            unsafe { qk_circuit_gate(self.circuit, gate, qubits.as_ptr(), std::ptr::null()) }
        } else {
            unsafe { qk_circuit_gate(self.circuit, gate, qubits.as_ptr(), params.as_ptr()) }
        };
        qk_to_qiskit_error(retval)
    }
    /// Apply a double-CNOT gate.
    pub fn dcx(&mut self, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_DCX, &[qubit1, qubit2], &[])
    }
    /// Apply an echoed cross-resonance gate.
    pub fn ecr(&mut self, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_ECR, &[qubit1, qubit2], &[])
    }
    /// Apply a Hadamard gate.
    pub fn h(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_H, &[qubit], &[])
    }
    /// Apply an Identity gate.
    pub fn id(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_I, &[qubit], &[])
    }
    /// Apply an iSWAP gate.
    pub fn iswap(&mut self, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_ISwap, &[qubit1, qubit2], &[])
    }
    /// Apply a Phase gate, a single-qubit rotation about the Z axis.
    pub fn p(&mut self, theta: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_Phase, &[qubit], &[theta])
    }
    /// Apply an RGate
    pub fn r(&mut self, theta: f64, phi: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_R, &[qubit], &[theta, phi])
    }
    /// Apply a simplified 3-controlled Toffoli gate.
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
        self.gate(
            qiskit_sys::QkGate_QkGate_RCCX,
            &[control_qubit1, control_qubit2, target_qubit],
            &[],
        )
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
        self.gate(qiskit_sys::QkGate_QkGate_RX, &[qubit], &[theta])
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
        self.gate(qiskit_sys::QkGate_QkGate_RXX, &[qubit1, qubit2], &[theta])
    }
    /// Apply a single-qubit rotation about the Y axis.
    pub fn ry(&mut self, theta: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_RY, &[qubit], &[theta])
    }
    /// Apply a 2-qubit rotation about YY.
    pub fn ryy(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_RY, &[qubit1, qubit2], &[theta])
    }
    /// Apply a single-qubit rotation about the Z axis.
    pub fn rz(&mut self, phi: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_RZ, &[qubit], &[phi])
    }
    /// Apply a 2-qubit rotation about ZX.
    pub fn rzx(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_RZX, &[qubit1, qubit2], &[theta])
    }
    /// Apply a 2-qubit rotation about ZX.
    pub fn rzz(&mut self, theta: f64, qubit1: u32, qubit2: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_RZZ, &[qubit1, qubit2], &[theta])
    }
    /// Apply a single qubit S gate.
    pub fn s(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_S, &[qubit], &[])
    }
    /// Apply a single qubit S-adjoint gate.
    pub fn sdg(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_Sdg, &[qubit], &[])
    }
    /// Apply a single-qubit Sqrt(X) gate.
    pub fn sx(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_SX, &[qubit], &[])
    }
    /// Apply an inverse single-qubit Sqrt(X) gate.
    pub fn sxdg(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_SXdg, &[qubit], &[])
    }
    /// Apply a single qubit T gate.
    pub fn t(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_T, &[qubit], &[])
    }
    /// Apply a single qubit T-adjoint gate.
    pub fn tdg(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_Tdg, &[qubit], &[])
    }
    /// Apply a generic single-qubit rotation.
    pub fn u(&mut self, theta: f64, phi: f64, lam: f64, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_U, &[qubit], &[theta, phi, lam])
    }
    /// Apply a single-qubit Pauli-X gate.
    pub fn x(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_X, &[qubit], &[])
    }
    /// Apply a single-qubit Pauli-Y gate.
    pub fn y(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_Y, &[qubit], &[])
    }
    /// Apply a single-qubit Pauli-Z gate.
    pub fn z(&mut self, qubit: u32) -> QiskitError {
        self.gate(qiskit_sys::QkGate_QkGate_Z, &[qubit], &[])
    }
    /// Apply a controlled-X gate.
    pub fn cx(&mut self, control_qubit: u32, target_qubit: u32) -> QiskitError {
        self.gate(
            qiskit_sys::QkGate_QkGate_CX,
            &[control_qubit, target_qubit],
            &[],
        )
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

fn bitterm_to_qkbitterm(bitterm: char) -> qiskit_sys::QkBitTerm {
    match bitterm {
        'X' => qiskit_sys::QkBitTerm_QkBitTerm_X,
        'Y' => qiskit_sys::QkBitTerm_QkBitTerm_Y,
        'Z' => qiskit_sys::QkBitTerm_QkBitTerm_Z,
        '+' => qiskit_sys::QkBitTerm_QkBitTerm_Plus,
        '-' => qiskit_sys::QkBitTerm_QkBitTerm_Minus,
        'r' => qiskit_sys::QkBitTerm_QkBitTerm_Right,
        'l' => qiskit_sys::QkBitTerm_QkBitTerm_Left,
        '0' => qiskit_sys::QkBitTerm_QkBitTerm_Zero,
        '1' => qiskit_sys::QkBitTerm_QkBitTerm_One,
        _ => panic!("Invalid bitterm: {}", bitterm.to_string()),
    }
}

/// An observable over Pauli bases that stores its data in a qubit-sparse format.
pub struct Observable {
    observable: *mut qiskit_sys::QkObs,
}

/// A complex, double-precision number representation.
pub type Complex64 = qiskit_sys::QkComplex64;

impl Observable {
    /// Create a new observable
    ///
    /// # Example
    ///
    /// ```
    /// use qiskit_rs::{Observable, Complex64};
    //
    /// let num_qubits = 100;
    /// let coeffs = [Complex64{re: 1.0, im: -1.0}, Complex64{re: 1.0, im: -1.0}];
    /// let bits = ['0', '1', '+', '-'];
    /// let indices = [0, 1, 98, 99];
    /// let boundaries = [0, 2, 4];
    ///
    /// let obs = Observable::new(num_qubits, &coeffs, &bits, &indices, &boundaries);
    /// ```
    pub fn new(
        num_qubits: u32,
        coeffs: &[Complex64],
        bit_terms: &[char],
        indices: &[u32],
        boundaries: &[usize],
    ) -> Observable {
        // Input shape checks (see https://quantum.cloud.ibm.com/docs/en/api/qiskit-c/qk-obs#representation)
        assert!(bit_terms.len() == indices.len());
        assert!(coeffs.len() + 1 == boundaries.len());

        let mut coeffs: Vec<qiskit_sys::QkComplex64> = Vec::from(coeffs);
        let mut bit_terms: Vec<u8> = bit_terms
            .into_iter()
            .map(|x| bitterm_to_qkbitterm(*x))
            .collect();
        let mut indices = Vec::from(indices);
        let mut boundaries = Vec::from(boundaries);
        Observable {
            observable: unsafe {
                qiskit_sys::qk_obs_new(
                    num_qubits,
                    coeffs.len().try_into().unwrap(),
                    bit_terms.len().try_into().unwrap(),
                    coeffs.as_mut_ptr(),
                    bit_terms.as_mut_ptr(),
                    indices.as_mut_ptr(),
                    boundaries.as_mut_ptr(),
                )
            },
        }
    }
    /// Construct the zero observable (without any terms).
    pub fn zero(num_qubits: u32) -> Observable {
        Observable {
            observable: unsafe { qiskit_sys::qk_obs_zero(num_qubits) },
        }
    }
    /// Construct the identity observable.
    pub fn identity(num_qubits: u32) -> Observable {
        Observable {
            observable: unsafe { qiskit_sys::qk_obs_identity(num_qubits) },
        }
    }
    /// Add two observables.
    pub fn add(&self, obs: &Observable) -> Observable {
        Observable {
            observable: unsafe { qiskit_sys::qk_obs_add(self.observable, obs.observable) },
        }
    }
    /// Multiply the observable by a complex coefficient.
    pub fn multiply(&self, coeff: &Complex64) -> Observable {
        Observable {
            observable: unsafe { qiskit_sys::qk_obs_multiply(self.observable, coeff) },
        }
    }
    /// Compose (multiply) two observables.
    pub fn compose(&self, obs: &Observable) -> Observable {
        Observable {
            observable: unsafe { qiskit_sys::qk_obs_compose(self.observable, obs.observable) },
        }
    }
    /// Compose (multiply) two observables according to a custom qubit order.
    pub fn compose_map(&self, obs: &Observable, qargs: &[u32]) -> Observable {
        Observable {
            observable: unsafe {
                qiskit_sys::qk_obs_compose_map(self.observable, obs.observable, qargs.as_ptr())
            },
        }
    }
    /// Get the number of terms in the observable.
    pub fn num_terms(&self) -> usize {
        unsafe { qiskit_sys::qk_obs_num_terms(self.observable) }
    }
    /// Get the number of qubits the observable is defined on.
    pub fn num_qubits(&self) -> u32 {
        unsafe { qiskit_sys::qk_obs_num_qubits(self.observable) }
    }
    /// Get the list of coefficients of an observable
    ///
    /// ```
    /// use qiskit_rs::{Observable, Complex64};
    ///
    /// let num_qubits = 100;
    /// let coeffs = [Complex64{re: 1.0, im: -1.0}, Complex64{re: 1.0, im: -1.0}];
    /// let bits = ['0', '1', '+', '-'];
    /// let indices = [0, 1, 98, 99];
    /// let boundaries = [0, 2, 4];
    ///
    /// let obs = Observable::new(num_qubits, &coeffs, &bits, &indices, &boundaries);
    ///
    /// for coef in obs.coeffs() {
    ///     print!("Coef: {{ re: {}, im: {} }}", coef.re, coef.im);
    /// }
    /// ```
    pub fn coeffs(&self) -> std::slice::Iter<'_, qiskit_sys::QkComplex64> {
        let num_coeffs: usize = unsafe { qiskit_sys::qk_obs_num_terms(self.observable) };
        let coeffs_ptr = unsafe { qiskit_sys::qk_obs_coeffs(self.observable) };
        let slice = unsafe { std::slice::from_raw_parts(coeffs_ptr, num_coeffs) };
        slice.iter()
    }
    /// Get the list of indices of an observable
    ///
    /// ```
    /// use qiskit_rs::{Observable, Complex64};
    /// let num_qubits = 100;
    /// let coeffs = [Complex64{re: 1.0, im: -1.0}, Complex64{re: 1.0, im: -1.0}];
    /// let bits = ['0', '1', '+', '-'];
    /// let indices = [0, 1, 98, 99];
    /// let boundaries = [0, 2, 4];
    /// let obs = Observable::new(num_qubits, &coeffs, &bits, &indices, &boundaries);
    /// for idx in obs.indices() {
    ///     print!("Index: {}", idx);
    /// }
    /// ```
    pub fn indices(&self) -> std::slice::Iter<'_, u32> {
        let num_indices: usize = unsafe { qiskit_sys::qk_obs_len(self.observable) };
        let indices_ptr = unsafe { qiskit_sys::qk_obs_indices(self.observable) };
        let slice = unsafe { std::slice::from_raw_parts(indices_ptr, num_indices) };
        slice.iter()
    }
    /// Get the number of bit terms/indices in the observable.
    pub fn len(&self) -> usize {
        unsafe { qiskit_sys::qk_obs_len(self.observable) }
    }
    /// Copy the observable
    pub fn copy(&self) -> Observable {
        Observable {
            observable: unsafe { qiskit_sys::qk_obs_copy(self.observable) },
        }
    }
    /// Compare two observables for equality.
    pub fn equal(&self, obs: &Observable) -> bool {
        unsafe { qiskit_sys::qk_obs_equal(self.observable, obs.observable) }
    }
    /// Return a string representation of the observable
    pub fn str(&self) -> String {
        let obs_str = unsafe { qiskit_sys::qk_obs_str(self.observable) };
        // Clone C string into String, which implements Drop
        let retval = String::from(unsafe { CStr::from_ptr(obs_str) }.to_str().unwrap());
        unsafe { qiskit_sys::qk_str_free(obs_str) };
        retval
    }
    /// Apply a new qubit layout to the observable.
    ///
    /// # Example
    ///
    /// ```
    /// use qiskit_rs::Observable;
    ///
    /// let identity = Observable::identity(2);
    ///
    /// // The number of qubits the observable acts on can be extended by
    /// // setting a larger num_qubits than the current observable has.
    /// let identity = identity.apply_layout(&[10, 9, 8, 7]);
    /// ```
    pub fn apply_layout(&self, layout: &[u32]) -> Observable {
        let new = self.copy();
        let _ = unsafe {
            qiskit_sys::qk_obs_apply_layout(new.observable, layout.as_ptr(), layout.len() as u32)
        };
        new
    }
    /// Calculate the canonical representation of the observable.
    ///
    /// # Example
    ///
    /// ```
    /// use qiskit_rs::Observable;
    ///
    /// let identity = Observable::identity(100);
    /// let two = identity.add(&identity);
    /// let canonical: Observable = identity.canonicalize(1e-6);
    /// ```
    pub fn canonicalize(&self, tolerance: f64) -> Observable {
        let new = self.copy();
        unsafe { qiskit_sys::qk_obs_canonicalize(new.observable, tolerance) };
        new
    }
}

impl Drop for Observable {
    fn drop(&mut self) {
        unsafe { qiskit_sys::qk_obs_free(self.observable) };
    }
}

#[cfg(test)]
mod tests {
    use super::{Complex64, Observable, QuantumCircuit};
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
    #[test]
    fn test_observables() {
        let num_qubits = 100;
        let coeffs = [
            Complex64 { re: 1.0, im: -1.0 },
            Complex64 { re: 1.0, im: -1.0 },
        ];
        let bits = ['0', '1', '+', '-'];
        let indices = [0, 1, 98, 99];
        let boundaries = [0, 2, 4];

        let obs = Observable::new(num_qubits, &coeffs, &bits, &indices, &boundaries);
        assert_eq!(obs.num_terms(), 2);

        let obs_b = Observable::new(num_qubits, &coeffs, &bits, &indices, &boundaries);
        assert!(obs.equal(&obs_b));

        assert_eq!(obs.str(), obs_b.str());
    }
    #[test]
    fn test_iterate_coeffs() {
        let num_qubits = 100;
        let coeffs = [
            Complex64 { re: 1.0, im: -1.0 },
            Complex64 { re: 1.0, im: -1.0 },
        ];
        let bits = ['0', '1', '+', '-'];
        let indices = [0, 1, 98, 99];
        let boundaries = [0, 2, 4];

        let obs = Observable::new(num_qubits, &coeffs, &bits, &indices, &boundaries);

        for (i, coef) in obs.coeffs().enumerate() {
            assert_eq!(coef.re, coeffs[i].re);
            assert_eq!(coef.im, coeffs[i].im);
        }
    }
    #[test]
    fn test_iterate_indices() {
        let num_qubits = 100;
        let coeffs = [
            Complex64 { re: 1.0, im: -1.0 },
            Complex64 { re: 1.0, im: -1.0 },
        ];
        let bits = ['0', '1', '+', '-'];
        let indices = [0, 1, 98, 99];
        let boundaries = [0, 2, 4];
        let obs = Observable::new(num_qubits, &coeffs, &bits, &indices, &boundaries);
        for (i, idx) in obs.indices().enumerate() {
            assert_eq!(*idx, indices[i]);
        }
    }
}
