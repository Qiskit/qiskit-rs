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

use std::ffi::CString;

/// A quantum register.
pub struct QuantumRegister {
    pub(super) register: *mut qiskit_sys::QkQuantumRegister,
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
    pub(super) register: *mut qiskit_sys::QkClassicalRegister,
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
