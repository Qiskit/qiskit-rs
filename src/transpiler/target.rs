// This code is part of Qiskit Rust bindings.
//
// (C) Copyright IBM 2026
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

/// A mapping of instructions and properties representing the particular constraints of a backend.
pub struct Target {
    pub(crate) target: *mut qiskit_sys::QkTarget,
}

impl Target {
    pub fn new(num_qubits: u32) -> Target {
        let target = unsafe { qiskit_sys::qk_target_new(num_qubits) };
        Target { target: target }
    }
}

impl Drop for Target {
    fn drop(&mut self) {
        unsafe { qiskit_sys::qk_target_free(self.target) };
    }
}
