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

use std::ffi::CStr;

use qiskit_sys::QkTranspileResult;

use crate::QuantumCircuit;
use crate::Target;

pub struct Layout {
    pub layout: *mut qiskit_sys::QkTranspileLayout,
}

pub type TranspilerOptions = qiskit_sys::QkTranspileOptions;

pub fn transpile(
    qc: QuantumCircuit,
    target: Target,
    options: TranspilerOptions,
) -> (QuantumCircuit, Layout) {
    let mut qk_transpile_result: qiskit_sys::QkTranspileResult = QkTranspileResult {
        circuit: std::ptr::null_mut(),
        layout: std::ptr::null_mut(),
    };
    let mut error: *mut *mut i8 = std::ptr::null_mut();
    let error_code = unsafe {
        qiskit_sys::qk_transpile(
            qc.circuit,
            target.target,
            &options,
            &mut qk_transpile_result,
            error,
        )
    };
    if error_code != qiskit_sys::QkExitCode_QkExitCode_Success {
        let error_str = unsafe { *error };
        let error_msg = unsafe { CStr::from_ptr(error_str) };
        let error_msg: &str = error_msg.to_str().unwrap();
        panic!("ERROR: {}", error_msg);
    }
    (
        QuantumCircuit {
            circuit: qk_transpile_result.circuit,
        },
        Layout {
            layout: qk_transpile_result.layout,
        },
    )
}
