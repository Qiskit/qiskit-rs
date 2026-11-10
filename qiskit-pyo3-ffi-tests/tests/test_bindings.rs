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

use pyo3::prelude::*;
use qiskit_pyo3_ffi as ffi;

fn import(py: Python<'_>) {
    ffi::qk_import(py).expect("`import qiskit` failed; is Qiskit installed in this Python env?");
}

/// Convert the version string (e.g. "2.5.1") to the HEX form, 0x020501F0.
/// Supports final releases and single-digit RC releases (i.e. does not support rc10 and higher)
fn version_string_to_hex(version: &str) -> u32 {
    let mut iter = version.chars();
    // Get the <major>.<minor>.<bugfix> parts
    let major = iter.next().unwrap().to_digit(10).unwrap();
    iter.next();
    let minor = iter.next().unwrap().to_digit(10).unwrap();
    iter.next();
    let bugfix = iter.next().unwrap().to_digit(10).unwrap();

    // Get the suffix
    let suffix_str = iter.collect::<String>();
    let suffix = if suffix_str.is_empty() {
        // final version
        0xF0
    } else if let Some(rc) = suffix_str.strip_prefix("rc")
        && rc.len() == 1
    {
        0xC0 | rc.chars().next().unwrap().to_digit(10).unwrap()
    } else {
        todo!("Unsupported suffix {suffix_str}");
    };

    (major << 24) | (minor << 16) | (bugfix << 8) | suffix
}

#[test]
fn test_circuit_capsule() {
    Python::attach(|py| {
        import(py);
        let qiskit_mod = py.import("qiskit").expect("Qiskit import worked already");
        let qiskit_version = qiskit_mod
            .getattr("__version__")
            .expect("Failed to get qiskit.__version__")
            .extract::<String>()
            .expect("Failed to cast version to string");
        let qiskit_version_int = version_string_to_hex(&qiskit_version);
        unsafe {
            assert_eq!(ffi::qk_api_version(), qiskit_version_int);
            let qc = ffi::qk_circuit_new(3, 1);
            assert!(!qc.is_null());
            assert_eq!(ffi::qk_circuit_num_qubits(qc), 3);
            assert_eq!(ffi::qk_circuit_num_clbits(qc), 1);
            ffi::qk_circuit_free(qc);
        }
    });
}

#[test]
fn test_transpile_capsule() {
    Python::attach(|py| {
        import(py);
        unsafe {
            let target = ffi::qk_target_new(5);
            assert!(!target.is_null());
            assert_eq!(ffi::qk_target_num_qubits(target), 5);
            ffi::qk_target_free(target);
        }
    });
}

#[test]
fn test_quantum_info_capsule() {
    Python::attach(|py| {
        import(py);
        unsafe {
            let obs = ffi::qk_obs_identity(4);
            assert!(!obs.is_null());
            assert_eq!(ffi::qk_obs_num_qubits(obs), 4);
            assert_eq!(ffi::qk_obs_num_terms(obs), 1);
            ffi::qk_obs_free(obs);
        }
    });
}
