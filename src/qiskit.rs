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

impl From<qiskit_sys::QkExitCode> for QiskitError {
    fn from(value: qiskit_sys::QkExitCode) -> Self {
        match value {
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
            qiskit_sys::QkExitCode_QkExitCode_TargetInvalidInstKey => {
                QiskitError::TargetInvalidInstKey
            }
            _ => panic!("Invalid option for QiskitError"),
        }
    }
}
