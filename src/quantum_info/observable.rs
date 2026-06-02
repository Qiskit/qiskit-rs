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
use std::fmt;

use crate::Complex64;

/// An enum that represents single qubit alphabet terms
#[derive(Clone, Copy, PartialEq, Debug)]
#[repr(u8)]
pub enum BitTerm {
    /// Pauli X operator
    X = qiskit_sys::QkBitTerm_QkBitTerm_X,
    /// Pauli Y operator
    Y = qiskit_sys::QkBitTerm_QkBitTerm_Y,
    /// Pauli Z operator
    Z = qiskit_sys::QkBitTerm_QkBitTerm_Z,
    /// The projector ∣+⟩⟨+∣∣+⟩⟨+∣ to the positive X eigenstate
    Plus = qiskit_sys::QkBitTerm_QkBitTerm_Plus,
    /// The projector ∣−⟩⟨−∣∣−⟩⟨−∣ to the negative X eigenstate
    Minus = qiskit_sys::QkBitTerm_QkBitTerm_Minus,
    /// The projector ∣r⟩⟨r∣∣r⟩⟨r∣ to the positive Y eigenstate
    Right = qiskit_sys::QkBitTerm_QkBitTerm_Right,
    /// The projector ∣l⟩⟨l∣∣l⟩⟨l∣ to the negative Y eigenstate
    Left = qiskit_sys::QkBitTerm_QkBitTerm_Left,
    /// The projector ∣0⟩⟨0∣∣0⟩⟨0∣ to the positive Z eigenstate
    Zero = qiskit_sys::QkBitTerm_QkBitTerm_Zero,
    /// The projector ∣1⟩⟨1∣∣1⟩⟨1∣ to the negative Z eigenstate
    One = qiskit_sys::QkBitTerm_QkBitTerm_One,
}

impl TryFrom<u8> for BitTerm {
    type Error = ();
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            qiskit_sys::QkBitTerm_QkBitTerm_X => Ok(BitTerm::X),
            qiskit_sys::QkBitTerm_QkBitTerm_Y => Ok(BitTerm::Y),
            qiskit_sys::QkBitTerm_QkBitTerm_Z => Ok(BitTerm::Z),
            qiskit_sys::QkBitTerm_QkBitTerm_Plus => Ok(BitTerm::Plus),
            qiskit_sys::QkBitTerm_QkBitTerm_Minus => Ok(BitTerm::Minus),
            qiskit_sys::QkBitTerm_QkBitTerm_Right => Ok(BitTerm::Right),
            qiskit_sys::QkBitTerm_QkBitTerm_Left => Ok(BitTerm::Left),
            qiskit_sys::QkBitTerm_QkBitTerm_Zero => Ok(BitTerm::Zero),
            qiskit_sys::QkBitTerm_QkBitTerm_One => Ok(BitTerm::One),
            _ => Err(()),
        }
    }
}

/// A list of BitTerms. BitTermList implements a custom iterator.
pub struct BitTermList(Vec<BitTerm>);

impl BitTermList {
    /// Create a BitTermList
    pub fn new() -> BitTermList {
        BitTermList(Vec::new())
    }
    /// Add a BitTerm to a BitTermList
    pub fn add(&mut self, elem: BitTerm) {
        self.0.push(elem);
    }
    /// Get the length of a BitTermList
    pub fn len(self) -> usize {
        self.0.len()
    }
    /// Convert a BitTermList to a slice
    pub fn as_slice(&self) -> &[BitTerm] {
        self.0.as_slice()
    }
}

impl FromIterator<BitTerm> for BitTermList {
    fn from_iter<T: IntoIterator<Item = BitTerm>>(iter: T) -> Self {
        let mut bit_term_list = BitTermList::new();

        for i in iter {
            bit_term_list.add(i);
        }

        bit_term_list
    }
}

/// An observable over Pauli bases that stores its data in a qubit-sparse format.
pub struct Observable {
    observable: *mut qiskit_sys::QkObs,
}

impl Observable {
    /// Create a new observable
    ///
    /// # Example
    ///
    /// ```
    /// use qiskit_rs::{Observable, Complex64, BitTerm};
    //
    /// let num_qubits = 100;
    /// let coeffs = [Complex64{re: 1.0, im: -1.0}, Complex64{re: 1.0, im: -1.0}];
    /// let bits = [BitTerm::Zero, BitTerm::One, BitTerm::Plus, BitTerm::Minus];
    /// let indices = [0, 1, 98, 99];
    /// let boundaries = [0, 2, 4];
    ///
    /// let obs = Observable::new(num_qubits, &coeffs, &bits, &indices, &boundaries);
    /// ```
    pub fn new(
        num_qubits: u32,
        coeffs: &[Complex64],
        bit_terms: &[BitTerm],
        indices: &[u32],
        boundaries: &[usize],
    ) -> Observable {
        // Input shape checks (see https://quantum.cloud.ibm.com/docs/en/api/qiskit-c/qk-obs#representation)
        assert!(bit_terms.len() == indices.len());
        assert!(coeffs.len() + 1 == boundaries.len());

        let mut coeffs: Vec<qiskit_sys::QkComplex64> = Vec::from(coeffs);
        let mut bit_terms: Vec<u8> = bit_terms.iter().map(|bitterm| *(bitterm) as u8).collect();
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
    /// use qiskit_rs::{Observable, Complex64, BitTerm};
    ///
    /// let num_qubits = 100;
    /// let coeffs = [Complex64{re: 1.0, im: -1.0}, Complex64{re: 1.0, im: -1.0}];
    /// let bits = [BitTerm::Zero, BitTerm::One, BitTerm::Plus, BitTerm::Minus];
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
    /// use qiskit_rs::{Observable, Complex64, BitTerm};
    /// let num_qubits = 100;
    /// let coeffs = [Complex64{re: 1.0, im: -1.0}, Complex64{re: 1.0, im: -1.0}];
    /// let bits = [BitTerm::Zero, BitTerm::One, BitTerm::Plus, BitTerm::Minus];
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
    /// Get the list of boundaries of an observable
    ///
    /// ```
    /// use qiskit_rs::{Observable, Complex64, BitTerm};
    /// let num_qubits = 100;
    /// let coeffs = [Complex64{re: 1.0, im: -1.0}, Complex64{re: 1.0, im: -1.0}];
    /// let bits = [BitTerm::Zero, BitTerm::One, BitTerm::Plus, BitTerm::Minus];
    /// let indices = [0, 1, 98, 99];
    /// let boundaries = [0, 2, 4];
    /// let obs = Observable::new(num_qubits, &coeffs, &bits, &indices, &boundaries);
    /// for b in obs.boundaries() {
    ///     print!("Boundary: {}", b);
    /// }
    /// ```
    pub fn boundaries(&self) -> std::slice::Iter<'_, usize> {
        let num_boundaries: usize = unsafe { qiskit_sys::qk_obs_num_terms(self.observable) } + 1;
        let boundaries_ptr = unsafe { qiskit_sys::qk_obs_boundaries(self.observable) };
        let slice = unsafe { std::slice::from_raw_parts(boundaries_ptr, num_boundaries) };
        slice.iter()
    }
    /// Get the list of bitterms of an observable
    ///
    /// ```
    /// use qiskit_rs::{Observable, Complex64, BitTerm};
    /// let num_qubits = 100;
    /// let coeffs = [Complex64{re: 1.0, im: -1.0}, Complex64{re: 1.0, im: -1.0}];
    /// let bits = [BitTerm::Zero, BitTerm::One, BitTerm::Plus, BitTerm::Minus];
    /// let indices = [0, 1, 98, 99];
    /// let boundaries = [0, 2, 4];
    /// let obs = Observable::new(num_qubits, &coeffs, &bits, &indices, &boundaries);
    /// for b in obs.bit_terms().as_slice().iter() {
    ///     print!("Bit Terms: {:?}", b);
    /// }
    /// ```
    pub fn bit_terms(&self) -> BitTermList {
        let num_bit_terms: usize = unsafe { qiskit_sys::qk_obs_len(self.observable) };
        let bit_terms_ptr = unsafe { qiskit_sys::qk_obs_bit_terms(self.observable) };
        let slice: &[u8] = unsafe { std::slice::from_raw_parts(bit_terms_ptr, num_bit_terms) };
        let mut qk_bit_terms = BitTermList::new();
        for item in slice {
            let bitterm = BitTerm::try_from(*item).unwrap();
            qk_bit_terms.add(bitterm);
        }
        qk_bit_terms
    }
    /// Get the number of bit terms/indices in the observable.
    pub fn len(&self) -> usize {
        unsafe { qiskit_sys::qk_obs_len(self.observable) }
    }

    /// Compare two observables for equality.
    pub fn equal(&self, obs: &Observable) -> bool {
        unsafe { qiskit_sys::qk_obs_equal(self.observable, obs.observable) }
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
        let new = self.clone();
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
        let new = self.clone();
        unsafe { qiskit_sys::qk_obs_canonicalize(new.observable, tolerance) };
        new
    }
}

impl fmt::Display for Observable {
    /// Return a string representation of the observable
    ///
    /// # Example
    ///
    /// ```
    /// use qiskit_rs::Observable;
    ///
    /// let identity = Observable::identity(100);
    /// println!("{}", identity);
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let obs_str = unsafe { qiskit_sys::qk_obs_str(self.observable) };
        // Clone C string into String, which implements Drop
        let retval = String::from(unsafe { CStr::from_ptr(obs_str) }.to_str().unwrap());
        unsafe { qiskit_sys::qk_str_free(obs_str) };
        write!(f, "{}", retval)
    }
}

impl Clone for Observable {
    /// Copy the observable
    fn clone(&self) -> Observable {
        Observable {
            observable: unsafe { qiskit_sys::qk_obs_copy(self.observable) },
        }
    }
}

impl Drop for Observable {
    fn drop(&mut self) {
        unsafe { qiskit_sys::qk_obs_free(self.observable) };
    }
}

#[cfg(test)]
mod observable_tests {
    use crate::quantum_info::observable::BitTerm;

    use super::{Complex64, Observable};
    #[test]
    fn test_new_observable() {
        let num_qubits = 100;
        let coeffs = [
            Complex64 { re: 1.0, im: -1.0 },
            Complex64 { re: 1.0, im: -1.0 },
        ];
        let bits = [BitTerm::Zero, BitTerm::One, BitTerm::Plus, BitTerm::Minus];
        let indices = [0, 1, 98, 99];
        let boundaries = [0, 2, 4];

        let obs = Observable::new(num_qubits, &coeffs, &bits, &indices, &boundaries);
        assert_eq!(obs.num_terms(), 2);
        assert_eq!(obs.coeffs().len(), coeffs.len());
        for (i, coef) in obs.coeffs().enumerate() {
            assert_eq!(coef.re, coeffs[i].re);
            assert_eq!(coef.im, coeffs[i].im);
        }
        assert_eq!(obs.indices().len(), indices.len());
        for (i, idx) in obs.indices().enumerate() {
            assert_eq!(*idx, indices[i]);
        }
        assert_eq!(obs.boundaries().len(), boundaries.len());
        for (i, idx) in obs.boundaries().enumerate() {
            assert_eq!(*idx, boundaries[i]);
        }

        assert_eq!(obs.bit_terms().len(), bits.len());
        for (i, idx) in obs.bit_terms().as_slice().iter().enumerate() {
            assert_eq!(*idx, bits[i]);
        }

        assert_eq!(obs.coeffs().len() + 1, obs.boundaries().len());
        assert_eq!(obs.bit_terms().len(), obs.indices().len());

        let obs_b = Observable::new(num_qubits, &coeffs, &bits, &indices, &boundaries);
        assert!(obs.equal(&obs_b));
    }
}
