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

use qiskit_rs::quantum_info::BitTerm;
use qiskit_rs::{Complex64, Observable};

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
