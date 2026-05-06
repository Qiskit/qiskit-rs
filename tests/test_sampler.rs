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

use qiskit_rs::{QuantumCircuit, StatevectorSampler};

#[test]
fn test_sampler_bell_state() {
    let mut qc = QuantumCircuit::new(2, 2);
    qc.h(0);
    qc.cx(0, 1);
    qc.measure(0, 0);
    qc.measure(1, 1);

    let sampler = StatevectorSampler::with_seed(42);
    let result = sampler.run(&mut qc, 4096);
    let counts = result.counts();

    // Bell state produces only |00⟩ and |11⟩
    for key in counts.keys() {
        assert!(key == "00" || key == "11", "Unexpected: {}", key);
    }
    assert!(counts.len() == 2);
}

#[test]
fn test_sampler_x_deterministic() {
    let mut qc = QuantumCircuit::new(1, 1);
    qc.x(0);
    qc.measure(0, 0);

    let sampler = StatevectorSampler::with_seed(1);
    let result = sampler.run(&mut qc, 100);

    assert_eq!(*result.counts().get("1").unwrap(), 100);
}

#[test]
fn test_sampler_probabilities() {
    let mut qc = QuantumCircuit::new(1, 1);
    qc.h(0);
    qc.measure(0, 0);

    let sampler = StatevectorSampler::with_seed(77);
    let result = sampler.run(&mut qc, 10000);
    let probs = result.probabilities();

    let p0 = *probs.get("0").unwrap_or(&0.0);
    let p1 = *probs.get("1").unwrap_or(&0.0);
    assert!((p0 - 0.5).abs() < 0.05, "p(0) = {} too far from 0.5", p0);
    assert!((p1 - 0.5).abs() < 0.05, "p(1) = {} too far from 0.5", p1);
}

#[test]
fn test_sampler_ghz_10qubit() {
    let n = 10;
    let mut qc = QuantumCircuit::new(n, n);
    qc.h(0);
    for i in 0..(n - 1) {
        qc.cx(i, i + 1);
    }
    for i in 0..n {
        qc.measure(i, i);
    }

    let sampler = StatevectorSampler::with_seed(123);
    let result = sampler.run(&mut qc, 2048);
    let counts = result.counts();

    let zeros: String = (0..n).map(|_| '0').collect();
    let ones: String = (0..n).map(|_| '1').collect();

    for key in counts.keys() {
        assert!(
            *key == zeros || *key == ones,
            "GHZ-{}: unexpected outcome '{}'",
            n,
            key
        );
    }
}

#[test]
fn test_sampler_rz_phase() {
    // RZ should not change measurement probabilities on |0⟩
    let mut qc = QuantumCircuit::new(1, 1);
    qc.rz(1.234, 0);
    qc.measure(0, 0);

    let sampler = StatevectorSampler::with_seed(9);
    let result = sampler.run(&mut qc, 100);
    // |0⟩ is eigenstate of RZ, so always measures 0
    assert_eq!(*result.counts().get("0").unwrap(), 100);
}
