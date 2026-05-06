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

//! Statevector-based quantum circuit sampler.
//!
//! Provides a local simulator backend that executes [`QuantumCircuit`]s via
//! full statevector simulation and returns shot-based measurement results.
//!
//! # Example
//!
//! ```
//! use qiskit_rs::{QuantumCircuit, StatevectorSampler, SamplerResult};
//!
//! let mut qc = QuantumCircuit::new(2, 2);
//! qc.h(0);
//! qc.cx(0, 1);
//! qc.measure(0, 0);
//! qc.measure(1, 1);
//!
//! let sampler = StatevectorSampler::new();
//! let result = sampler.run(&mut qc, 1024);
//! let counts = result.counts();
//! assert!(counts.contains_key("00") || counts.contains_key("11"));
//! ```

use crate::QuantumCircuit;
use std::collections::HashMap;

// ============================================================================
// Complex number
// ============================================================================

/// Minimal complex number for internal simulation use.
#[derive(Clone, Copy, Debug)]
struct Cx {
    re: f64,
    im: f64,
}

impl Cx {
    const ZERO: Cx = Cx { re: 0.0, im: 0.0 };
    const ONE: Cx = Cx { re: 1.0, im: 0.0 };
    const I: Cx = Cx { re: 0.0, im: 1.0 };

    const fn new(re: f64, im: f64) -> Self {
        Cx { re, im }
    }

    fn from_polar(r: f64, theta: f64) -> Self {
        Cx {
            re: r * theta.cos(),
            im: r * theta.sin(),
        }
    }

    fn norm_sq(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
}

impl std::ops::Add for Cx {
    type Output = Cx;
    fn add(self, rhs: Cx) -> Cx {
        Cx {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl std::ops::Mul for Cx {
    type Output = Cx;
    fn mul(self, rhs: Cx) -> Cx {
        Cx {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl std::ops::Mul<f64> for Cx {
    type Output = Cx;
    fn mul(self, rhs: f64) -> Cx {
        Cx {
            re: self.re * rhs,
            im: self.im * rhs,
        }
    }
}

impl std::ops::Neg for Cx {
    type Output = Cx;
    fn neg(self) -> Cx {
        Cx {
            re: -self.re,
            im: -self.im,
        }
    }
}

// ============================================================================
// Gate matrices
// ============================================================================

/// Single-qubit gate matrix stored row-major: [a00, a01, a10, a11].
type Gate1Q = [Cx; 4];

/// Two-qubit gate matrix stored row-major (4×4 = 16 entries).
/// Basis order relative to (q1, q2): |00⟩, |01⟩, |10⟩, |11⟩.
type Gate2Q = [Cx; 16];

const O: Cx = Cx::ZERO;
const I1: Cx = Cx::ONE;
const IM: Cx = Cx::I;

fn gate_id() -> Gate1Q {
    [I1, O, O, I1]
}

fn gate_x() -> Gate1Q {
    [O, I1, I1, O]
}

fn gate_y() -> Gate1Q {
    [O, -IM, IM, O]
}

fn gate_z() -> Gate1Q {
    [I1, O, O, -I1]
}

fn gate_h() -> Gate1Q {
    let s = Cx::new(std::f64::consts::FRAC_1_SQRT_2, 0.0);
    [s, s, s, -s]
}

fn gate_s() -> Gate1Q {
    [I1, O, O, IM]
}

fn gate_sdg() -> Gate1Q {
    [I1, O, O, -IM]
}

fn gate_t() -> Gate1Q {
    [I1, O, O, Cx::from_polar(1.0, std::f64::consts::FRAC_PI_4)]
}

fn gate_tdg() -> Gate1Q {
    [I1, O, O, Cx::from_polar(1.0, -std::f64::consts::FRAC_PI_4)]
}

fn gate_sx() -> Gate1Q {
    let a = Cx::new(0.5, 0.5);
    let b = Cx::new(0.5, -0.5);
    [a, b, b, a]
}

fn gate_sxdg() -> Gate1Q {
    let a = Cx::new(0.5, -0.5);
    let b = Cx::new(0.5, 0.5);
    [a, b, b, a]
}

fn gate_p(lam: f64) -> Gate1Q {
    [I1, O, O, Cx::from_polar(1.0, lam)]
}

fn gate_rx(theta: f64) -> Gate1Q {
    let c = Cx::new((theta / 2.0).cos(), 0.0);
    let s = Cx::new(0.0, -(theta / 2.0).sin());
    [c, s, s, c]
}

fn gate_ry(theta: f64) -> Gate1Q {
    let c = Cx::new((theta / 2.0).cos(), 0.0);
    let sn = Cx::new(-(theta / 2.0).sin(), 0.0);
    let sp = Cx::new((theta / 2.0).sin(), 0.0);
    [c, sn, sp, c]
}

fn gate_rz(phi: f64) -> Gate1Q {
    [
        Cx::from_polar(1.0, -phi / 2.0),
        O,
        O,
        Cx::from_polar(1.0, phi / 2.0),
    ]
}

fn gate_u(theta: f64, phi: f64, lam: f64) -> Gate1Q {
    let c = (theta / 2.0).cos();
    let s = (theta / 2.0).sin();
    [
        Cx::new(c, 0.0),
        -Cx::from_polar(s, lam),
        Cx::from_polar(s, phi),
        Cx::from_polar(c, phi + lam),
    ]
}

fn gate_r(theta: f64, phi: f64) -> Gate1Q {
    let c = Cx::new((theta / 2.0).cos(), 0.0);
    let s = (theta / 2.0).sin();
    [
        c,
        Cx::new(0.0, -1.0) * Cx::from_polar(s, -phi),
        Cx::new(0.0, -1.0) * Cx::from_polar(s, phi),
        c,
    ]
}

// --- Two-qubit gates ---

fn gate_cx() -> Gate2Q {
    //  |00⟩→|00⟩, |01⟩→|01⟩, |10⟩→|11⟩, |11⟩→|10⟩
    [
        I1, O, O, O, // row 0
        O, I1, O, O, // row 1
        O, O, O, I1, // row 2
        O, O, I1, O, // row 3
    ]
}

fn gate_iswap() -> Gate2Q {
    [
        I1, O, O, O, O, O, IM, O, O, IM, O, O, O, O, O, I1,
    ]
}

fn gate_ecr() -> Gate2Q {
    let s = Cx::new(std::f64::consts::FRAC_1_SQRT_2, 0.0);
    let si = Cx::new(0.0, std::f64::consts::FRAC_1_SQRT_2);
    [
        O, O, s, si, O, O, si, s, s, -si, O, O, -si, s, O, O,
    ]
}

fn gate_rxx(theta: f64) -> Gate2Q {
    let c = Cx::new((theta / 2.0).cos(), 0.0);
    let s = Cx::new(0.0, -(theta / 2.0).sin());
    [
        c, O, O, s, O, c, s, O, O, s, c, O, s, O, O, c,
    ]
}

fn gate_ryy(theta: f64) -> Gate2Q {
    let c = Cx::new((theta / 2.0).cos(), 0.0);
    let s = Cx::new(0.0, (theta / 2.0).sin());
    [
        c, O, O, s, O, c, -s, O, O, -s, c, O, s, O, O, c,
    ]
}

fn gate_rzz(theta: f64) -> Gate2Q {
    let p = Cx::from_polar(1.0, theta / 2.0);
    let m = Cx::from_polar(1.0, -theta / 2.0);
    [
        m, O, O, O, O, p, O, O, O, O, p, O, O, O, O, m,
    ]
}

fn gate_rzx(theta: f64) -> Gate2Q {
    let c = Cx::new((theta / 2.0).cos(), 0.0);
    let s = Cx::new(0.0, (theta / 2.0).sin());
    [
        c, O, -s, O, O, c, O, s, -s, O, c, O, O, s, O, c,
    ]
}

// ============================================================================
// Statevector
// ============================================================================

/// Internal quantum state representation as a dense complex vector of length 2^n.
struct Statevector {
    num_qubits: u32,
    state: Vec<Cx>,
}

impl Statevector {
    /// Initialise to |00…0⟩.
    fn new(num_qubits: u32) -> Self {
        let size = 1usize << num_qubits;
        let mut state = vec![Cx::ZERO; size];
        state[0] = Cx::ONE;
        Statevector { num_qubits, state }
    }

    /// Apply a single-qubit gate to the given qubit.
    fn apply_1q(&mut self, gate: Gate1Q, qubit: u32) {
        let size = 1usize << self.num_qubits;
        let bit = 1usize << qubit;
        for i in 0..size {
            if i & bit == 0 {
                let j = i | bit;
                let a = self.state[i];
                let b = self.state[j];
                self.state[i] = gate[0] * a + gate[1] * b;
                self.state[j] = gate[2] * a + gate[3] * b;
            }
        }
    }

    /// Apply a two-qubit gate. Basis order: |q1 q2⟩.
    fn apply_2q(&mut self, gate: Gate2Q, q1: u32, q2: u32) {
        let size = 1usize << self.num_qubits;
        let b1 = 1usize << q1;
        let b2 = 1usize << q2;
        for i in 0..size {
            if i & b1 == 0 && i & b2 == 0 {
                let idx = [i, i | b2, i | b1, i | b1 | b2];
                let s = [
                    self.state[idx[0]],
                    self.state[idx[1]],
                    self.state[idx[2]],
                    self.state[idx[3]],
                ];
                for row in 0..4 {
                    let r = row * 4;
                    self.state[idx[row]] =
                        gate[r] * s[0] + gate[r + 1] * s[1] + gate[r + 2] * s[2] + gate[r + 3] * s[3];
                }
            }
        }
    }

    /// Apply a general n-qubit gate (n ≤ 4). Used for RCCX / RC3X.
    fn apply_nq(&mut self, gate: &[Cx], qubits: &[u32]) {
        let k = qubits.len();
        let gs = 1usize << k;
        let size = 1usize << self.num_qubits;
        for i in 0..size {
            if qubits.iter().any(|&q| i & (1 << q) != 0) {
                continue;
            }
            let indices: Vec<usize> = (0..gs)
                .map(|g| {
                    let mut idx = i;
                    for (bit, &q) in qubits.iter().enumerate() {
                        if g & (1 << bit) != 0 {
                            idx |= 1 << q;
                        }
                    }
                    idx
                })
                .collect();
            let old: Vec<Cx> = indices.iter().map(|&idx| self.state[idx]).collect();
            for row in 0..gs {
                let mut sum = Cx::ZERO;
                for col in 0..gs {
                    sum = sum + gate[row * gs + col] * old[col];
                }
                self.state[indices[row]] = sum;
            }
        }
    }

    /// Return Born-rule probabilities for each computational basis state.
    fn probabilities(&self) -> Vec<f64> {
        self.state.iter().map(|c| c.norm_sq()).collect()
    }
}

// ============================================================================
// RCCX / RC3X gate matrices
// ============================================================================

/// Build the RCCX (simplified Toffoli) 8×8 matrix.
/// Acts as CCX but with relative phases on states where both controls are not 1.
fn gate_rccx() -> Vec<Cx> {
    // Decomposition-equivalent unitary for the simplified Toffoli.
    // |110⟩ ↔ |111⟩ swap; |101⟩ picks up a phase of -i; |100⟩ picks up +i.
    let mut m = vec![Cx::ZERO; 64];
    // Diagonal entries (identity on most states)
    m[0 * 8 + 0] = I1; // |000⟩
    m[1 * 8 + 1] = I1; // |001⟩
    m[2 * 8 + 2] = I1; // |010⟩
    m[3 * 8 + 3] = I1; // |011⟩
    m[4 * 8 + 4] = I1; // |100⟩
    m[5 * 8 + 5] = -IM; // |101⟩ phase
    // Swap |110⟩ ↔ |111⟩
    m[6 * 8 + 7] = I1; // |110⟩ → |111⟩
    m[7 * 8 + 6] = I1; // |111⟩ → |110⟩
    m
}

/// Build the RC3X (simplified 3-controlled Toffoli) 16×16 matrix.
fn gate_rc3x() -> Vec<Cx> {
    let mut m = vec![Cx::ZERO; 256];
    // Identity on most basis states
    for i in 0..16 {
        m[i * 16 + i] = I1;
    }
    // Clear the entries we override
    m[13 * 16 + 13] = Cx::ZERO;
    m[14 * 16 + 14] = Cx::ZERO;
    m[15 * 16 + 15] = Cx::ZERO;
    // |1101⟩ (13) picks up phase -i
    m[13 * 16 + 13] = -IM;
    // Swap |1110⟩ (14) ↔ |1111⟩ (15)
    m[14 * 16 + 15] = I1;
    m[15 * 16 + 14] = I1;
    m
}

// ============================================================================
// SamplerResult
// ============================================================================

/// The result of sampling a quantum circuit.
#[derive(Debug, Clone)]
pub struct SamplerResult {
    counts: HashMap<String, u32>,
}

impl SamplerResult {
    /// Return the measurement counts.
    ///
    /// Keys are bitstrings with the most-significant classical bit on the left.
    pub fn counts(&self) -> &HashMap<String, u32> {
        &self.counts
    }

    /// Return the measurement counts as probabilities (normalised to 1).
    pub fn probabilities(&self) -> HashMap<String, f64> {
        let total: u32 = self.counts.values().sum();
        self.counts
            .iter()
            .map(|(k, &v)| (k.clone(), v as f64 / total as f64))
            .collect()
    }
}

// ============================================================================
// StatevectorSampler
// ============================================================================

/// A statevector-based backend that simulates quantum circuits locally.
///
/// Performs full statevector simulation, computing Born-rule probabilities and
/// sampling measurement outcomes for a given number of shots.
///
/// # Example
///
/// ```
/// use qiskit_rs::{QuantumCircuit, StatevectorSampler};
///
/// let mut qc = QuantumCircuit::new(2, 2);
/// qc.h(0);
/// qc.cx(0, 1);
/// qc.measure(0, 0);
/// qc.measure(1, 1);
///
/// let sampler = StatevectorSampler::new();
/// let result = sampler.run(&mut qc, 1024);
/// println!("{:?}", result.counts());
/// ```
pub struct StatevectorSampler {
    seed: Option<u64>,
}

impl StatevectorSampler {
    /// Create a new sampler with a random seed.
    pub fn new() -> Self {
        StatevectorSampler { seed: None }
    }

    /// Create a new sampler with a fixed seed for reproducibility.
    pub fn with_seed(seed: u64) -> Self {
        StatevectorSampler { seed: Some(seed) }
    }

    /// Execute the circuit and return shot-based measurement results.
    ///
    /// # Panics
    ///
    /// Panics if `num_qubits > 24` (would require > 256 MiB of memory).
    pub fn run(&self, circuit: &mut QuantumCircuit, shots: u32) -> SamplerResult {
        let n_qubits = circuit.num_qubits();
        let n_clbits = circuit.num_clbits();
        assert!(
            n_qubits <= 24,
            "StatevectorSampler supports at most 24 qubits (requested {})",
            n_qubits
        );

        let mut sv = Statevector::new(n_qubits);
        let mut qubit_to_clbit: HashMap<u32, u32> = HashMap::new();

        // Walk circuit instructions
        for inst in circuit.instructions() {
            match inst.name {
                // --- single-qubit gates ---
                "id" => sv.apply_1q(gate_id(), inst.qubits[0]),
                "x" => sv.apply_1q(gate_x(), inst.qubits[0]),
                "y" => sv.apply_1q(gate_y(), inst.qubits[0]),
                "z" => sv.apply_1q(gate_z(), inst.qubits[0]),
                "h" => sv.apply_1q(gate_h(), inst.qubits[0]),
                "s" => sv.apply_1q(gate_s(), inst.qubits[0]),
                "sdg" => sv.apply_1q(gate_sdg(), inst.qubits[0]),
                "t" => sv.apply_1q(gate_t(), inst.qubits[0]),
                "tdg" => sv.apply_1q(gate_tdg(), inst.qubits[0]),
                "sx" => sv.apply_1q(gate_sx(), inst.qubits[0]),
                "sxdg" => sv.apply_1q(gate_sxdg(), inst.qubits[0]),
                "p" => {
                    let lam = unsafe { qiskit_sys::qk_param_as_real(inst.params[0]) };
                    sv.apply_1q(gate_p(lam), inst.qubits[0]);
                }
                "rx" => {
                    let theta = unsafe { qiskit_sys::qk_param_as_real(inst.params[0]) };
                    sv.apply_1q(gate_rx(theta), inst.qubits[0]);
                }
                "ry" => {
                    let theta = unsafe { qiskit_sys::qk_param_as_real(inst.params[0]) };
                    sv.apply_1q(gate_ry(theta), inst.qubits[0]);
                }
                "rz" => {
                    let phi = unsafe { qiskit_sys::qk_param_as_real(inst.params[0]) };
                    sv.apply_1q(gate_rz(phi), inst.qubits[0]);
                }
                "u" => {
                    let theta = unsafe { qiskit_sys::qk_param_as_real(inst.params[0]) };
                    let phi = unsafe { qiskit_sys::qk_param_as_real(inst.params[1]) };
                    let lam = unsafe { qiskit_sys::qk_param_as_real(inst.params[2]) };
                    sv.apply_1q(gate_u(theta, phi, lam), inst.qubits[0]);
                }
                "r" => {
                    let theta = unsafe { qiskit_sys::qk_param_as_real(inst.params[0]) };
                    let phi = unsafe { qiskit_sys::qk_param_as_real(inst.params[1]) };
                    sv.apply_1q(gate_r(theta, phi), inst.qubits[0]);
                }
                // --- two-qubit gates ---
                "cx" => sv.apply_2q(gate_cx(), inst.qubits[0], inst.qubits[1]),
                "dcx" => {
                    sv.apply_2q(gate_cx(), inst.qubits[0], inst.qubits[1]);
                    sv.apply_2q(gate_cx(), inst.qubits[1], inst.qubits[0]);
                }
                "ecr" => sv.apply_2q(gate_ecr(), inst.qubits[0], inst.qubits[1]),
                "iswap" => sv.apply_2q(gate_iswap(), inst.qubits[0], inst.qubits[1]),
                "rxx" => {
                    let theta = unsafe { qiskit_sys::qk_param_as_real(inst.params[0]) };
                    sv.apply_2q(gate_rxx(theta), inst.qubits[0], inst.qubits[1]);
                }
                "ryy" => {
                    let theta = unsafe { qiskit_sys::qk_param_as_real(inst.params[0]) };
                    sv.apply_2q(gate_ryy(theta), inst.qubits[0], inst.qubits[1]);
                }
                "rzz" => {
                    let theta = unsafe { qiskit_sys::qk_param_as_real(inst.params[0]) };
                    sv.apply_2q(gate_rzz(theta), inst.qubits[0], inst.qubits[1]);
                }
                "rzx" => {
                    let theta = unsafe { qiskit_sys::qk_param_as_real(inst.params[0]) };
                    sv.apply_2q(gate_rzx(theta), inst.qubits[0], inst.qubits[1]);
                }
                // --- three / four-qubit gates ---
                "rccx" => {
                    sv.apply_nq(&gate_rccx(), inst.qubits);
                }
                "rcccx" | "rc3x" => {
                    sv.apply_nq(&gate_rc3x(), inst.qubits);
                }
                // --- measurement ---
                "measure" => {
                    qubit_to_clbit.insert(inst.qubits[0], inst.clbits[0]);
                }
                other => panic!("StatevectorSampler: unsupported gate '{}'", other),
            }
        }

        // Sample from the probability distribution
        let probs = sv.probabilities();
        let counts = sample(&probs, shots, n_clbits, &qubit_to_clbit, self.seed);

        SamplerResult {
            counts,
        }
    }
}

impl Default for StatevectorSampler {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Sampling
// ============================================================================

/// Simple xorshift64 PRNG (avoids external dependency).
struct Xorshift64(u64);

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        Xorshift64(if seed == 0 { 0xdeadbeef } else { seed })
    }

    /// Return a pseudo-random f64 in [0, 1).
    fn next_f64(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 as f64) / (u64::MAX as f64)
    }
}

/// Sample `shots` bitstrings from the probability distribution, mapping
/// qubit indices to classical bit indices via `qubit_to_clbit`.
fn sample(
    probs: &[f64],
    shots: u32,
    n_clbits: u32,
    qubit_to_clbit: &HashMap<u32, u32>,
    seed: Option<u64>,
) -> HashMap<String, u32> {
    // Build cumulative distribution
    let mut cdf: Vec<f64> = Vec::with_capacity(probs.len());
    let mut acc = 0.0;
    for &p in probs {
        acc += p;
        cdf.push(acc);
    }
    // Normalise to guarantee last element ≈ 1.0
    if let Some(last) = cdf.last_mut() {
        *last = 1.0;
    }

    // Seed: use provided seed or derive one from time
    let s = seed.unwrap_or_else(|| {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos() as u64)
            .unwrap_or(42)
    });
    let mut rng = Xorshift64::new(s);

    let mut counts: HashMap<String, u32> = HashMap::new();
    for _ in 0..shots {
        let r = rng.next_f64();
        let state_idx = match cdf.binary_search_by(|v| v.partial_cmp(&r).unwrap()) {
            Ok(i) => i,
            Err(i) => i,
        };

        // Map quantum state bits → classical bitstring
        let mut clbits = vec!['0'; n_clbits as usize];
        for (&qubit, &clbit) in qubit_to_clbit {
            if state_idx & (1 << qubit) != 0 {
                // Classical bits are displayed MSB-first
                clbits[(n_clbits - 1 - clbit) as usize] = '1';
            }
        }
        let key: String = clbits.into_iter().collect();
        *counts.entry(key).or_insert(0) += 1;
    }

    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bell_state() {
        let mut qc = QuantumCircuit::new(2, 2);
        qc.h(0);
        qc.cx(0, 1);
        qc.measure(0, 0);
        qc.measure(1, 1);

        let sampler = StatevectorSampler::with_seed(12345);
        let result = sampler.run(&mut qc, 4096);
        let counts = result.counts();

        // Bell state: only "00" and "11" should appear
        for (key, _) in counts {
            assert!(
                key == "00" || key == "11",
                "Unexpected outcome: {}",
                key
            );
        }
        // Both should have roughly equal counts (~2048 each)
        let c00 = *counts.get("00").unwrap_or(&0);
        let c11 = *counts.get("11").unwrap_or(&0);
        assert!(c00 > 1500, "00 count too low: {}", c00);
        assert!(c11 > 1500, "11 count too low: {}", c11);
    }

    #[test]
    fn test_x_gate() {
        let mut qc = QuantumCircuit::new(1, 1);
        qc.x(0);
        qc.measure(0, 0);

        let sampler = StatevectorSampler::with_seed(99);
        let result = sampler.run(&mut qc, 100);
        let counts = result.counts();

        // X|0⟩ = |1⟩, so we should always measure "1"
        assert_eq!(*counts.get("1").unwrap_or(&0), 100);
        assert_eq!(*counts.get("0").unwrap_or(&0), 0);
    }

    #[test]
    fn test_hadamard_statistics() {
        let mut qc = QuantumCircuit::new(1, 1);
        qc.h(0);
        qc.measure(0, 0);

        let sampler = StatevectorSampler::with_seed(42);
        let result = sampler.run(&mut qc, 10000);
        let counts = result.counts();

        // Should be roughly 50/50
        let c0 = *counts.get("0").unwrap_or(&0);
        let c1 = *counts.get("1").unwrap_or(&0);
        assert!(c0 > 4000 && c0 < 6000, "0 count out of range: {}", c0);
        assert!(c1 > 4000 && c1 < 6000, "1 count out of range: {}", c1);
    }

    #[test]
    fn test_ghz_3qubit() {
        let mut qc = QuantumCircuit::new(3, 3);
        qc.h(0);
        qc.cx(0, 1);
        qc.cx(1, 2);
        qc.measure(0, 0);
        qc.measure(1, 1);
        qc.measure(2, 2);

        let sampler = StatevectorSampler::with_seed(7);
        let result = sampler.run(&mut qc, 4096);
        let counts = result.counts();

        // GHZ: only "000" and "111"
        for (key, _) in counts {
            assert!(
                key == "000" || key == "111",
                "Unexpected GHZ outcome: {}",
                key
            );
        }
    }

    #[test]
    fn test_deterministic_seed() {
        let mut qc = QuantumCircuit::new(2, 2);
        qc.h(0);
        qc.cx(0, 1);
        qc.measure(0, 0);
        qc.measure(1, 1);

        let r1 = StatevectorSampler::with_seed(555).run(&mut qc, 1000);
        let r2 = StatevectorSampler::with_seed(555).run(&mut qc, 1000);
        assert_eq!(r1.counts(), r2.counts());
    }
}
