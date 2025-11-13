> [!WARNING]
> This code is experimental and its API may change without prior warning. Use it at your own risk.

# qiskit-rs

This library exposes the [C API for Qiskit](https://docs.quantum.ibm.com/api/qiskit-c) in Rust.

## Installation

```bash
cargo add --git https://github.com/Qiskit/qiskit-rs qiskit
```

## Example

Create a simple bell state circuit in qiskit-rs

```
use qiskit_rs::QuantumCircuit;

// Initialize a circuit with 2 quantum registers and 2 classical registers
let mut qc = QuantumCircuit::new(2, 2);

// Apply circuit operations
qc.h(0);
qc.cx(0, 1);

// Apply measurements
qc.measure(0, 0);
qc.measure(1, 1);
```

## License

[Apache License 2.0](https://github.com/Qiskit/qiskit/blob/main/LICENSE.txt)
