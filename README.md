# qiskit-rs

Rust bindings for qiskit


## Installation

qiskit-rs needs the qiskit c api to function. There are two installation methods
built in to qiskit-rs:

- Path (Manually specified path): Uses qiskit c api binary or source from a path
    - Set envvar `QISKIT_CEXT_INSTALL_METHOD=path`
    - Set envvar `QISKIT_CEXT_PATH` to the path of the qiskit source directory
- Clone (no path specified): Automatically clones and builds the qiskit c api from source. 
    - Set envvar `QISKIT_CEXT_INSTALL_METHOD=clone` to use the clone method. WARNING, cloning and building from source is very slow.

A crate will be published soon. For now, qiskit-rs can be installed as follows:

```
cargo install --path=path/to/qiskit-rs
```

## Getting Started

Create a simple GHZ state circuit with qiskit-rs:

```
use qiskit_rs::QuantumCircuit;

fn main() {
    let num_qubits = 10;
    let mut qc = QuantumCircuit::new(num_qubits, num_qubits);
    qc.h(0);
    for i in 0..(num_qubits - 1) {
        qc.cx(i, i + 1);
    }
    for i in 0..num_qubits {
        qc.measure(i, i);
    }
}
```
