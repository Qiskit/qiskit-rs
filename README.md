# qiskit-rs

Rust bindings for qiskit


## Installation

A crate will be published soon. For now, qiskit-rs can be installed as follows:

```
cargo install --path=path/to/qiskit-rs
```

qiskit-rs needs the qiskit c api to function. There are two installation methods
built in to qiskit-rs:

- Path (Manually specified path): Uses qiskit c api binary or source from a path
    Set envvar `QISKIT_CEXT_PATH` to the path of the qiskit source directory
- Clone (no path specified): Automatically clones and builds the qiskit c api from source. 
    Set envvar `QISKIT_CEXT_CLONE=1` to use the clone method. WARNING, cloning and building from source is very slow.
