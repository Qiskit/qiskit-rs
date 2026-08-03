# qiskit-sys

qiskit-sys is a *-sys package within qiskit-rs. It mainly provides the following functionality:
* It links to system library `libqiskit` (see [Install the Qiskit C API](https://quantum.cloud.ibm.com/docs/en/guides/install-c-api)).
* It provides declarations for types and functions in `bindings.rs`, but does **not** add higher level abstractions.

qiskit-rs implements safe, high-level abstractions on top of qiskit-sys.

For more details about *-sys packages, see [*-sys Packages](https://doc.rust-lang.org/cargo/reference/build-scripts.html#-sys-packages).

## Installation Methods (Optional)

There are two installation methods, please select one. By default,
the `Clone` method is used.

* Clone (no path specified)

  Automatically clones and builds the qiskit c api from source. WARNING, cloning and building from source is very slow. 
  
  To specify this method, set the `QISKIT_CEXT_INSTALL_METHOD` to `clone`:
  
  ```bash
  export QISKIT_CEXT_INSTALL_METHOD="clone"
  ```

* Path (Manually specified path)

  Uses the qiskit c api binary from a path. 
  
  To specify this method, set the `QISKIT_CEXT_INSTALL_METHOD` to `path` and set `QISKIT_CEXT_PATH` to the root of the qiskit repository where the binary is installed:
  
  ```bash
  export QISKIT_CEXT_INSTALL_METHOD="path"
  export QISKIT_CEXT_PATH="path/to/qiskit-cext-dir"
  ```

## Installation

To build qiskit-sys,

1. (Optional) specify an [installation method](#installation-method-optional)
2. Build the qiskit-sys crate. From the root of qiskit-rs, run:

```bash
cargo build --package qiskit-sys
```

If the build command completes successfully, a `bindings.rs` file will be created in `$OUT_DIR/bindings.rs`, where `OUT_DIR` is an [environment variable set by cargo](https://doc.rust-lang.org/cargo/reference/environment-variables.html#environment-variables-cargo-sets-for-crates).
