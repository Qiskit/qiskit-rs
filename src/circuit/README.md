# `circuit`

The main module for circuit representation in `qiskit-rs`.

The purpose of this module is to establish the designated functionality for
circuit creation and interaction provided by the [Qiskit C API](https://docs.quantum.ibm.com/api/qiskit-c)
while making it safer to run and imitating its behavior in Python.

Here's what the current model looks like:

## `QuantumCircuit`

As it stands, it represents the C native ``QkCircuit`` struct which contains a list of instructions
comprising the circuit instance. You may refer to the local documentation as to how to use it further.
