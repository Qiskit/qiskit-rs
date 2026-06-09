# `quantum_info`

This module contains a collection of quantum operators and mathematical
tools that can be used to represent quantum states or processing the
output of quantum workflows.

``quantum_info`` mostly focuses on representing the current functionality
provided by the [Qiskit C API](https://docs.quantum.ibm.com/api/qiskit-c).
The main features represented are:

## `Observable`

This struct represents the functionality of ``QkObs`` representing an
observable over Pauli bases, storing its data in qubit-sparse format.

## `BitTerm`

An enumeration of the possible Pauli term operators.
