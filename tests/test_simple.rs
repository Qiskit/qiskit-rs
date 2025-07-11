use qiskit_rs::qk_circuit_new;

#[test]
fn test_new_circuit() {
    let _circuit = unsafe { qk_circuit_new(2, 2) };
}