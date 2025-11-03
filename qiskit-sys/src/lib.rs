#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_circuit() {
        unsafe {
            let qc = qk_circuit_new(0, 0);
            assert_eq!(0, qk_circuit_num_qubits(qc));
            assert_eq!(0, qk_circuit_num_clbits(qc));
            assert_eq!(0, qk_circuit_num_instructions(qc));
            let mut op_counts = qk_circuit_count_ops(qc);
            assert_eq!(0, op_counts.len);

            qk_circuit_free(qc);
            qk_opcounts_clear(&mut op_counts);
        }
    }
}
