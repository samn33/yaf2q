//! Module for transformation from the constant value to a Qubit Operator

use crate::binary_matrix::BinaryMatrix;
use crate::qubit_operator::{PauliOperator, PauliProduct, QubitOperator};

/// Transform a constant value to QubitOperator
///
/// # Arguments
///
/// * `constant` - constant value
///
#[allow(unused_variables, unused_mut)]
pub fn transform_constant(em: &BinaryMatrix, constant: f64) -> Result<QubitOperator, String> {
    let num_qubits = em.size_of_row();
    let mut qo = QubitOperator::new(num_qubits);

    let mut pp = PauliProduct::new();
    let _ = pp.insert(&PauliOperator::I);
    let _ = qo.add_term(&pp, constant.into());

    Ok(qo)
}
