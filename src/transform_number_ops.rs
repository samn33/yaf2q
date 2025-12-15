//! Module for transformation from the number operators to a Qubit Operator

use crate::binary_matrix::BinaryMatrix;
use crate::fermion_operator::NumberOperators;
use crate::qubit_operator::{PauliOperator, PauliProduct, QubitOperator};

/// Transform number operators to Qubit Operator
///
/// # Arguments
///
/// * `em` - encoding matrix
/// * `ops` - number operators
///
#[allow(unused_variables, unused_mut)]
pub fn transform_number_ops(
    em: &BinaryMatrix,
    ops: &NumberOperators,
) -> Result<QubitOperator, String> {
    let num_qubits = em.size_of_row();
    let mut qo = QubitOperator::new(num_qubits);

    let em_inv = match em.inverse() {
        Ok(m) => m,
        _ => {
            return Err("Inverse matrix calculation failed.".to_string());
        }
    };

    let mut e = BinaryMatrix::zero(1, num_qubits).unwrap();
    for (idx, coef) in ops.iter() {
        for i in 0..num_qubits {
            e.elements[0][i] = 0;
        }
        e.elements[0][*idx] = 1;
        e = e.mul(&em_inv).unwrap();

        let mut pp = PauliProduct::new();
        for i in 0..num_qubits {
            match e.elements[0][i] {
                0 => {
                    let _ = pp.insert(&PauliOperator::I);
                }
                1 => {
                    let _ = pp.insert(&PauliOperator::Z(i));
                }
                _ => {
                    return Err("Invalid elements is found.".to_string());
                }
            };
        }
        let _ = qo.add_term(&pp, (*coef * (-0.5)).into());

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        let _ = qo.add_term(&pp, (*coef * 0.5).into());
    }

    Ok(qo)
}

mod common_test_data {
    use super::*;

    #[allow(dead_code)]
    pub fn em_jordan_wigner() -> BinaryMatrix {
        BinaryMatrix::any(&vec![
            vec![1, 0, 0, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![0, 0, 0, 1],
        ])
        .unwrap()
    }

    #[allow(dead_code)]
    pub fn em_parity() -> BinaryMatrix {
        BinaryMatrix::any(&vec![
            vec![1, 0, 0, 0],
            vec![1, 1, 0, 0],
            vec![1, 1, 1, 0],
            vec![1, 1, 1, 1],
        ])
        .unwrap()
    }

    #[allow(dead_code)]
    pub fn em_bravyi_kitaev() -> BinaryMatrix {
        BinaryMatrix::any(&vec![
            vec![1, 0, 0, 0],
            vec![1, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![1, 1, 1, 1],
        ])
        .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::common_test_data::*;
    use super::*;
    use crate::fermion_operator::CrAnOperator;
    use crate::fermion_operator::FermionOperator;
    use float_cmp::approx_eq;

    #[test]
    fn transform_number_ops_success_0() {
        // Jordan-Wigner
	// (0, 1, 2, 3)
	// (1.0, 2.0, 3.0, 4.0)
        let em = em_jordan_wigner();
        let mut fo = FermionOperator::new();
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(1),
            2.0,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(2),
            3.0,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(3),
            4.0,
        );
        let number_ops = fo.number_operators().unwrap();
        let qo = transform_number_ops(&em, &number_ops).unwrap();

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 5.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -1.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -1.5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -2.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        // Parity
        let em = em_parity();
        let mut fo = FermionOperator::new();
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(1),
            2.0,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(2),
            3.0,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(3),
            4.0,
        );
        let number_ops = fo.number_operators().unwrap();
        let qo = transform_number_ops(&em, &number_ops).unwrap();

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 5.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -1.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -1.5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -2.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        // Bravyi-Kitaev
        let em = em_bravyi_kitaev();
        let mut fo = FermionOperator::new();
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(1),
            2.0,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(2),
            3.0,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(3),
            4.0,
        );
        let number_ops = fo.number_operators().unwrap();
        let qo = transform_number_ops(&em, &number_ops).unwrap();

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 5.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -1.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -1.5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -2.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
    }
}

//Result of Openfermion

// [0]
//fermion hamiltonian
//1.0 [0^ 0] +
//2.0 [1^ 1] +
//3.0 [2^ 2] +
//4.0 [3^ 3]
//jw
//(5+0j) [] +
//(-0.5+0j) [Z0] +
//(-1+0j) [Z1] +
//(-1.5+0j) [Z2] +
//(-2+0j) [Z3]
//bk
//(5+0j) [] +
//(-0.5+0j) [Z0] +
//(-1+0j) [Z0 Z1] +
//(-2+0j) [Z1 Z2 Z3] +
//(-1.5+0j) [Z2]
