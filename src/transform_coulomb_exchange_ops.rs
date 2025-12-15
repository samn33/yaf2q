//! Module for transformation from the coulomb/exchange operators to a Qubit Operator

use crate::binary_matrix::BinaryMatrix;
use crate::fermion_operator::CoulombExchangeOperators;
use crate::qubit_operator::{PauliOperator, PauliProduct, QubitOperator};

/// Transform coulomb/exchange operators to QubitOperator
///
/// # Arguments
///
/// * `em` - encoding matrix
/// * `ops` - coulomb/exchange operators
///
#[allow(unused_variables, unused_mut)]
pub fn transform_coulomb_exchange_ops(
    em: &BinaryMatrix,
    ops: &CoulombExchangeOperators,
) -> Result<QubitOperator, String> {
    let num_qubits = em.size_of_row();
    let mut qo = QubitOperator::new(num_qubits);

    let em_inv = match em.inverse() {
        Ok(m) => m,
        _ => {
            return Err("Inverse matrix calculation failed.".to_string());
        }
    };

    let mut e_i = BinaryMatrix::zero(1, num_qubits).unwrap();
    let mut e_j = BinaryMatrix::zero(1, num_qubits).unwrap();

    for ((i, j), coef) in ops.iter() {
        for k in 0..num_qubits {
            e_i.elements[0][k] = 0;
            e_j.elements[0][k] = 0;
        }
        e_i.elements[0][*i] = 1;
        e_j.elements[0][*j] = 1;
        e_i = e_i.mul(&em_inv).unwrap();
        e_j = e_j.mul(&em_inv).unwrap();

        // qo_i
        let mut qo_i = QubitOperator::new(num_qubits);
        let mut pp = PauliProduct::new();
        for n in 0..num_qubits {
            match e_i.elements[0][n] {
                0 => {
                    let _ = pp.insert(&PauliOperator::I);
                }
                1 => {
                    let _ = pp.insert(&PauliOperator::Z(n));
                }
                _ => {
                    return Err("Invalid elements is found.".to_string());
                }
            };
        }
	
        let _ = qo_i.add_term(&pp, (*coef * (-0.5)).into());
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        let _ = qo_i.add_term(&pp, (*coef * 0.5).into());

        // qo_j
        let mut qo_j = QubitOperator::new(num_qubits);
        let mut pp = PauliProduct::new();
        for n in 0..num_qubits {
            match e_j.elements[0][n] {
                0 => {
                    let _ = pp.insert(&PauliOperator::I);
                }
                1 => {
                    let _ = pp.insert(&PauliOperator::Z(n));
                }
                _ => {
                    return Err("Invalid elements is found.".to_string());
                }
            };
        }

        let _ = qo_j.add_term(&pp, (-0.5).into());
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        let _ = qo_j.add_term(&pp, (0.5).into());

        let _ = qo_j.mul(&qo_i);
        let _ = qo.add(&qo_j);
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
    use crate::fermion_operator::{CrAnOperator, FermionOperator};
    use float_cmp::approx_eq;

    #[test]
    fn transform_coulomb_exchange_ops_success_0() {
        // Fermion Operator
	//(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
	//(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            0.5,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(2),
            CrAnOperator::Annihilation(0),
            0.5,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(3),
            CrAnOperator::Annihilation(0),
            0.5,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(2),
            CrAnOperator::Annihilation(1),
            0.5,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(3),
            CrAnOperator::Annihilation(1),
            0.5,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(3),
            CrAnOperator::Annihilation(2),
            0.5,
        );

        // Jordan-Wigner
        let em = em_jordan_wigner();
        let coulomb_exchange_ops = fo.coulomb_exchange_operators().unwrap();
        let qo = transform_coulomb_exchange_ops(&em, &coulomb_exchange_ops).unwrap();

        assert_eq!(qo.len(), 11);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.75));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.375));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.375));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.375));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.375));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        // Bravyi-Kitaev
        let em = em_bravyi_kitaev();
        let coulomb_exchange_ops = fo.coulomb_exchange_operators().unwrap();
        let qo = transform_coulomb_exchange_ops(&em, &coulomb_exchange_ops).unwrap();

        assert_eq!(qo.len(), 11);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.75));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.375));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.375));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.375));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.375));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
    }

    #[test]
    fn transform_coulomb_exchange_ops_success_1() {
        // Fermion Operator
	//(0,1),(1,2),(2,3)
	//(0.5, 0.2, 0.1)
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            0.5,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(2),
            CrAnOperator::Annihilation(1),
            0.2,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(3),
            CrAnOperator::Annihilation(2),
            0.1,
        );

        // Jordan-Wigner
        let em = em_jordan_wigner();
        let coulomb_exchange_ops = fo.coulomb_exchange_operators().unwrap();
        let qo = transform_coulomb_exchange_ops(&em, &coulomb_exchange_ops).unwrap();

        assert_eq!(qo.len(), 8);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.19999999999999998));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.175));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.05));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.07500000000000001));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.025));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.025));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
    }
}

//Result of Openfermion

// [0]
//fermion hamiltonian
//0.5 [0^ 1^ 1 0] +
//0.5 [0^ 2^ 2 0] +
//0.5 [0^ 3^ 3 0] +
//0.5 [1^ 2^ 2 1] +
//0.5 [1^ 3^ 3 1] +
//0.5 [2^ 3^ 3 2]
//jw
//(0.75+0j) [] +
//(-0.375+0j) [Z0] +
//(0.125+0j) [Z0 Z1] +
//(0.125+0j) [Z0 Z2] +
//(0.125+0j) [Z0 Z3] +
//(-0.375+0j) [Z1] +
//(0.125+0j) [Z1 Z2] +
//(0.125+0j) [Z1 Z3] +
//(-0.375+0j) [Z2] +
//(0.125+0j) [Z2 Z3] +
//(-0.375+0j) [Z3]
//bk
//(0.75+0j) [] +
//(-0.375+0j) [Z0] +
//(-0.375+0j) [Z0 Z1] +
//(0.125+0j) [Z0 Z1 Z2] +
//(0.125+0j) [Z0 Z1 Z2 Z3] +
//(0.125+0j) [Z0 Z2] +
//(0.125+0j) [Z0 Z2 Z3] +
//(0.125+0j) [Z1] +
//(-0.375+0j) [Z1 Z2 Z3] +
//(0.125+0j) [Z1 Z3] +
//(-0.375+0j) [Z2]

// [1]
//fermion hamiltonian
//0.5 [0^ 1^ 1 0] +
//0.2 [1^ 2^ 2 1] +
//0.1 [2^ 3^ 3 2]
//jw
//(0.19999999999999998+0j) [] +
//(-0.125+0j) [Z0] +
//(0.125+0j) [Z0 Z1] +
//(-0.175+0j) [Z1] +
//(0.05+0j) [Z1 Z2] +
//(-0.07500000000000001+0j) [Z2] +
//(0.025+0j) [Z2 Z3] +
//(-0.025+0j) [Z3]
//bk
//(0.19999999999999998+0j) [] +
//(-0.125+0j) [Z0] +
//(-0.175+0j) [Z0 Z1] +
//(0.05+0j) [Z0 Z1 Z2] +
//(0.125+0j) [Z1] +
//(-0.025+0j) [Z1 Z2 Z3] +
//(0.025+0j) [Z1 Z3] +
//(-0.07500000000000001+0j) [Z2]
