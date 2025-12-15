//! Module for transformation from the double excitation operators to a Qubit Operator

use crate::binary_matrix::BinaryMatrix;
use crate::fermion_operator::DoubleExcitationOperators;
use crate::qubit_operator::{PauliOperator, PauliProduct, QubitOperator};

/// Transform double excitation operators to QubitOperator
///
/// # Arguments
///
/// * `em` - encoding Matrix
/// * `ops` - double excitation operators
///
#[allow(unused_variables, unused_mut)]
pub fn transform_double_excitation_ops(
    em: &BinaryMatrix,
    ops: &DoubleExcitationOperators,
) -> Result<QubitOperator, String> {
    let num_qubits = em.size_of_row();
    let mut qo = QubitOperator::new(num_qubits);

    let em_inv = match em.inverse() {
        Ok(m) => m,
        _ => {
            return Err("Inverse matrix calculation failed.".to_string());
        }
    };

    let mut e_i = BinaryMatrix::zero(num_qubits, 1).unwrap();
    let mut e_j = BinaryMatrix::zero(num_qubits, 1).unwrap();
    let mut e_k = BinaryMatrix::zero(num_qubits, 1).unwrap();
    let mut e_l = BinaryMatrix::zero(num_qubits, 1).unwrap();
    let mut one_i = BinaryMatrix::zero(num_qubits, 1).unwrap();
    let mut one_j = BinaryMatrix::zero(num_qubits, 1).unwrap();
    let mut one_k = BinaryMatrix::zero(num_qubits, 1).unwrap();
    let mut one_l = BinaryMatrix::zero(num_qubits, 1).unwrap();

    for ((i, j, k, l), coef) in ops.iter() {
	let a = match (i,j) {(i,j) if i < j => 0, _ => 1};
	let b = match (i,k) {(i,k) if i < k => 0, _ => 1};
	let c = match (l,j) {(l,j) if l < j => 0, _ => 1};
	let d = match (l,k) {(l,k) if l < k => 0, _ => 1};
	let e = match (l,i) {(l,i) if l < i => 0, _ => 1};
	let f = match (k,j) {(k,j) if k < j => 0, _ => 1};
	let mut factor = -1.0;
	if (a + b + c + d + e + f) % 2 == 0 {
	    factor = 1.0;
	}

        for n in 0..num_qubits {
            e_i.elements[n][0] = 0;
            e_j.elements[n][0] = 0;
            e_k.elements[n][0] = 0;
            e_l.elements[n][0] = 0;
            one_i.elements[n][0] = 0;
            one_j.elements[n][0] = 0;
            one_k.elements[n][0] = 0;
            one_l.elements[n][0] = 0;
        }
        e_i.elements[*i][0] = 1;
        e_j.elements[*j][0] = 1;
        e_k.elements[*k][0] = 1;
        e_l.elements[*l][0] = 1;

        // v_ijkl_t
        let mut e_ijkl = e_i.clone();
        e_ijkl = e_ijkl.add(&e_j).unwrap();
        e_ijkl = e_ijkl.add(&e_k).unwrap();
        e_ijkl = e_ijkl.add(&e_l).unwrap();
        let v_ijkl = em.mul(&e_ijkl).unwrap();
        let v_ijkl_t = v_ijkl.transpose().unwrap();

        // e_i_t
        let mut e_i_t = e_i.transpose().unwrap();
        e_i_t = e_i_t.mul(&em_inv).unwrap();

        // e_j_t
        let mut e_j_t = e_j.transpose().unwrap();
        e_j_t = e_j_t.mul(&em_inv).unwrap();

        // e_k_t
        let mut e_k_t = e_k.transpose().unwrap();
        e_k_t = e_k_t.mul(&em_inv).unwrap();

        // e_l_t
        let mut e_l_t = e_l.transpose().unwrap();
        e_l_t = e_l_t.mul(&em_inv).unwrap();

        // one_ijkl_t
        for n in 0..num_qubits {
            match i {
                i if *i != 0 && n <= *i - 1 => {
                    one_i.elements[n][0] = 1;
                }
                _ => {}
            }
            match j {
                j if *j != 0 && n <= *j - 1 => {
                    one_j.elements[n][0] = 1;
                }
                _ => {}
            }
            match k {
                k if *k != 0 && n <= *k - 1 => {
                    one_k.elements[n][0] = 1;
                }
                _ => {}
            }
            match l {
                l if *l != 0 && n <= *l - 1 => {
                    one_l.elements[n][0] = 1;
                }
                _ => {}
            }
        }
        let mut one_ijkl = one_i.clone();
        let one_ijkl = one_ijkl.add(&one_j).unwrap();
        let one_ijkl = one_ijkl.add(&one_k).unwrap();
        let one_ijkl = one_ijkl.add(&one_l).unwrap();
        let mut one_ijkl_t = one_ijkl.transpose().unwrap();
        one_ijkl_t = one_ijkl_t.mul(&em_inv).unwrap();

        // Qubit Operator: qo_l
        let mut qo_l = QubitOperator::new(num_qubits);
        let mut pp = PauliProduct::new();
        for n in 0..num_qubits {
            match v_ijkl_t.elements[0][n] {
                0 => {
                    let _ = pp.insert(&PauliOperator::I);
                }
                1 => {
                    let _ = pp.insert(&PauliOperator::X(n));
                }
                _ => {
                    return Err("Invalid elements is found.".to_string());
                }
            };
        }
        let _ = qo_l.add_term(&pp, (1.0).into());

        // Qubit Operator: qo_r
        let mut qo_r = QubitOperator::new(num_qubits);
        let mut pp = PauliProduct::new();
        for n in 0..num_qubits {
            match v_ijkl_t.elements[0][n] {
                0 => {
                    let _ = pp.insert(&PauliOperator::I);
                }
                1 => {
                    let _ = pp.insert(&PauliOperator::X(n));
                }
                _ => {
                    return Err("Invalid elements is found.".to_string());
                }
            };
        }
       let _ = qo_r.add_term(&pp, (1.0).into());

        // Qubit Operator: qo_m_0
        let mut qo_m_0 = QubitOperator::new(num_qubits);
        let mut pp = PauliProduct::new();
        for n in 0..num_qubits {
            match e_i_t.elements[0][n] {
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
        let _ = qo_m_0.add_term(&pp, (0.5).into());
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        let _ = qo_m_0.add_term(&pp, (0.5).into());

        // Qubit Operator: qo_m_1
        let mut qo_m_1 = QubitOperator::new(num_qubits);
        let mut pp = PauliProduct::new();
        for n in 0..num_qubits {
            match e_j_t.elements[0][n] {
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
        let _ = qo_m_1.add_term(&pp, (0.5).into());
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        let _ = qo_m_1.add_term(&pp, (0.5).into());

        // Qubit Operator: qo_m_2
        let mut qo_m_2 = QubitOperator::new(num_qubits);
        let mut pp = PauliProduct::new();
        for n in 0..num_qubits {
            match e_k_t.elements[0][n] {
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
        let _ = qo_m_2.add_term(&pp, (-0.5).into());
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        let _ = qo_m_2.add_term(&pp, (0.5).into());

        // Qubit Operator: qo_m_3
        let mut qo_m_3 = QubitOperator::new(num_qubits);
        let mut pp = PauliProduct::new();
        for n in 0..num_qubits {
            match e_l_t.elements[0][n] {
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
        let _ = qo_m_3.add_term(&pp, (-0.5).into());
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        let _ = qo_m_3.add_term(&pp, (0.5).into());

        // Qubit Operator: qo_m_4
        let mut qo_m_4 = QubitOperator::new(num_qubits);
        let mut pp = PauliProduct::new();
        for n in 0..num_qubits {
            match one_ijkl_t.elements[0][n] {
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
        let _ = qo_m_4.add_term(&pp, (*coef * factor).into());
	
        // Qubit Operator: qo_m = qo_m_4 * qo_m_3 * qo_m_2 * qo_m_1 * qo_m_0
        let mut qo_m = qo_m_4.clone();
        let _ = qo_m.mul(&qo_m_3);
        let _ = qo_m.mul(&qo_m_2);
        let _ = qo_m.mul(&qo_m_1);
        let _ = qo_m.mul(&qo_m_0);

        // Qubit Operator: qo_r = qo_r * qo_m
        let _ = qo_r.mul(&qo_m);

        // Qubit Operator: qo_m = qo_m * qo_l
        let _ = qo_m.mul(&qo_l);

        // Qubit Operator: qo = qo_r + qo_m
        let _ = qo.add(&qo_r);
        let _ = qo.add(&qo_m);
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
    pub fn em_jordan_wigner_6() -> BinaryMatrix {
        BinaryMatrix::any(&vec![
            vec![1, 0, 0, 0, 0, 0],
            vec![0, 1, 0, 0, 0, 0],
            vec![0, 0, 1, 0, 0, 0],
            vec![0, 0, 0, 1, 0, 0],
            vec![0, 0, 0, 0, 1, 0],
            vec![0, 0, 0, 0, 0, 1],
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
    fn transform_double_excitation_ops_success_0() {
        // Fermion Operator
        // (0,1,4,5),
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(4),
            CrAnOperator::Annihilation(5),
            0.5,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(5),
            CrAnOperator::Creation(4),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            0.5,
        );

        // Jordan-Wigner
        let em = em_jordan_wigner_6();
        let double_excitation_ops = fo.double_excitation_operators().unwrap();
        let qo = transform_double_excitation_ops(&em, &double_excitation_ops).unwrap();

        assert_eq!(qo.len(), 8);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::X(4));
        let _ = pp.insert(&PauliOperator::X(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Y(4));
        let _ = pp.insert(&PauliOperator::Y(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::X(4));
        let _ = pp.insert(&PauliOperator::Y(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Y(4));
        let _ = pp.insert(&PauliOperator::X(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::X(4));
        let _ = pp.insert(&PauliOperator::Y(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Y(4));
        let _ = pp.insert(&PauliOperator::X(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::X(4));
        let _ = pp.insert(&PauliOperator::X(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Y(4));
        let _ = pp.insert(&PauliOperator::Y(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
    }
    
    #[test]
    fn transform_double_excitation_ops_success_1() {
        // Fermion Operator
        // (0,1,3,2),(0,3,1,2)
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(3),
            CrAnOperator::Annihilation(2),
            0.5,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            0.5,
        );

        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(2),
            0.5,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(3),
            CrAnOperator::Annihilation(0),
            0.5,
        );

        // Jordan-Wigner
        let em = em_jordan_wigner();
        let double_excitation_ops = fo.double_excitation_operators().unwrap();
        let qo = transform_double_excitation_ops(&em, &double_excitation_ops).unwrap();

        assert_eq!(qo.len(), 4);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::X(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::X(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        // Bravyi-Kitaev
        let em = em_bravyi_kitaev();
        let double_excitation_ops = fo.double_excitation_operators().unwrap();
        let qo = transform_double_excitation_ops(&em, &double_excitation_ops).unwrap();

        assert_eq!(qo.len(), 4);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.125));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
    }

    #[test]
    fn transform_double_excitation_ops_success_2() {
        // Fermion Operator
        // (0,1,3,2),(0,3,1,2)
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(3),
            CrAnOperator::Annihilation(2),
            0.5,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            0.5,
        );

        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(2),
            0.2,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(3),
            CrAnOperator::Annihilation(0),
            0.2,
        );

        // Jordan-Wigner
        let em = em_jordan_wigner();
        let double_excitation_ops = fo.double_excitation_operators().unwrap();
        let qo = transform_double_excitation_ops(&em, &double_excitation_ops).unwrap();

        assert_eq!(qo.len(), 8);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::X(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.037500000000000006));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.0875));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.037500000000000006));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::X(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0875));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0875));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::X(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.037500000000000006));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::X(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.0875));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.037500000000000006));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        // Bravyi-Kitaev
        let em = em_bravyi_kitaev();
        let double_excitation_ops = fo.double_excitation_operators().unwrap();
        let qo = transform_double_excitation_ops(&em, &double_excitation_ops).unwrap();
	
        assert_eq!(qo.len(), 8);
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0875));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0875));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.037500000000000006));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.037500000000000006));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0875));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0875));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.037500000000000006));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.037500000000000006));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
    }

    #[test]
    fn transform_double_excitation_ops_success_3() {
        // Fermion Operator
        // (0,1,5,2),(0,3,4,2)
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(5),
            CrAnOperator::Annihilation(2),
            0.5,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Creation(5),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            0.5,
        );

        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(4),
            CrAnOperator::Annihilation(2),
            0.2,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Creation(4),
            CrAnOperator::Annihilation(3),
            CrAnOperator::Annihilation(0),
            0.2,
        );

        // Jordan-Wigner
        let em = em_jordan_wigner_6();
        let double_excitation_ops = fo.double_excitation_operators().unwrap();
        let qo = transform_double_excitation_ops(&em, &double_excitation_ops).unwrap();

        assert_eq!(qo.len(), 16);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = pp.insert(&PauliOperator::Z(4));
        let _ = pp.insert(&PauliOperator::X(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = pp.insert(&PauliOperator::Z(4));
        let _ = pp.insert(&PauliOperator::Y(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = pp.insert(&PauliOperator::Z(4));
        let _ = pp.insert(&PauliOperator::Y(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = pp.insert(&PauliOperator::Z(4));
        let _ = pp.insert(&PauliOperator::X(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::X(3));
        let _ = pp.insert(&PauliOperator::X(4));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.025));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        let _ = pp.insert(&PauliOperator::Y(4));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.025));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::X(3));
        let _ = pp.insert(&PauliOperator::Y(4));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.025));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        let _ = pp.insert(&PauliOperator::X(4));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.025));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = pp.insert(&PauliOperator::Z(4));
        let _ = pp.insert(&PauliOperator::Y(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = pp.insert(&PauliOperator::Z(4));
        let _ = pp.insert(&PauliOperator::X(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = pp.insert(&PauliOperator::Z(4));
        let _ = pp.insert(&PauliOperator::X(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = pp.insert(&PauliOperator::Z(4));
        let _ = pp.insert(&PauliOperator::Y(5));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.0625));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::X(3));
        let _ = pp.insert(&PauliOperator::Y(4));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.025));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        let _ = pp.insert(&PauliOperator::X(4));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.025));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::X(3));
        let _ = pp.insert(&PauliOperator::X(4));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.025));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
	
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        let _ = pp.insert(&PauliOperator::Y(4));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.025));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
    }
}

//Result of Openfermion

// [0]
//fermion hamiltonian
//0.5 [0^ 1^ 4 5] +
//0.5 [5^ 4^ 1 0]
//jw
//(-0.0625+0j) [X0 X1 X4 X5] +
//(0.0625+0j) [X0 X1 Y4 Y5] +
//(-0.0625+0j) [X0 Y1 X4 Y5] +
//(-0.0625+0j) [X0 Y1 Y4 X5] +
//(-0.0625+0j) [Y0 X1 X4 Y5] +
//(-0.0625+0j) [Y0 X1 Y4 X5] +
//(0.0625+0j) [Y0 Y1 X4 X5] +
//(-0.0625+0j) [Y0 Y1 Y4 Y5]
//bk
//(-0.0625+0j) [X0 Z1 X4] +
//(-0.0625+0j) [X0 Z1 X4 Z5] +
//(-0.0625+0j) [X0 X4] +
//(-0.0625+0j) [X0 X4 Z5] +
//(-0.0625+0j) [Y0 Z1 Y4] +
//(-0.0625+0j) [Y0 Z1 Y4 Z5] +
//(-0.0625+0j) [Y0 Y4] +
//(-0.0625+0j) [Y0 Y4 Z5]

// [1]
//fermion hamiltonian
//0.5 [0^ 1^ 3 2] +
//0.5 [0^ 3^ 1 2] +
//0.5 [2^ 1^ 3 0] +
//0.5 [2^ 3^ 1 0]
//jw
//(-0.125+0j) [X0 X1 Y2 Y3] +
//(0.125+0j) [X0 Y1 Y2 X3] +
//(0.125+0j) [Y0 X1 X2 Y3] +
//(-0.125+0j) [Y0 Y1 X2 X3]
//bk
//(0.125+0j) [X0 Z1 X2] +
//(0.125+0j) [X0 Z1 X2 Z3] +
//(0.125+0j) [Y0 Z1 Y2] +
//(0.125+0j) [Y0 Z1 Y2 Z3]

// [2]
//fermion hamiltonian
//0.5 [0^ 1^ 3 2] +
//0.2 [0^ 3^ 1 2] +
//0.2 [2^ 1^ 3 0] +
//0.5 [2^ 3^ 1 0]
//jw
//(0.037500000000000006+0j) [X0 X1 X2 X3] +
//(-0.0875+0j) [X0 X1 Y2 Y3] +
//(0.037500000000000006+0j) [X0 Y1 X2 Y3] +
//(0.0875+0j) [X0 Y1 Y2 X3] +
//(0.0875+0j) [Y0 X1 X2 Y3] +
//(0.037500000000000006+0j) [Y0 X1 Y2 X3] +
//(-0.0875+0j) [Y0 Y1 X2 X3] +
//(0.037500000000000006+0j) [Y0 Y1 Y2 Y3]
//bk
//(0.0875+0j) [X0 Z1 X2] +
//(0.0875+0j) [X0 Z1 X2 Z3] +
//(0.037500000000000006+0j) [X0 X2] +
//(0.037500000000000006+0j) [X0 X2 Z3] +
//(0.0875+0j) [Y0 Z1 Y2] +
//(0.0875+0j) [Y0 Z1 Y2 Z3] +
//(0.037500000000000006+0j) [Y0 Y2] +
//(0.037500000000000006+0j) [Y0 Y2 Z3]

// [3]
//fermion hamiltonian
//0.5 [0^ 1^ 5 2] +
//0.2 [0^ 3^ 4 2] +
//0.2 [2^ 4^ 3 0] +
//0.5 [2^ 5^ 1 0]
//jw
//(0.0625+0j) [X0 X1 X2 Z3 Z4 X5] +
//(-0.0625+0j) [X0 X1 Y2 Z3 Z4 Y5] +
//(0.0625+0j) [X0 Y1 X2 Z3 Z4 Y5] +
//(0.0625+0j) [X0 Y1 Y2 Z3 Z4 X5] +
//(0.025+0j) [X0 Z1 X2 X3 X4] +
//(0.025+0j) [X0 Z1 X2 Y3 Y4] +
//(-0.025+0j) [X0 Z1 Y2 X3 Y4] +
//(0.025+0j) [X0 Z1 Y2 Y3 X4] +
//(0.0625+0j) [Y0 X1 X2 Z3 Z4 Y5] +
//(0.0625+0j) [Y0 X1 Y2 Z3 Z4 X5] +
//(-0.0625+0j) [Y0 Y1 X2 Z3 Z4 X5] +
//(0.0625+0j) [Y0 Y1 Y2 Z3 Z4 Y5] +
//(0.025+0j) [Y0 Z1 X2 X3 Y4] +
//(-0.025+0j) [Y0 Z1 X2 Y3 X4] +
//(0.025+0j) [Y0 Z1 Y2 X3 X4] +
//(0.025+0j) [Y0 Z1 Y2 Y3 Y4]
//bk
//(-0.025+0j) [X0 X1 X2 Y3 Y4 X5] +
//(-0.025+0j) [X0 X1 Y2 Y3 X4 X5] +
//(0.025+0j) [X0 Y1 X2 X3 Y4 X5] +
//(0.025+0j) [X0 Y1 Y2 X3 X4 X5] +
//(-0.0625+0j) [X0 Z1 X2 Y3 Y5] +
//(-0.0625+0j) [X0 Z1 Y2 Y3 Z4 X5] +
//(-0.0625+0j) [X0 X2 Y3 Y5] +
//(-0.0625+0j) [X0 Y2 Y3 Z4 X5] +
//(0.025+0j) [Y0 X1 X2 Y3 X4 X5] +
//(-0.025+0j) [Y0 X1 Y2 Y3 Y4 X5] +
//(-0.025+0j) [Y0 Y1 X2 X3 X4 X5] +
//(0.025+0j) [Y0 Y1 Y2 X3 Y4 X5] +
//(0.0625+0j) [Y0 Z1 X2 Y3 Z4 X5] +
//(-0.0625+0j) [Y0 Z1 Y2 Y3 Y5] +
//(0.0625+0j) [Y0 X2 Y3 Z4 X5] +
//(-0.0625+0j) [Y0 Y2 Y3 Y5]
