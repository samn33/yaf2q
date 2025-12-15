//! Module for transformation from the Excitation operators to a Qubit Operator

use crate::binary_matrix::BinaryMatrix;
use crate::fermion_operator::ExcitationOperators;
use crate::qubit_operator::{PauliOperator, PauliProduct, QubitOperator};

/// Transform excitation operators to QubitOperator
///
/// # Arguments
///
/// * `em` - encoding matrix
/// * `ops` - excitation operators
///
#[allow(unused_variables, unused_mut)]
pub fn transform_excitation_ops(
    em: &BinaryMatrix,
    ops: &ExcitationOperators,
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
    let mut one_i = BinaryMatrix::zero(num_qubits, 1).unwrap();
    let mut one_j = BinaryMatrix::zero(num_qubits, 1).unwrap();

    for ((i, j), coef) in ops.iter() {
        for n in 0..num_qubits {
            e_i.elements[n][0] = 0;
            e_j.elements[n][0] = 0;
            one_i.elements[n][0] = 0;
            one_j.elements[n][0] = 0;
        }
        e_i.elements[*i][0] = 1;
        e_j.elements[*j][0] = 1;

        // v_ij_t
        let e_ij = e_i.add(&e_j).unwrap();
        let v_ij = em.mul(&e_ij).unwrap();
        let v_ij_t = v_ij.transpose().unwrap();

        // e_i_t
        let mut e_i_t = e_i.transpose().unwrap();
        e_i_t = e_i_t.mul(&em_inv).unwrap();

        // e_j_t
        let mut e_j_t = e_j.transpose().unwrap();
        e_j_t = e_j_t.mul(&em_inv).unwrap();

        // one_ij_t
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
        }
        let one_ij = one_i.add(&one_j).unwrap();
        let mut one_ij_t = one_ij.transpose().unwrap();
        one_ij_t = one_ij_t.mul(&em_inv).unwrap();

        // Qubit Operator: qo_l
        let mut qo_l = QubitOperator::new(num_qubits);
        let mut pp = PauliProduct::new();
        for n in 0..num_qubits {
            match v_ij_t.elements[0][n] {
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
        let _ = qo_l.add_term(&pp, (*coef).into());

        // Qubit Operator: qo_r
        let mut qo_r = QubitOperator::new(num_qubits);
        let mut pp = PauliProduct::new();
        for n in 0..num_qubits {
            match v_ij_t.elements[0][n] {
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
        let _ = qo_r.add_term(&pp, (*coef).into());

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
        let _ = qo_m_1.add_term(&pp, (-0.5).into());
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        let _ = qo_m_1.add_term(&pp, (0.5).into());

        // Qubit Operator: qo_m_2
        let mut qo_m_2 = QubitOperator::new(num_qubits);
        let mut pp = PauliProduct::new();
        for n in 0..num_qubits {
            match one_ij_t.elements[0][n] {
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
        let _ = qo_m_2.add_term(&pp, (1.0).into());

        // Qubit Operator: qo_m = qo_m_0 * qo_m_1 * qo_m_2
        let mut qo_m = qo_m_2.clone();
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
    fn transform_excitation_ops_success_0() {
        // Fermion Operator
        let mut fo = FermionOperator::new();
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(1),
            0.5,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(0),
            0.5,
        );

        let _ = fo.insert_second_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(2),
            0.5,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(1),
            0.5,
        );

        let _ = fo.insert_second_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(3),
            0.5,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(2),
            0.5,
        );

        let _ = fo.insert_second_order(
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(0),
            0.5,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(3),
            0.5,
        );

        // Jordan-Wigner
        let em = em_jordan_wigner();
        let excitation_ops = fo.excitation_operators().unwrap();
        let qo = transform_excitation_ops(&em, &excitation_ops).unwrap();

        assert_eq!(qo.len(), 8);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::X(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::X(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::X(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        // Bravyi-Kitaev
        let em = em_bravyi_kitaev();
        let excitation_ops = fo.excitation_operators().unwrap();
        let qo = transform_excitation_ops(&em, &excitation_ops).unwrap();

        assert_eq!(qo.len(), 8);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::X(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, -0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
    }

    #[test]
    fn transform_excitation_ops_success_1() {
        // Fermion Operator
	// (0,1), (1,3)
	// (0.5, 0.2)
        let mut fo = FermionOperator::new();
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(1),
            0.5,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(0),
            0.5,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(3),
            0.2,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(1),
            0.2,
        );

        // Jordan-Wigner
        let em = em_jordan_wigner();
        let excitation_ops = fo.excitation_operators().unwrap();
        let qo = transform_excitation_ops(&em, &excitation_ops).unwrap();

        assert_eq!(qo.len(), 4);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::X(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
    }

    #[test]
    fn transform_excitation_ops_success_2() {
        // Fermion Operator
	// (0,2), (1,3), (2,4)
	// (0.5, 0.2, 0.1)
        let mut fo = FermionOperator::new();
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(2),
            0.5,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(0),
            0.5,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(3),
            0.2,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(3),
            CrAnOperator::Annihilation(1),
            0.2,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(4),
            0.1,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(4),
            CrAnOperator::Annihilation(2),
            0.1,
        );

        // Jordan-Wigner
        let em = em_jordan_wigner_6();
        let excitation_ops = fo.excitation_operators().unwrap();
        let qo = transform_excitation_ops(&em, &excitation_ops).unwrap();

        assert_eq!(qo.len(), 6);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.25));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::X(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.1));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = pp.insert(&PauliOperator::X(4));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.05));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = pp.insert(&PauliOperator::Y(4));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 0.05));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
    }
}

//Result of Openfermion

// [0]
//fermion hamiltonian
//0.5 [0^ 1] +
//0.5 [0^ 3] +
//0.5 [1^ 0] +
//0.5 [1^ 2] +
//0.5 [2^ 1] +
//0.5 [2^ 3] +
//0.5 [3^ 0] +
//0.5 [3^ 2]
//jw
//(0.25+0j) [X0 X1] +
//(0.25+0j) [X0 Z1 Z2 X3] +
//(0.25+0j) [Y0 Y1] +
//(0.25+0j) [Y0 Z1 Z2 Y3] +
//(0.25+0j) [X1 X2] +
//(0.25+0j) [Y1 Y2] +
//(0.25+0j) [X2 X3] +
//(0.25+0j) [Y2 Y3]
//bk
//(0.25+0j) [X0] +
//(-0.25+0j) [X0 X1 Z3] +
//(-0.25+0j) [X0 Z1] +
//(-0.25+0j) [Y0 Y1 Z2] +
//(0.25+0j) [Z0 Y1 Y2] +
//(0.25+0j) [X1 X2] +
//(-0.25+0j) [Z1 X2 Z3] +
//(0.25+0j) [X2]

// [1]
//fermion hamiltonian
//0.5 [0^ 1] +
//0.5 [1^ 0] +
//0.2 [1^ 3] +
//0.2 [3^ 1]
//jw
//(0.25+0j) [X0 X1] +
//(0.25+0j) [Y0 Y1] +
//(0.1+0j) [X1 Z2 X3] +
//(0.1+0j) [Y1 Z2 Y3]
//bk
//(0.25+0j) [X0] +
//(-0.25+0j) [X0 Z1] +
//(-0.1+0j) [Z0 X1 Z3] +
//(0.1+0j) [X1 Z2]

// [2]
//fermion hamiltonian
//0.5 [0^ 2] +
//0.2 [1^ 3] +
//0.5 [2^ 0] +
//0.1 [2^ 4] +
//0.2 [3^ 1] +
//0.1 [4^ 2]
//jw
//(0.25+0j) [X0 Z1 X2] +
//(0.25+0j) [Y0 Z1 Y2] +
//(0.1+0j) [X1 Z2 X3] +
//(0.1+0j) [Y1 Z2 Y3] +
//(0.05+0j) [X2 Z3 X4] +
//(0.05+0j) [Y2 Z3 Y4]
//bk
//(0.25+0j) [X0 Y1 Y2] +
//(-0.25+0j) [Y0 Y1 X2] +
//(-0.1+0j) [Z0 X1 Z3] +
//(0.1+0j) [X1 Z2] +
//(0.05+0j) [Z1 X2 Y3 Y4] +
//(-0.05+0j) [Z1 Y2 Y3 X4]
