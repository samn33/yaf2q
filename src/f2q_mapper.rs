//! Module for Fermion to Qubit Operator

use crate::binary_matrix::BinaryMatrix;
use crate::fermion_operator::FermionOperator;
use crate::qubit_operator::QubitOperator;
use crate::ternary_tree_spec::{TernaryTree, TernayTreeSpec};
use crate::transform_constant::transform_constant;
use crate::transform_coulomb_exchange_ops::transform_coulomb_exchange_ops;
use crate::transform_double_excitation_ops::transform_double_excitation_ops;
use crate::transform_excitation_ops::transform_excitation_ops;
use crate::transform_number_excitation_ops::transform_number_excitation_ops;
use crate::transform_number_ops::transform_number_ops;

/// Enumeration of mapper kinds
///
/// # Fields
///
/// * `JordanWigner` - Jordan-Wigner transform
/// * `Parity` - parity transform
/// * `BravyiKitaev` - Bravyi-Kitaev transform
///
#[derive(Debug, Clone, PartialEq)]
pub enum F2QKind {
    JordanWigner,
    Parity,
    BravyiKitaev,
}

/// Enumeration of fermion to qubit mapper
///
/// # Fields
///
/// * `NamedMapper` - conventional named mapper (Jordan-Wigner, Parity, and Bravyi-Kitaev transform)
/// * `TernaryTreeMapper` - ternary tree mapper (fermion to qubit mapper dedined using ternary tree)
///
#[derive(Debug, Clone, PartialEq)]
pub enum F2QMapper {
    NamedMapper(F2QKind, usize),
    TernaryTreeMapper(TernayTreeSpec),
}

/// Return the encoding matrix of the F2QMapper
///
/// # Arguments
///
/// * `mapper` - Fermion to Qubit Operator mapper
///
pub fn f2q_encoding_matrix(mapper: &F2QMapper) -> Result<BinaryMatrix, String> {
    // Ternary tree object for the kind of fermion to qubit mapping method
    let tree = match mapper {
        F2QMapper::NamedMapper(f2q_kind, num_qubits) => match f2q_kind {
            F2QKind::JordanWigner => Some(TernaryTree::create_jordan_wigner(*num_qubits).unwrap()),
            F2QKind::Parity => Some(TernaryTree::create_parity(*num_qubits).unwrap()),
            F2QKind::BravyiKitaev => Some(TernaryTree::create_bravyi_kitaev(*num_qubits).unwrap()),
        },
        F2QMapper::TernaryTreeMapper(ttspec) => {
            Some(TernaryTree::create_ternary_tree(&ttspec).unwrap())
        }
    };
    let tree = tree.unwrap();

    // Encoding matrix for the ternary tree
    let em = tree.encoding_matrix().unwrap();

    Ok(em)
}

/// Transform the fermion operator to the qubit operator using the F2QMapper
///
/// # Arguments
///
/// * `mapper` - fermion to qubit mapper
/// * `fo` - fermion operator
///
pub fn f2q_operator(mapper: &F2QMapper, fo: &FermionOperator) -> Result<QubitOperator, String> {
    // encodint matrix
    let em = f2q_encoding_matrix(&mapper).unwrap();

    // Number Operators
    let number_ops = fo.number_operators().unwrap();

    // Coulomb/Exchange Operators
    let coulomb_exchange_ops = fo.coulomb_exchange_operators().unwrap();

    // Excitation Operators
    let excitation_ops = fo.excitation_operators().unwrap();

    // Number Excitation Operators
    let number_excitation_ops = fo.number_excitation_operators().unwrap();

    // Double Excitation Operators
    let double_excitation_ops = fo.double_excitation_operators().unwrap();

    // Qubit Operator
    let mut qo = QubitOperator::new(fo.num_orbits);

    // Qubit Operator for constant
    let qo_constant = transform_constant(&em, fo.constant).unwrap();
    qo.add(&qo_constant)
        .expect("Something went wrong adding the Qubit Operator.");

    // Qubit Operator for Number operations
    let qo_number_ops = transform_number_ops(&em, &number_ops).unwrap();
    qo.add(&qo_number_ops)
        .expect("Something went wrong adding the Qubit Operator.");

    // Qubit Operator for Coulomb/Exchange operations
    let qo_coulomb_exchange_ops =
        transform_coulomb_exchange_ops(&em, &coulomb_exchange_ops).unwrap();
    qo.add(&qo_coulomb_exchange_ops)
        .expect("Something went wrong adding the Qubit Operator.");

    // Qubit Operator for Excitation operations
    let qo_excitation_ops = transform_excitation_ops(&em, &excitation_ops).unwrap();
    qo.add(&qo_excitation_ops)
        .expect("Something went wrong adding the Qubit Operator.");

    // Qubit Operator for Number Excitation operations
    let qo_number_excitation_ops =
        transform_number_excitation_ops(&em, &number_excitation_ops).unwrap();
    qo.add(&qo_number_excitation_ops)
        .expect("Something went wrong adding the Qubit Operator.");

    // Qubit Operator for Double Excitation operations
    let qo_double_excitation_ops =
        transform_double_excitation_ops(&em, &double_excitation_ops).unwrap();
    qo.add(&qo_double_excitation_ops)
        .expect("Something went wrong adding the Qubit Operator.");

    Ok(qo)
}

mod common_test_data {

    #[allow(dead_code)]
    pub fn h2_str() -> &'static str {
        "\n
() 0.8141187860307693
((0, 1), (0, 0)) -1.309509868464871
((1, 1), (1, 0)) -1.309509868464871
((2, 1), (2, 0)) -0.4100263808117848
((3, 1), (3, 0)) -0.4100263808117848
((0, 1), (0, 1), (0, 0), (0, 0)) 0.3459522026149021
((0, 1), (0, 1), (2, 0), (2, 0)) 0.08815922580510362
((0, 1), (1, 1), (1, 0), (0, 0)) 0.3459522026149021
((0, 1), (1, 1), (3, 0), (2, 0)) 0.08815922580510362
((0, 1), (2, 1), (0, 0), (2, 0)) 0.08815922580510362
((0, 1), (2, 1), (2, 0), (0, 0)) 0.33984195696523023
((0, 1), (3, 1), (1, 0), (2, 0)) 0.08815922580510362
((0, 1), (3, 1), (3, 0), (0, 0)) 0.33984195696523023
((1, 1), (0, 1), (0, 0), (1, 0)) 0.3459522026149021
((1, 1), (0, 1), (2, 0), (3, 0)) 0.08815922580510362
((1, 1), (1, 1), (1, 0), (1, 0)) 0.3459522026149021
((1, 1), (1, 1), (3, 0), (3, 0)) 0.08815922580510362
((1, 1), (2, 1), (0, 0), (3, 0)) 0.08815922580510362
((1, 1), (2, 1), (2, 0), (1, 0)) 0.33984195696523023
((1, 1), (3, 1), (1, 0), (3, 0)) 0.08815922580510362
((1, 1), (3, 1), (3, 0), (1, 0)) 0.33984195696523023
((2, 1), (0, 1), (0, 0), (2, 0)) 0.3398419569652301
((2, 1), (0, 1), (2, 0), (0, 0)) 0.08815922580510362
((2, 1), (1, 1), (1, 0), (2, 0)) 0.3398419569652301
((2, 1), (1, 1), (3, 0), (0, 0)) 0.08815922580510362
((2, 1), (2, 1), (0, 0), (0, 0)) 0.08815922580510362
((2, 1), (2, 1), (2, 0), (2, 0)) 0.35733555551906787
((2, 1), (3, 1), (1, 0), (0, 0)) 0.08815922580510362
((2, 1), (3, 1), (3, 0), (2, 0)) 0.35733555551906787
((3, 1), (0, 1), (0, 0), (3, 0)) 0.3398419569652301
((3, 1), (0, 1), (2, 0), (1, 0)) 0.08815922580510362
((3, 1), (1, 1), (1, 0), (3, 0)) 0.3398419569652301
((3, 1), (1, 1), (3, 0), (1, 0)) 0.08815922580510362
((3, 1), (2, 1), (0, 0), (1, 0)) 0.08815922580510362
((3, 1), (2, 1), (2, 0), (3, 0)) 0.35733555551906787
((3, 1), (3, 1), (1, 0), (1, 0)) 0.08815922580510362
((3, 1), (3, 1), (3, 0), (3, 0)) 0.35733555551906787
"
    }
}

#[cfg(test)]
mod tests {
    use super::common_test_data::*;
    use super::*;
    use crate::fermion_operator::FermionOperator;
    use crate::qubit_operator::{PauliOperator, PauliProduct};
    use float_cmp::approx_eq;

    #[test]
    fn f2q_success() {
        let fo = FermionOperator::from_str(h2_str()).unwrap();

        // Jordan-Wigner

        let mapper = F2QMapper::NamedMapper(F2QKind::JordanWigner, 4);
        let qo = f2q_operator(&mapper, &fo).unwrap();

        //I: 0.03775110394645538
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.03775110394645538,
            epsilon = 1.0e-12
        ));

        //X[0]X[1]Y[2]Y[3]: -0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            -0.04407961290255181,
            epsilon = 1.0e-12
        ));

        //X[0]Y[1]Y[2]X[3]: 0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::X(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.04407961290255181,
            epsilon = 1.0e-12
        ));

        //Y[0]X[1]X[2]Y[3]: 0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.04407961290255181,
            epsilon = 1.0e-12
        ));

        //Y[0]Y[1]X[2]X[3]:-0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::X(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            -0.04407961290255181,
            epsilon = 1.0e-12
        ));

        //Z[0]: 0.186016488862306
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.186016488862306,
            epsilon = 1.0e-12
        ));

        //Z[0]Z[1]: 0.17297610130745106
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.17297610130745106,
            epsilon = 1.0e-12
        ));

        //Z[0]Z[2]: 0.12584136558006329
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.12584136558006329,
            epsilon = 1.0e-12
        ));

        //Z[0]Z[3]: 0.1699209784826151
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.1699209784826151,
            epsilon = 1.0e-12
        ));

        //Z[1]: 0.18601648886230598
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.18601648886230598,
            epsilon = 1.0e-12
        ));

        //Z[1]Z[2]: 0.1699209784826151
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.1699209784826151,
            epsilon = 1.0e-12
        ));

        //Z[1]Z[3]: 0.12584136558006329
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.12584136558006329,
            epsilon = 1.0e-12
        ));

        //Z[2]: -0.2694169314163199
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            -0.2694169314163199,
            epsilon = 1.0e-12
        ));

        //Z[2]Z[3]:0.17866777775953394
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.17866777775953394,
            epsilon = 1.0e-12
        ));

        //Z[3]:-0.2694169314163199
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            -0.2694169314163199,
            epsilon = 1.0e-12
        ));

        // Bravyi-Kitaev

        let mapper = F2QMapper::NamedMapper(F2QKind::BravyiKitaev, 4);
        let qo = f2q_operator(&mapper, &fo).unwrap();

        //I: 0.03775110394645538
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.03775110394645538,
            epsilon = 1.0e-12
        ));

        //X[0]Z[1]X[2] :0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.04407961290255181,
            epsilon = 1.0e-12
        ));

        //X[0]Z[1]X[2]Z[3]: 0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.04407961290255181,
            epsilon = 1.0e-12
        ));

        //Y[0]Z[1]Y[2]: 0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.04407961290255181,
            epsilon = 1.0e-12
        ));

        //Y[0]Z[1]Y[2]Z[3]: 0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.04407961290255181,
            epsilon = 1.0e-12
        ));

        //Z[0]: 0.186016488862306
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.186016488862306,
            epsilon = 1.0e-12
        ));

        //Z[0]Z[1]: 0.18601648886230598
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.18601648886230598,
            epsilon = 1.0e-12
        ));

        //Z[0]Z[1]Z[2]: 0.1699209784826151
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.1699209784826151,
            epsilon = 1.0e-12
        ));

        //Z[0]Z[1]Z[2]Z[3]: 0.1699209784826151
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.1699209784826151,
            epsilon = 1.0e-12
        ));

        //Z[0]Z[2]: 0.12584136558006329
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.12584136558006329,
            epsilon = 1.0e-12
        ));

        //Z[0]Z[2]Z[3]: 0.12584136558006329
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.12584136558006329,
            epsilon = 1.0e-12
        ));

        //Z[1]: 0.17297610130745106
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.17297610130745106,
            epsilon = 1.0e-12
        ));

        //Z[1]Z[2]Z[3]: -0.2694169314163199
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            -0.2694169314163199,
            epsilon = 1.0e-12
        ));

        //Z[1]Z[3]: 0.17866777775953394
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(3));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            0.17866777775953394,
            epsilon = 1.0e-12
        ));

        //Z[2]: -0.2694169314163199
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        assert!(approx_eq!(
            f64,
            qo.get_coef(&pp).unwrap().re,
            -0.2694169314163199,
            epsilon = 1.0e-12
        ));
    }

}
