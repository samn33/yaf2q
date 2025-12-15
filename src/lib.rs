pub mod binary_matrix;
pub mod f2q_mapper;
pub mod fermion_operator;
pub mod qubit_operator;
pub mod ternary_tree_spec;
pub mod transform_constant;
pub mod transform_coulomb_exchange_ops;
pub mod transform_double_excitation_ops;
pub mod transform_excitation_ops;
pub mod transform_number_excitation_ops;
pub mod transform_number_ops;
pub mod util;

use num_complex::Complex;
use pyo3::prelude::*;
use std::collections::HashMap;

use crate::f2q_mapper::{f2q_encoding_matrix, f2q_operator, F2QKind, F2QMapper};
use crate::fermion_operator::FermionOperator;
use crate::qubit_operator::PauliProduct;
use crate::ternary_tree_spec::TernayTreeSpec;

/// Make fermion to qubit mapper
///
/// # Arguments
///
/// * `kind` - kind of fermion to qubit transformation (jordan-wigner/parity/bravyi-kitaev)
/// * `num_qubits` - number of qubits
/// * `ttspec_str` - ternary tree generator (string)
///
fn make_mapper(
    kind: Option<&str>,
    num_qubits: usize,
    ttspec_str: Option<&str>,
) -> Result<F2QMapper, String> {
    let mapper = match ttspec_str {
        Some(t) => match kind {
            None => F2QMapper::TernaryTreeMapper(TernayTreeSpec::from_str(t).unwrap()),
            _ => {
                return Err("Invalid kind string is specified.".to_string());
            }
        },
        None => match kind {
            Some("jordan-wigner") => F2QMapper::NamedMapper(F2QKind::JordanWigner, num_qubits),
            Some("parity") => F2QMapper::NamedMapper(F2QKind::Parity, num_qubits),
            Some("bravyi-kitaev") => F2QMapper::NamedMapper(F2QKind::BravyiKitaev, num_qubits),
            _ => {
                return Err("Invalid kind string is specified.".to_string());
            }
        },
    };

    Ok(mapper)
}

/// Get the encoding matrix
///
/// # Arguments
///
/// * `kind` - kind of fermion to qubit transformation (jordan-wigner/parity/bravyi-kitaev)
/// * `ttspec_str` - ternary tree generator (string)
///
#[pyfunction]
fn yaf2q_encoding_matrix(
    kind: Option<&str>,
    ttspec_str: Option<&str>,
    num_qubits: usize,
) -> PyResult<Vec<Vec<i32>>> {
    let mapper = make_mapper(kind, num_qubits, ttspec_str).unwrap();
    let em_tmp = f2q_encoding_matrix(&mapper).unwrap();

    let mut em: Vec<Vec<i32>> = vec![];
    for vec in em_tmp.elements.iter() {
        let mut row: Vec<i32> = vec![];
        for b in vec.iter() {
            row.push(*b as i32);
        }
        em.push(row);
    }

    Ok(em)
}

/// Get the inverse of encoding matrix
///
/// # Arguments
///
/// * `kind` - kind of fermion to qubit transformation (jordan-wigner/parity/bravyi-kitaev)
/// * `ttspec_str` - ternary tree generator (string)
///
#[pyfunction]
fn yaf2q_encoding_matrix_inv(
    kind: Option<&str>,
    ttspec_str: Option<&str>,
    num_qubits: usize,
) -> PyResult<Vec<Vec<i32>>> {
    let mapper = make_mapper(kind, num_qubits, ttspec_str).unwrap();
    let em_tmp = f2q_encoding_matrix(&mapper).unwrap().inverse().unwrap();

    let mut em: Vec<Vec<i32>> = vec![];
    for vec in em_tmp.elements.iter() {
        let mut row: Vec<i32> = vec![];
        for b in vec.iter() {
            row.push(*b as i32);
        }
        em.push(row);
    }

    Ok(em)
}

/// Transform the fermion operator (string) to qubit operator
///
/// # Arguments
///
/// * `kind` - kind of fermion to qubit transformation (jordan-wigner/parity/bravyi-kitaev)
/// * `ttspec_str` - ternary tree generator (string)
/// * `fh_str` - fermion operator (string)
///
#[pyfunction]
fn yaf2q_operator(
    kind: Option<&str>,
    ttspec_str: Option<&str>,
    fh_str: &str,
) -> PyResult<HashMap<String, Complex<f64>>> {
    let fh = FermionOperator::from_str(fh_str).unwrap();

    let mapper = make_mapper(kind, fh.num_orbits, ttspec_str).unwrap();

    let qo_tmp = f2q_operator(&mapper, &fh).unwrap();

    let mut qo = HashMap::<String, Complex<f64>>::new();
    for (prod, coef) in qo_tmp.terms.iter() {
        let pp = PauliProduct {
            factor: 0,
            prod: prod.clone(),
        };
        qo.insert(pp.to_string(), *coef);
    }

    Ok(qo)
}

/// A Python module implemented in Rust.
#[pymodule]
fn yaf2q(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(yaf2q_encoding_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(yaf2q_encoding_matrix_inv, m)?)?;
    m.add_function(wrap_pyfunction!(yaf2q_operator, m)?)?;
    Ok(())
}
