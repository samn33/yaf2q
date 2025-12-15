//! Module for Qubit Operator

use float_cmp::approx_eq;
use ndarray::linalg::kron;
use ndarray::{arr2, Array2};
use num_complex::Complex;
use std::cmp::max;
use std::collections::{BTreeSet, HashMap};
use std::fmt;
use std::fs::File;
use std::io::prelude::*;

/// Enumeration of Pauli Operator (X,Y,Z)
///
/// # Fields
///
/// * `I` - Identity operator
/// * `X` - Pauli X operator
/// * `Y` - Pauli Y operator
/// * `Z` - Pauli Z operator
///
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum PauliOperator {
    I,
    X(usize),
    Y(usize),
    Z(usize),
}

/// Struct of pauli product
///
/// # Fields
///
/// * `num_qubits` - number of qubits (maximum of the qubit indices + 1)
/// * `prod` - Pauli product
/// * `factor` - phase factor (integer) of the Pauli Product (0->1.0, 1->1.0i, 2->-1.0, 3->-1.0i)
///
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct PauliProduct {
    pub factor: u8,
    pub prod: BTreeSet<PauliOperator>,
}

impl PauliProduct {
    /// Return a Pauli product
    ///
    pub fn new() -> Self {
        Self {
            factor: 0,
            prod: BTreeSet::new(),
        }
    }

    /// Get the length of the PauliProduct
    ///
    pub fn len(&self) -> usize {
        self.prod.len()
    }

    /// Check if the Pauli operator is in the PauliProduct
    ///
    /// # Arguments
    ///
    /// * `pauli_op` - pauli operator
    ///
    pub fn contains(&self, pauli_op: &PauliOperator) -> bool {
        for p in &self.prod {
            if p == pauli_op {
                return true;
            }
        }
        false
    }

    /// Check if the pauli operator is exist or not
    ///
    /// # Arguments
    ///
    /// * `pauli_op` - pauli operator
    ///
    pub fn contains_index(&self, pauli_op: &PauliOperator) -> bool {
        let idx = match pauli_op {
            PauliOperator::X(i) => i,
            PauliOperator::Y(i) => i,
            PauliOperator::Z(i) => i,
            _ => {
                return false;
            }
        };
        for p in &self.prod {
            match p {
                PauliOperator::X(i) => {
                    if i == idx {
                        return true;
                    } else {
                    }
                }
                PauliOperator::Y(i) => {
                    if i == idx {
                        return true;
                    } else {
                    }
                }
                PauliOperator::Z(i) => {
                    if i == idx {
                        return true;
                    } else {
                    }
                }
                _ => {}
            };
        }
        false
    }

    /// Check if the PauliProduct is duplicate with another one or not
    ///
    /// # Arguments
    ///
    /// * `other` - pauli product
    ///
    pub fn duplicate_index(&self, other: &Self) -> bool {
        for p in &other.prod {
            if self.contains_index(&p) {
                return true;
            }
        }
        false
    }

    /// Insert the PauliOperator to the PauliProduct
    ///
    /// # Arguments
    ///
    /// * `pauli_op` - pauli operator
    ///
    pub fn insert(&mut self, pauli_op: &PauliOperator) -> Result<(), String> {
        if !self.prod.contains(pauli_op) {
            match pauli_op {
                PauliOperator::X(idx) => match self.pauli_operator(*idx) {
                    Some(PauliOperator::Y(idx)) => {
                        self.prod.remove(&PauliOperator::Y(idx));
                        self.prod.insert(PauliOperator::Z(idx));
                        self.factor = (self.factor + 3) % 4;
                    }
                    Some(PauliOperator::Z(idx)) => {
                        self.prod.remove(&PauliOperator::Z(idx));
                        self.prod.insert(PauliOperator::Y(idx));
                        self.factor = (self.factor + 1) % 4;
                    }
                    _ => {
                        self.prod.insert(*pauli_op);
                    }
                },
                PauliOperator::Y(idx) => match self.pauli_operator(*idx) {
                    Some(PauliOperator::X(idx)) => {
                        self.prod.remove(&PauliOperator::X(idx));
                        self.prod.insert(PauliOperator::Z(idx));
                        self.factor = (self.factor + 1) % 4;
                    }
                    Some(PauliOperator::Z(idx)) => {
                        self.prod.remove(&PauliOperator::Z(idx));
                        self.prod.insert(PauliOperator::X(idx));
                        self.factor = (self.factor + 3) % 4;
                    }
                    _ => {
                        self.prod.insert(*pauli_op);
                    }
                },
                PauliOperator::Z(idx) => match self.pauli_operator(*idx) {
                    Some(PauliOperator::X(idx)) => {
                        self.prod.remove(&PauliOperator::X(idx));
                        self.prod.insert(PauliOperator::Y(idx));
                        self.factor = (self.factor + 3) % 4;
                    }
                    Some(PauliOperator::Y(idx)) => {
                        self.prod.remove(&PauliOperator::Y(idx));
                        self.prod.insert(PauliOperator::X(idx));
                        self.factor = (self.factor + 1) % 4;
                    }
                    _ => {
                        self.prod.insert(*pauli_op);
                    }
                },
                _ => {
                    self.prod.insert(*pauli_op);
                }
            }
        }
	else if self.prod.len() > 1 {
            self.prod.remove(pauli_op);
        }
	else if self.prod.len() == 1 {
            self.prod.remove(pauli_op);
            self.prod.insert(PauliOperator::I);
        }
	else {
            return Err("Something went wrong inserting a pauli operator.".to_string());
	}

        // remove identity operator
        if self.prod.len() > 1 {
            self.prod.remove(&PauliOperator::I);
        }

        Ok(())
    }

    /// Return a PauliOperator corresponding to the index
    ///
    /// # Arguments
    ///
    /// * `index` - index of pauli operator
    ///
    fn pauli_operator(&self, index: usize) -> Option<PauliOperator> {
        for p in self.prod.iter() {
            match p {
                PauliOperator::X(idx) if *idx == index => {
                    return Some(PauliOperator::X(*idx));
                }
                PauliOperator::Y(idx) if *idx == index => {
                    return Some(PauliOperator::Y(*idx));
                }
                PauliOperator::Z(idx) if *idx == index => {
                    return Some(PauliOperator::Z(*idx));
                }
                _ => {
                    continue;
                }
            };
        }
        None
    }

    /// Return a matrix representation of the PauliProduct
    ///
    /// # Arguments
    ///
    /// * `index` - index of pauli operator
    ///
    pub fn matrix(&self, num_qubits_whole: usize) -> Result<Array2<Complex<f64>>, String> {
        let mut mat: Array2<Complex<f64>> = arr2(&[[Complex::new(1.0, 0.0)]]);
        for i in 0..num_qubits_whole {
            match self.pauli_operator(i) {
                Some(pauli) => match pauli {
                    PauliOperator::X(_) => {
                        let p: Array2<Complex<f64>> = arr2(&[
                            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
                            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
                        ]);
                        mat = kron(&mat, &p);
                    }
                    PauliOperator::Y(_) => {
                        let p: Array2<Complex<f64>> = arr2(&[
                            [Complex::new(0.0, 0.0), Complex::new(0.0, -1.0)],
                            [Complex::new(0.0, 1.0), Complex::new(0.0, 0.0)],
                        ]);
                        mat = kron(&mat, &p);
                    }
                    PauliOperator::Z(_) => {
                        let p: Array2<Complex<f64>> = arr2(&[
                            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
                            [Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0)],
                        ]);
                        mat = kron(&mat, &p);
                    }
                    _ => {
                        return Err("Unknown Pauli Operator.".to_string());
                    }
                },
                _ => {
                    let p: Array2<Complex<f64>> = arr2(&[
                        [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
                        [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
                    ]);
                    mat = kron(&(mat.clone()), &p);
                }
            }
        }
        Ok(mat)
    }
}

/// Struct of Qubit Operator
///
/// # Fields
///
/// * `num_qubits` - number of qubits
/// * `terms` - terms of Qubit Operator (key:set of pauli operators, value: (num_qubits,coefficient))
///
//#[pyclass]
#[derive(Debug, Clone, PartialEq)]
pub struct QubitOperator {
    pub num_qubits: usize,
    pub terms: HashMap<BTreeSet<PauliOperator>, Complex<f64>>,
}

impl QubitOperator {
    /// Return a Qubit Operator (empty)
    ///
    pub fn new(num_qubits: usize) -> Self {
        Self {
            num_qubits: num_qubits,
            terms: HashMap::new(),
        }
    }

    /// Get a length of the Qubit Operator (number of terms)
    ///
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    /// Get a number of pauli operators of the Qubit Operator
    ///
    pub fn num_paulis(&self) -> usize {
        let mut num_paulis = 0;
        for (prod, _) in self.terms.iter() {
            num_paulis += prod.len();
        }
        num_paulis
    }

    pub fn num_duplicates(&self) -> usize {
        let mut num_duplicates = 0;
        for (i, (prod_i, _)) in self.terms.iter().enumerate() {
            let pp_i = PauliProduct {
                factor: 0,
                prod: prod_i.clone(),
            };
            for (j, (prod_j, _)) in self.terms.iter().enumerate() {
                if i < j {
                    let pp_j = PauliProduct {
                        factor: 0,
                        prod: prod_j.clone(),
                    };
                    num_duplicates += match pp_i.duplicate_index(&pp_j) {
                        true => 1,
                        false => 0,
                    };
                }
            }
        }
        num_duplicates
    }

    /// Get a maximum number of the degrees of the Qubit Operator
    ///
    pub fn max_degree(&self) -> usize {
        let mut max_degree = 0;
        for (prod, _) in self.terms.iter() {
            if prod.len() > max_degree {
                max_degree = prod.len();
            }
        }
        max_degree
    }

    /// Get an average number of the degrees of the Qubit Operator
    ///
    pub fn ave_degree(&self) -> f64 {
        let mut ave_degree = 0.0_f64;
        for (prod, _) in self.terms.iter() {
            ave_degree += prod.len() as f64;
        }
        ave_degree /= self.len() as f64;
        ave_degree
    }

    /// Add a term (Pauli Product with some coefficient) to the Qubit Operator
    ///
    /// # Arguments
    ///
    /// * `pauli_prod` - Pauli product
    /// * `coef` - coefficient for the Pauli product
    ///
    pub fn add_term(
        &mut self,
        pauli_prod: &PauliProduct,
        coef: Complex<f64>,
    ) -> Result<(), String> {
        let factor_0 = Complex::<f64>::new(1.0, 0.0);
        let factor_1 = Complex::<f64>::new(0.0, 1.0);
        let factor_2 = Complex::<f64>::new(-1.0, 0.0);
        let factor_3 = Complex::<f64>::new(0.0, -1.0);

        let prod = &pauli_prod.prod;
        let factor = &pauli_prod.factor;
        let coefficient = coef
            * match factor {
                0 => factor_0,
                1 => factor_1,
                2 => factor_2,
                _ => factor_3,
            };

        if self.terms.contains_key(prod) {
            let c = self.terms.get(prod).unwrap();
            if approx_eq!(f64, (c + coefficient).re, 0.0) && approx_eq!(f64, (c + coefficient).im, 0.0) {
                self.terms.remove(prod);
	    }
	    else {
                self.terms.insert(prod.clone(), c + coefficient);
            }
        }
	else {
            self.terms.insert(prod.clone(), coefficient);
        }
        Ok(())
    }

    /// Add a term (Pauli Product with some coefficient) to the Qubit Operator
    ///
    /// # Arguments
    ///
    /// * `other` - Qubit Operator
    ///
    pub fn add(&mut self, other: &Self) -> Result<(), String> {
        self.num_qubits = max(self.num_qubits, other.num_qubits);
        for (key, value) in other.terms.iter() {
            let coef = *value;
            let prod = key.clone();
            let pp = PauliProduct {
                factor: 0,
                prod: prod,
            };
            let _ = self.add_term(&pp, coef);
        }
        Ok(())
    }

    /// Subtract a term (Pauli Product with some coefficient) to the Qubit Operator
    ///
    /// # Arguments
    ///
    /// * `other` - Qubit Operator
    ///
    pub fn sub(&mut self, other: &Self) -> Result<(), String> {
        self.num_qubits = max(self.num_qubits, other.num_qubits);
        for (key, value) in other.terms.iter() {
            let coef = *value * (-1.0);
            let prod = key.clone();
            let pp = PauliProduct {
                factor: 0,
                prod: prod,
            };
            let _ = self.add_term(&pp, coef);
        }
        Ok(())
    }

    /// Multiply a term (Pauli Product with some coefficient) to the Qubit Operator
    ///
    /// # Arguments
    ///
    /// * `other` - Qubit Operator
    ///
    pub fn mul(&mut self, other: &Self) -> Result<(), String> {
	let num_qubits = max(self.num_qubits, other.num_qubits);
        let mut qo = QubitOperator::new(num_qubits);
        for (key_self, value_self) in self.terms.iter() {
            let coef_self = *value_self;
            let prod_self = key_self.clone();
            let pp_self = PauliProduct {
                factor: 0,
                prod: prod_self,
            };
            for (key_other, value_other) in other.terms.iter() {
                let coef_other = *value_other;
                let prod_other = key_other.clone();
                let mut pp = pp_self.clone();
                for p_other in prod_other.iter() {
                    let _ = pp.insert(&p_other);
                }
                let _ = qo.add_term(&pp, coef_self * coef_other);
            }
        }
        *self = qo.clone();
        Ok(())
    }

    /// Get a coefficient for the term of the Qubit Operator
    ///
    /// # Arguments
    ///
    /// * `pauli_prod` - Pauli product
    /// * `coef` - coefficient for the Pauli product
    ///
    pub fn get_coef(&self, pauli_prod: &PauliProduct) -> Option<Complex<f64>> {
        match self.terms.get(&pauli_prod.prod) {
            Some(value) => Some(*value),
            _ => None,
        }
    }

    /// Save the Qubit Operator to a file
    ///
    /// # Arguments
    ///
    /// * `qo_path` - Qubit Operator file path
    ///
    pub fn save(&self, qo_path: &str) -> Result<(), String> {
        let mut f = File::create(qo_path).expect("File not created.");
        f.write_all(self.to_string().as_bytes())
            .expect("Something went wrong writing the file.");
        Ok(())
    }

    /// Check if the Qubit Operator is hermitian or not
    ///
    pub fn is_hermitian(&self) -> bool {
        for (_, value) in self.terms.iter() {
            let coef = value;
            if !approx_eq!(f64, coef.im, 0.0) {
                return false;
            }
        }
        return true;
    }

    /// Get a Matrix representation of the Qubit Operator
    ///
    pub fn matrix(&self) -> Result<Array2<Complex<f64>>, String> {
        let mat_size = 1 << self.num_qubits;
        let mut mat: Array2<Complex<f64>> = Array2::zeros((mat_size, mat_size));
        for (key, value) in self.terms.iter() {
            let coef = *value;
            let prod = key.clone();
            let pp = PauliProduct {
                factor: 0,
                prod: prod,
            };
            mat = &mat + (coef * &(pp.matrix(self.num_qubits)).unwrap());
        }
        Ok(mat)
    }
}

impl fmt::Display for PauliOperator {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self {
            Self::I => {
                let _ = write!(f, "");
            }
            Self::X(idx) => {
                let _ = write!(f, "X{} ", idx);
            }
            Self::Y(idx) => {
                let _ = write!(f, "Y{} ", idx);
            }
            Self::Z(idx) => {
                let _ = write!(f, "Z{} ", idx);
            }
        }
        write!(f, "")
    }
}

impl fmt::Display for PauliProduct {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for p in &self.prod {
            let _ = write!(f, "{}", p);
        }
        write!(f, "")
    }
}

impl fmt::Display for QubitOperator {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (key, value) in &self.terms {
            let _ = writeln!(
                f,
                "{}: {}",
                PauliProduct {
                    factor: 0,
                    prod: key.clone()
                },
                *value
            );
        }
        let _ = writeln!(f, "---");
        let _ = writeln!(f, "number of terms      : {}", self.len());
        let _ = writeln!(f, "number of paulis     : {}", self.num_paulis());
        let _ = writeln!(f, "number of duplicates : {}", self.num_duplicates());
        let _ = writeln!(f, "maximum degree       : {}", self.max_degree());
        let _ = writeln!(f, "average degree       : {}", self.ave_degree());
        write!(f, "")
    }
}

mod common_test_data {
    use super::*;

    #[allow(dead_code)]
    pub fn simple_qo() -> QubitOperator {
        let mut qo = QubitOperator::new(5);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = qo.add_term(&pp, 1.0.into());

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(3));
        let _ = pp.insert(&PauliOperator::Z(4));
        let _ = qo.add_term(&pp, 1.0.into());

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        let _ = qo.add_term(&pp, 1.0.into());

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(4));
        let _ = qo.add_term(&pp, 1.0.into());

        qo
    }

    #[allow(dead_code)]
    pub fn h2_qo_jordan_wigner() -> QubitOperator {
        let mut qo = QubitOperator::new(4);

        //I: 0.03775110394645538
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        let _ = qo.add_term(&pp, 0.03775110394645538.into());

        //X[0]X[1]Y[2]Y[3]: -0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        let _ = qo.add_term(&pp, (-0.04407961290255181).into());

        //X[0]Y[1]Y[2]X[3]: 0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::X(3));
        let _ = qo.add_term(&pp, 0.04407961290255181.into());

        //Y[0]X[1]X[2]Y[3]: 0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Y(3));
        let _ = qo.add_term(&pp, 0.04407961290255181.into());

        //Y[0]Y[1]X[2]X[3]:-0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::X(3));
        let _ = qo.add_term(&pp, (-0.04407961290255181).into());

        //Z[0]: 0.186016488862306
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = qo.add_term(&pp, 0.186016488862306.into());

        //Z[0]Z[1]: 0.17297610130745106
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = qo.add_term(&pp, 0.17297610130745106.into());

        //Z[0]Z[2]: 0.12584136558006329
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = qo.add_term(&pp, 0.12584136558006329.into());

        //Z[0]Z[3]: 0.1699209784826151
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = qo.add_term(&pp, 0.1699209784826151.into());

        //Z[1]: 0.18601648886230598
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = qo.add_term(&pp, 0.18601648886230598.into());

        //Z[1]Z[2]: 0.1699209784826151
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = qo.add_term(&pp, 0.1699209784826151.into());

        //Z[1]Z[3]: 0.12584136558006329
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = qo.add_term(&pp, 0.12584136558006329.into());

        //Z[2]: -0.2694169314163199
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = qo.add_term(&pp, (-0.2694169314163199).into());

        //Z[2]Z[3]: 0.17866777775953394
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = qo.add_term(&pp, 0.17866777775953394.into());

        //Z[3]: -0.2694169314163199
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = qo.add_term(&pp, (-0.2694169314163199).into());

        qo
    }

    #[allow(dead_code)]
    pub fn h2_qo_bravyi_kitaev() -> QubitOperator {
        let mut qo = QubitOperator::new(4);

        //I: 0.03775110394645538
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::I);
        let _ = qo.add_term(&pp, 0.03775110394645538.into());

        //X[0]Z[1]X[2]: 0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = qo.add_term(&pp, 0.04407961290255181.into());

        //X[0]Z[1]X[2]Z[3]: 0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::X(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = qo.add_term(&pp, 0.04407961290255181.into());

        //Y[0]Z[1]Y[2]: 0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = qo.add_term(&pp, 0.04407961290255181.into());

        //Y[0]Z[1]Y[2]Z[3]: 0.04407961290255181
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = qo.add_term(&pp, 0.04407961290255181.into());

        //Z[0]: 0.186016488862306
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = qo.add_term(&pp, 0.186016488862306.into());

        //Z[0]Z[1]: 0.18601648886230598
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = qo.add_term(&pp, 0.18601648886230598.into());

        //Z[0]Z[1]Z[2]: 0.1699209784826151
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = qo.add_term(&pp, 0.1699209784826151.into());

        //Z[0]Z[1]Z[2]Z[3]: 0.1699209784826151
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = qo.add_term(&pp, 0.1699209784826151.into());

        //Z[0]Z[2]: 0.12584136558006329
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = qo.add_term(&pp, 0.12584136558006329.into());

        //Z[0]Z[2]Z[3]: 0.12584136558006329
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = qo.add_term(&pp, 0.12584136558006329.into());

        //Z[1]: 0.17297610130745106
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = qo.add_term(&pp, 0.17297610130745106.into());

        //Z[1]Z[2]Z[3]: -0.2694169314163199
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = qo.add_term(&pp, (-0.2694169314163199).into());

        //Z[1]Z[3]: 0.17866777775953394
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        let _ = pp.insert(&PauliOperator::Z(3));
        let _ = qo.add_term(&pp, 0.17866777775953394.into());

        //Z[2]: -0.2694169314163199
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(2));
        let _ = qo.add_term(&pp, (-0.2694169314163199).into());

        qo
    }
}

#[cfg(test)]
mod tests {
    use super::common_test_data::*;
    use super::*;

    #[test]
    fn len_success() {
        let qo = simple_qo();
        assert_eq!(qo.len(), 4);
    }

    #[test]
    fn num_paulis_success() {
        let qo = simple_qo();
        assert_eq!(qo.num_paulis(), 10);
    }

    #[test]
    fn max_degree_success() {
        let qo = simple_qo();
        assert_eq!(qo.max_degree(), 3);
    }

    #[test]
    fn ave_degree_success() {
        let qo = simple_qo();
        assert!(approx_eq!(f64, qo.ave_degree(), 2.5, epsilon = 1.0e-12));
    }

    #[test]
    fn num_duplicates_success() {
        let qo = simple_qo();
        assert_eq!(qo.num_duplicates(), 4);
    }

    #[test]
    fn pauli_product_success() {
        let mut pauli_prod = PauliProduct::new();
        let _ = pauli_prod.insert(&PauliOperator::X(1));
        let _ = pauli_prod.insert(&PauliOperator::Y(2));
        let _ = pauli_prod.insert(&PauliOperator::Z(0));

        assert_eq!(pauli_prod.len(), 3);
        assert_eq!(pauli_prod.contains(&PauliOperator::X(1)), true);
        assert_eq!(pauli_prod.contains(&PauliOperator::Y(2)), true);
        assert_eq!(pauli_prod.contains(&PauliOperator::Z(0)), true);
        assert_eq!(pauli_prod.contains(&PauliOperator::X(0)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Y(1)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Z(2)), false);
        assert_eq!(pauli_prod.factor, 0);

        let mut pauli_prod = PauliProduct::new();
        let _ = pauli_prod.insert(&PauliOperator::X(0));
        let _ = pauli_prod.insert(&PauliOperator::Y(1));
        let _ = pauli_prod.insert(&PauliOperator::Z(2));
        let _ = pauli_prod.insert(&PauliOperator::Z(2));

        assert_eq!(pauli_prod.len(), 2);
        assert_eq!(pauli_prod.contains(&PauliOperator::X(0)), true);
        assert_eq!(pauli_prod.contains(&PauliOperator::Y(1)), true);
        assert_eq!(pauli_prod.contains(&PauliOperator::Z(2)), false);
        assert_eq!(pauli_prod.factor, 0);

        let mut pauli_prod2 = PauliProduct::new();
        let _ = pauli_prod2.insert(&PauliOperator::X(0));
        let _ = pauli_prod2.insert(&PauliOperator::Y(1));

        assert_eq!(pauli_prod, pauli_prod2);

        // Y[2]Y[2]
        let mut pauli_prod3 = PauliProduct::new();
        let _ = pauli_prod3.insert(&PauliOperator::Y(2));
        let _ = pauli_prod3.insert(&PauliOperator::Y(2));
        assert_eq!(pauli_prod3.len(), 1);
        assert_eq!(pauli_prod3.contains(&PauliOperator::I), true);
        assert_eq!(pauli_prod3.contains(&PauliOperator::Y(2)), false);
        assert_eq!(pauli_prod.factor, 0);

        // X[0]Y[0]
        let mut pauli_prod = PauliProduct::new();
        let _ = pauli_prod.insert(&PauliOperator::X(0));
        let _ = pauli_prod.insert(&PauliOperator::Y(0));
        assert_eq!(pauli_prod.len(), 1);
        assert_eq!(pauli_prod.contains(&PauliOperator::X(0)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Y(0)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Z(0)), true);
        assert_eq!(pauli_prod.factor, 1);

        // X[0]Z[0]
        let mut pauli_prod = PauliProduct::new();
        let _ = pauli_prod.insert(&PauliOperator::X(0));
        let _ = pauli_prod.insert(&PauliOperator::Z(0));
        assert_eq!(pauli_prod.len(), 1);
        assert_eq!(pauli_prod.contains(&PauliOperator::X(0)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Y(0)), true);
        assert_eq!(pauli_prod.contains(&PauliOperator::Z(0)), false);
        assert_eq!(pauli_prod.factor, 3);

        // Y[0]X[0]
        let mut pauli_prod = PauliProduct::new();
        let _ = pauli_prod.insert(&PauliOperator::Y(0));
        let _ = pauli_prod.insert(&PauliOperator::X(0));
        assert_eq!(pauli_prod.len(), 1);
        assert_eq!(pauli_prod.contains(&PauliOperator::X(0)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Y(0)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Z(0)), true);
        assert_eq!(pauli_prod.factor, 3);

        // Y[0]Z[0]
        let mut pauli_prod = PauliProduct::new();
        let _ = pauli_prod.insert(&PauliOperator::Y(0));
        let _ = pauli_prod.insert(&PauliOperator::Z(0));
        assert_eq!(pauli_prod.len(), 1);
        assert_eq!(pauli_prod.contains(&PauliOperator::X(0)), true);
        assert_eq!(pauli_prod.contains(&PauliOperator::Y(0)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Z(0)), false);
        assert_eq!(pauli_prod.factor, 1);

        // Z[0]X[0]
        let mut pauli_prod = PauliProduct::new();
        let _ = pauli_prod.insert(&PauliOperator::Z(0));
        let _ = pauli_prod.insert(&PauliOperator::X(0));
        assert_eq!(pauli_prod.len(), 1);
        assert_eq!(pauli_prod.contains(&PauliOperator::X(0)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Y(0)), true);
        assert_eq!(pauli_prod.contains(&PauliOperator::Z(0)), false);
        assert_eq!(pauli_prod.factor, 1);

        // Z[0]Y[0]
        let mut pauli_prod = PauliProduct::new();
        let _ = pauli_prod.insert(&PauliOperator::Z(0));
        let _ = pauli_prod.insert(&PauliOperator::Y(0));
        assert_eq!(pauli_prod.len(), 1);
        assert_eq!(pauli_prod.contains(&PauliOperator::X(0)), true);
        assert_eq!(pauli_prod.contains(&PauliOperator::Y(0)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Z(0)), false);
        assert_eq!(pauli_prod.factor, 3);

        // X[0]Y[0]X[0]
        let mut pauli_prod = PauliProduct::new();
        let _ = pauli_prod.insert(&PauliOperator::X(0));
        let _ = pauli_prod.insert(&PauliOperator::Y(0));
        let _ = pauli_prod.insert(&PauliOperator::X(0));
        assert_eq!(pauli_prod.len(), 1);
        assert_eq!(pauli_prod.contains(&PauliOperator::X(0)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Y(0)), true);
        assert_eq!(pauli_prod.contains(&PauliOperator::Z(0)), false);
        assert_eq!(pauli_prod.factor, 2);

        // Z[1]X[0]Y[0]X[0]
        let mut pauli_prod = PauliProduct::new();
        let _ = pauli_prod.insert(&PauliOperator::Z(1));
        let _ = pauli_prod.insert(&PauliOperator::X(0));
        let _ = pauli_prod.insert(&PauliOperator::Y(0));
        let _ = pauli_prod.insert(&PauliOperator::X(0));
        assert_eq!(pauli_prod.len(), 2);
        assert_eq!(pauli_prod.contains(&PauliOperator::X(0)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Y(0)), true);
        assert_eq!(pauli_prod.contains(&PauliOperator::Z(0)), false);
        assert_eq!(pauli_prod.contains(&PauliOperator::Z(1)), true);
        assert_eq!(pauli_prod.factor, 2);

        // I
        let mut pauli_prod = PauliProduct::new();
        let _ = pauli_prod.insert(&PauliOperator::I);
        assert_eq!(pauli_prod.len(), 1);
        assert_eq!(pauli_prod.contains(&PauliOperator::I), true);
        assert_eq!(pauli_prod.contains(&PauliOperator::Y(2)), false);
        assert_eq!(pauli_prod.factor, 0);
    }

    #[test]
    fn add_term_success() {
        let mut qo = QubitOperator::new(5);
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Z(0));

        let mut pp2 = PauliProduct::new();
        let _ = pp2.insert(&PauliOperator::X(3));
        let _ = pp2.insert(&PauliOperator::Z(4));

        let _ = qo.add_term(&pp, 1.0.into());
        let _ = qo.add_term(&pp2, 2.0.into());

        assert_eq!(qo.len(), 2);
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 1.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp2).unwrap().re, 2.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp2).unwrap().im, 0.0));
        assert_eq!(qo.num_qubits, 5);

        let _ = qo.add_term(&pp, 1.0.into());
        assert_eq!(qo.len(), 2);
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 2.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp2).unwrap().re, 2.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp2).unwrap().im, 0.0));
        assert_eq!(qo.num_qubits, 5);

        let _ = qo.add_term(&pp2, (-2.0).into());
        assert_eq!(qo.len(), 1);
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().re, 2.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp).unwrap().im, 0.0));
        assert_eq!(qo.get_coef(&pp2), None);
        assert_eq!(qo.num_qubits, 5);

        let qo2 = qo.clone();
        assert_eq!(qo, qo2);
    }

    #[test]
    fn add_term_considering_factor_success() {
        let mut qo = QubitOperator::new(1);

        let mut pp_0 = PauliProduct::new();
        let _ = pp_0.insert(&PauliOperator::X(0));
        let _ = pp_0.insert(&PauliOperator::Y(0));
        let _ = qo.add_term(&pp_0, 1.0.into());

        assert_eq!(qo.len(), 1);
        assert_eq!(qo.num_qubits, 1);
        assert!(approx_eq!(f64, qo.get_coef(&pp_0).unwrap().re, 0.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp_0).unwrap().im, 1.0));

        let mut pp_1 = PauliProduct::new();
        let _ = pp_1.insert(&PauliOperator::Z(0));
        let _ = qo.add_term(&pp_1, 1.0.into());

        assert_eq!(qo.len(), 1);
        assert_eq!(qo.num_qubits, 1);
        assert!(approx_eq!(f64, qo.get_coef(&pp_1).unwrap().re, 1.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp_1).unwrap().im, 1.0));

        let mut pp_2 = PauliProduct::new();
        let _ = pp_2.insert(&PauliOperator::Y(0));
        let _ = pp_2.insert(&PauliOperator::X(0));
        let _ = qo.add_term(&pp_2, 1.0.into());

        assert_eq!(qo.len(), 1);
        assert_eq!(qo.num_qubits, 1);
        assert!(approx_eq!(f64, qo.get_coef(&pp_2).unwrap().re, 1.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp_2).unwrap().im, 0.0));

        let mut pp_3 = PauliProduct::new();
        let _ = pp_3.insert(&PauliOperator::Z(0));
        let _ = qo.add_term(&pp_3, (-1.0).into());

        assert_eq!(qo.len(), 0);
        assert_eq!(qo.num_qubits, 1);

        let mut pp_4 = PauliProduct::new();
        let _ = pp_4.insert(&PauliOperator::X(0));
        let _ = pp_4.insert(&PauliOperator::X(0));
        let _ = pp_4.insert(&PauliOperator::I);
        let _ = qo.add_term(&pp_4, 4.0.into());

        assert_eq!(qo.len(), 1);
        assert_eq!(qo.num_qubits, 1);
        assert!(approx_eq!(f64, qo.get_coef(&pp_4).unwrap().re, 4.0));
        assert!(approx_eq!(f64, qo.get_coef(&pp_4).unwrap().im, 0.0));
    }

    #[test]
    fn add_term_falilure() {
        let mut qo = QubitOperator::new(5);
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(1));
        let _ = pp.insert(&PauliOperator::Y(2));
        let _ = pp.insert(&PauliOperator::Z(0));

        let mut pp2 = PauliProduct::new();
        let _ = pp2.insert(&PauliOperator::X(3));
        let _ = pp2.insert(&PauliOperator::Z(4));

        let _ = qo.add_term(&pp, 1.0.into());

        assert_eq!(qo.get_coef(&pp2), None);
    }

    #[test]
    fn add_success() {
        let mut qo_a = QubitOperator::new(5);
        let mut pp1 = PauliProduct::new();
        let _ = pp1.insert(&PauliOperator::X(1));
        let _ = pp1.insert(&PauliOperator::Y(2));
        let _ = pp1.insert(&PauliOperator::Z(0));
        let mut pp2 = PauliProduct::new();
        let _ = pp2.insert(&PauliOperator::X(3));
        let _ = pp2.insert(&PauliOperator::Z(4));
        let _ = qo_a.add_term(&pp1, 1.0.into());
        let _ = qo_a.add_term(&pp2, 2.0.into());

        let mut qo_b = QubitOperator::new(5);
        let mut pp3 = PauliProduct::new();
        let _ = pp3.insert(&PauliOperator::X(0));
        let _ = pp3.insert(&PauliOperator::Y(2));
        let _ = pp3.insert(&PauliOperator::Z(1));
        let mut pp4 = PauliProduct::new();
        let _ = pp4.insert(&PauliOperator::X(3));
        let _ = pp4.insert(&PauliOperator::Z(4));
        let _ = qo_b.add_term(&pp3, 3.0.into());
        let _ = qo_b.add_term(&pp4, 4.0.into());

        let mut qo_ab = QubitOperator::new(5);
        let _ = qo_ab.add(&qo_a);
        let _ = qo_ab.add(&qo_b);

        assert_eq!(qo_ab.len(), 3);
        assert!(approx_eq!(f64, qo_ab.get_coef(&pp1).unwrap().re, 1.0));
        assert!(approx_eq!(f64, qo_ab.get_coef(&pp1).unwrap().im, 0.0));
        assert!(approx_eq!(f64, qo_ab.get_coef(&pp3).unwrap().re, 3.0));
        assert!(approx_eq!(f64, qo_ab.get_coef(&pp3).unwrap().im, 0.0));
        assert!(approx_eq!(f64, qo_ab.get_coef(&pp2).unwrap().re, 6.0));
        assert!(approx_eq!(f64, qo_ab.get_coef(&pp2).unwrap().im, 0.0));
        assert_eq!(qo_ab.num_qubits, 5);
    }

    #[test]
    fn mul_success() {
        // (X[0] + 2 Z[0])(3 X[1] + 4 Z[1])
        let mut qo_a = QubitOperator::new(2);
        let mut pp1 = PauliProduct::new();
        let _ = pp1.insert(&PauliOperator::X(0));
        let mut pp2 = PauliProduct::new();
        let _ = pp2.insert(&PauliOperator::Z(0));
        let _ = qo_a.add_term(&pp1, 1.0.into());
        let _ = qo_a.add_term(&pp2, 2.0.into());

        let mut qo_b = QubitOperator::new(2);
        let mut pp3 = PauliProduct::new();
        let _ = pp3.insert(&PauliOperator::X(1));
        let mut pp4 = PauliProduct::new();
        let _ = pp4.insert(&PauliOperator::Z(1));
        let _ = qo_b.add_term(&pp3, 3.0.into());
        let _ = qo_b.add_term(&pp4, 4.0.into());

        let _ = qo_a.mul(&qo_b);

        assert_eq!(qo_a.len(), 4);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::X(1));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().re, 3.0));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().re, 4.0));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::X(1));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().re, 6.0));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().re, 8.0));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().im, 0.0));

        // (I + Z[0])(X[0] - Z[1]) = X[0] - Z[1] + i Y[0] - Z[0]Z[1]
        let mut qo_a = QubitOperator::new(2);
        let mut pp1 = PauliProduct::new();
        let _ = pp1.insert(&PauliOperator::I);
        let mut pp2 = PauliProduct::new();
        let _ = pp2.insert(&PauliOperator::Z(0));
        let _ = qo_a.add_term(&pp1, 1.0.into());
        let _ = qo_a.add_term(&pp2, 1.0.into());

        let mut qo_b = QubitOperator::new(2);
        let mut pp3 = PauliProduct::new();
        let _ = pp3.insert(&PauliOperator::X(0));
        let mut pp4 = PauliProduct::new();
        let _ = pp4.insert(&PauliOperator::Z(1));
        let _ = qo_b.add_term(&pp3, 1.0.into());
        let _ = qo_b.add_term(&pp4, (-1.0).into());

        let _ = qo_a.mul(&qo_b);

        assert_eq!(qo_a.len(), 4);

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().re, 1.0));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().re, -1.0));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().im, 0.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Y(0));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().re, 0.0));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().im, 1.0));

        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = pp.insert(&PauliOperator::Z(1));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().re, -1.0));
        assert!(approx_eq!(f64, qo_a.get_coef(&pp).unwrap().im, 0.0));
    }

    #[test]
    fn matrix_success() {
        // Z[0]
        let mut qo_a = QubitOperator::new(1);
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::Z(0));
        let _ = qo_a.add_term(&pp, 1.0.into());
        let m_actual = qo_a.matrix().unwrap();
        let m_expect: Array2<Complex<f64>> = arr2(&[
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0)],
        ]);
        assert_eq!(m_actual, m_expect);

        // 2.0 X[0] Y[1]
        let mut qo_b = QubitOperator::new(2);
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = qo_b.add_term(&pp, 2.0.into());
        let m_actual = qo_b.matrix().unwrap();
        let m_expect: Array2<Complex<f64>> = arr2(&[
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -2.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 2.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -2.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 2.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
        ]);
        assert_eq!(m_actual, m_expect);

        // Z[0] + 2.0 X[0] Y[1]
        let _ = qo_a.add(&qo_b);
        let m_actual = qo_a.matrix().unwrap();
        let m_expect: Array2<Complex<f64>> = arr2(&[
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -2.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 2.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -2.0),
                Complex::new(-1.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 2.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(-1.0, 0.0),
            ],
        ]);
        assert_eq!(m_actual, m_expect);
    }

    #[test]
    fn is_hermitian_success() {
        let mut qo = QubitOperator::new(6);
        let mut pp = PauliProduct::new();
        let _ = pp.insert(&PauliOperator::X(0));
        let _ = pp.insert(&PauliOperator::Y(1));
        let _ = pp.insert(&PauliOperator::Z(2));

        let mut pp2 = PauliProduct::new();
        let _ = pp2.insert(&PauliOperator::X(3));
        let _ = pp2.insert(&PauliOperator::Z(4));

        let _ = qo.add_term(&pp, 1.0.into());
        let _ = qo.add_term(&pp2, 2.0.into());

        assert_eq!(qo.is_hermitian(), true);

        let mut pp3 = PauliProduct::new();
        let _ = pp3.insert(&PauliOperator::X(5));
        let _ = pp3.insert(&PauliOperator::Y(5));
        let _ = qo.add_term(&pp3, 3.0.into());

        assert_eq!(qo.is_hermitian(), false);
    }
}
