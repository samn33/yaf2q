//! Module for Fermion Operator

use float_cmp::approx_eq;
use regex::Regex;
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::prelude::*;

/// Enumeration of fermion operator (creation and annihilation operator)
///
/// # Fields
///
/// * `Creation` - creation operator
/// * `Annililation` - annihilation operator
///
#[derive(Debug, Clone, Copy, Eq, Hash, PartialEq)]
pub enum CrAnOperator {
    Creation(usize),
    Annihilation(usize),
}

pub type NumberOperators = HashMap<usize, f64>;
pub type CoulombExchangeOperators = HashMap<(usize, usize), f64>;
pub type ExcitationOperators = HashMap<(usize, usize), f64>;
pub type NumberExcitationOperators = HashMap<(usize, usize, usize), f64>;
pub type DoubleExcitationOperators = HashMap<(usize, usize, usize, usize), f64>;

/// Struct of Fermion Operator
///
/// # Fields
///
/// * `constant` - constant term
/// * `second_orders` - second order terms
/// * `fourth_orders` - fourth order terms
///
#[derive(Debug, Clone, PartialEq)]
pub struct FermionOperator {
    pub num_orbits: usize,
    pub constant: f64,
    pub second_orders: HashMap<(CrAnOperator, CrAnOperator), f64>,
    pub fourth_orders: HashMap<(CrAnOperator, CrAnOperator, CrAnOperator, CrAnOperator), f64>,
}

impl FermionOperator {
    /// Return a Fermion Operator (empty)
    ///
    pub fn new() -> Self {
        Self {
            num_orbits: 0,
            constant: 0.0,
            second_orders: HashMap::new(),
            fourth_orders: HashMap::new(),
        }
    }

    /// Insert the constant value to the FermionOperator
    ///
    /// # Arguments
    ///
    /// * `constant` - constant value
    ///
    pub fn insert_constant_value(&mut self, constant: f64) -> Result<(), String> {
        self.constant += constant;
        Ok(())
    }

    /// Insert the second order term to the FermionOperator
    ///
    /// # Arguments
    ///
    /// * `fop_i` - fermion operator for i-th orbit
    /// * `fop_j` - fermion operator for j-th orbit
    /// * `coef` - coefficient
    ///
    pub fn insert_second_order(
        &mut self,
        fop_i: CrAnOperator,
        fop_j: CrAnOperator,
        coef: f64,
    ) -> Result<(), String> {
        let fop_i = match fop_i {
            CrAnOperator::Creation(i) => CrAnOperator::Creation(i),
            CrAnOperator::Annihilation(_) => {
                return Err("The Fermion Operator is invalid.".to_string());
            }
        };
        let fop_j = match fop_j {
            CrAnOperator::Creation(_) => {
                return Err("The Fermion Operator is invalid.".to_string());
            }
            CrAnOperator::Annihilation(j) => CrAnOperator::Annihilation(j),
        };

        *self.second_orders.entry((fop_i, fop_j)).or_insert(0.0) += coef;
        self.num_orbits = match fop_i {
            CrAnOperator::Creation(i) if i + 1 > self.num_orbits => i + 1,
            _ => self.num_orbits,
        };
        self.num_orbits = match fop_j {
            CrAnOperator::Annihilation(j) if j + 1 > self.num_orbits => j + 1,
            _ => self.num_orbits,
        };
        Ok(())
    }

    /// Insert the fourth order term to the FermionOperator
    ///
    /// # Arguments
    ///
    /// * `fop_i` - fermion operator for i-th orbit
    /// * `fop_j` - fermion operator for j-th orbit
    /// * `fop_k` - fermion operator for k-th orbit
    /// * `fop_l` - fermion operator for l-th orbit
    /// * `coef` - coefficient
    ///
    pub fn insert_fourth_order(
        &mut self,
        fop_i: CrAnOperator,
        fop_j: CrAnOperator,
        fop_k: CrAnOperator,
        fop_l: CrAnOperator,
        coef: f64,
    ) -> Result<(), String> {
        let fop_i = match fop_i {
            CrAnOperator::Creation(i) => CrAnOperator::Creation(i),
            CrAnOperator::Annihilation(_) => {
                return Err("The Fermion Operator is invalid.".to_string());
            }
        };
        let fop_j = match fop_j {
            CrAnOperator::Creation(j) => CrAnOperator::Creation(j),
            CrAnOperator::Annihilation(_) => {
                return Err("The Fermion Operator is invalid.".to_string());
            }
        };
        let fop_k = match fop_k {
            CrAnOperator::Creation(_) => {
                return Err("The Fermion Operator is invalid.".to_string());
            }
            CrAnOperator::Annihilation(k) => CrAnOperator::Annihilation(k),
        };
        let fop_l = match fop_l {
            CrAnOperator::Creation(_) => {
                return Err("The Fermion Operator is invalid.".to_string());
            }
            CrAnOperator::Annihilation(l) => CrAnOperator::Annihilation(l),
        };

        match (fop_i, fop_j, fop_k, fop_l) {
            (
                CrAnOperator::Creation(i),
                CrAnOperator::Creation(j),
                CrAnOperator::Annihilation(k),
                CrAnOperator::Annihilation(l),
            ) => {
                match (i, j, k, l) {
                    (i, j, k, l) if i == j || k == l => {
                        return Ok(());
                    }

		    (i, j, k, l) if j == k && i != l => {
                        *self
                            .fourth_orders
                            .entry((fop_i, fop_j, fop_k, fop_l))
                            .or_insert(0.0) += coef;
                    }
                    (i, j, k, l) if i == l && j != k => {
                        *self
                            .fourth_orders
                            .entry((fop_j, fop_i, fop_l, fop_k))
                            .or_insert(0.0) += coef;
                    }
		    (i, j, k, l) if i == k && j != l => {
                        *self
                            .fourth_orders
                            .entry((fop_i, fop_j, fop_l, fop_k))
                            .or_insert(0.0) -= coef;
                    }
		    (i, j, k, l) if j == l && i != k => {
                        *self
                            .fourth_orders
                            .entry((fop_i, fop_j, fop_l, fop_k))
                            .or_insert(0.0) -= coef;
                    }
		    (i, j, k, l) if j == k && i == l && i < j => {
                        *self
                            .fourth_orders
                            .entry((fop_i, fop_j, fop_k, fop_l))
                            .or_insert(0.0) += coef;
                    }
                    (i, j, k, l) if j == k && i == l && i > j => {
                        *self
                            .fourth_orders
                            .entry((fop_j, fop_i, fop_l, fop_k))
                            .or_insert(0.0) += coef;
                    }

		    (i, j, k, l) if i == k && j == l && i < j => {
                        *self
                            .fourth_orders
                            .entry((fop_i, fop_j, fop_l, fop_k))
                            .or_insert(0.0) -= coef;
                    }
                    (i, j, k, l) if i == k && j == l && i > j => {
                        *self
                            .fourth_orders
                            .entry((fop_j, fop_i, fop_k, fop_l))
                            .or_insert(0.0) -= coef;
                    }

		    _ => {
                        self.fourth_orders
                            .insert((fop_i, fop_j, fop_k, fop_l), coef);
                    }
                };
                self.num_orbits = match fop_i {
                    CrAnOperator::Creation(i) if i + 1 > self.num_orbits => i + 1,
                    _ => self.num_orbits,
                };
                self.num_orbits = match fop_j {
                    CrAnOperator::Creation(j) if j + 1 > self.num_orbits => j + 1,
                    _ => self.num_orbits,
                };
                self.num_orbits = match fop_k {
                    CrAnOperator::Annihilation(k) if k + 1 > self.num_orbits => k + 1,
                    _ => self.num_orbits,
                };
                self.num_orbits = match fop_l {
                    CrAnOperator::Annihilation(l) if l + 1 > self.num_orbits => l + 1,
                    _ => self.num_orbits,
                };
                Ok(())
            }
            _ => Ok(()),
        }
    }

    /// Get the FermionOperator from the text
    ///
    /// # Arguments
    ///
    /// * `contents` - Text describing the Fermion Operator
    ///
    pub fn from_str(contents: &str) -> Result<Self, String> {
        let mut fo = Self::new();
        let re_constant = Regex::new(r"\(\)([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)").unwrap();
        let re_second_orders =
          Regex::new(r"\(\((\d+),(\d+)\),\((\d+),(\d+)\)\)([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)").unwrap();
        let re_fourth_orders = Regex::new(r"\(\((\d+),(\d+)\),\((\d+),(\d+)\),\((\d+),(\d+)\),\((\d+),(\d+)\)\)([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)").unwrap();

        for l in contents.lines() {
            let mut line = l.to_string();
            line.retain(|c| !c.is_whitespace());
            match re_constant.captures(&line) {
                Some(caps) => {
                    fo.constant = (caps[1]).parse().unwrap();
                }
                None => (),
            }
            match re_second_orders.captures(&line) {
                Some(caps) => {
                    let (i, i_dg, j, j_dg, coef) = (
                        caps[1].parse::<usize>().unwrap(),
                        caps[2].parse::<usize>().unwrap(),
                        caps[3].parse::<usize>().unwrap(),
                        caps[4].parse::<usize>().unwrap(),
                        caps[5].parse::<f64>().unwrap(),
                    );
                    let fop_i = match i_dg {
                        0 => CrAnOperator::Annihilation(i),
                        _ => CrAnOperator::Creation(i),
                    };
                    let fop_j = match j_dg {
                        0 => CrAnOperator::Annihilation(j),
                        _ => CrAnOperator::Creation(j),
                    };
                    let _ = fo.insert_second_order(fop_i, fop_j, coef);
                }
                None => (),
            }
            match re_fourth_orders.captures(&line) {
                Some(caps) => {
                    let (i, i_dg, j, j_dg, k, k_dg, l, l_dg, coef) = (
                        caps[1].parse::<usize>().unwrap(),
                        caps[2].parse::<usize>().unwrap(),
                        caps[3].parse::<usize>().unwrap(),
                        caps[4].parse::<usize>().unwrap(),
                        caps[5].parse::<usize>().unwrap(),
                        caps[6].parse::<usize>().unwrap(),
                        caps[7].parse::<usize>().unwrap(),
                        caps[8].parse::<usize>().unwrap(),
                        caps[9].parse::<f64>().unwrap(),
                    );
                    let fop_i = match i_dg {
                        0 => CrAnOperator::Annihilation(i),
                        _ => CrAnOperator::Creation(i),
                    };
                    let fop_j = match j_dg {
                        0 => CrAnOperator::Annihilation(j),
                        _ => CrAnOperator::Creation(j),
                    };
                    let fop_k = match k_dg {
                        0 => CrAnOperator::Annihilation(k),
                        _ => CrAnOperator::Creation(k),
                    };
                    let fop_l = match l_dg {
                        0 => CrAnOperator::Annihilation(l),
                        _ => CrAnOperator::Creation(l),
                    };
                    let _ = fo.insert_fourth_order(fop_i, fop_j, fop_k, fop_l, coef);
                }
                None => (),
            }
        }

        match fo.is_hermitian() {
            true => Ok(fo),
            false => Err("The Fermion Operator is not a hermitian.".to_string()),
        }
    }

    /// Load the FermionOperator from the file
    ///
    /// # Arguments
    ///
    /// * `fo_path` - file path describing the FermionOperator
    ///
    pub fn load(fo_path: &str) -> Result<Self, String> {
        let mut f = File::open(fo_path).expect("File not found.");
        let mut contents = String::new();
        f.read_to_string(&mut contents)
            .expect("Something went wrong reading the file.");
        Self::from_str(&contents)
    }

    
    /// Save the FermionOperator to the file
    ///
    /// # Arguments
    ///
    /// * `fo_path` - Fermion Operator file path
    ///
    pub fn save(&self, fo_path: &str) -> Result<(), String> {
        let mut f = File::create(fo_path).expect("File not created.");
        f.write_all(self.to_string().as_bytes())
            .expect("Something went wrong writeing the file.");
        Ok(())
    }

    /// Check if the FermionOperator hermitian or not
    ///
    pub fn is_hermitian(&self) -> bool {
        // second order terms
        for (op, coef) in &self.second_orders {
            let idx_0 = match op.0 {
                CrAnOperator::Annihilation(_) => None,
                CrAnOperator::Creation(idx) => Some(idx),
            };
            let idx_1 = match op.1 {
                CrAnOperator::Annihilation(idx) => Some(idx),
                CrAnOperator::Creation(_) => None,
            };
            match (idx_0, idx_1) {
                (Some(i), Some(j)) if i == j => {
                    continue;
                }
                (Some(i), Some(j)) if i != j => {
                    let op_inv = (CrAnOperator::Creation(j), CrAnOperator::Annihilation(i));
                    if self.second_orders.get(&op_inv) != None
                        && approx_eq!(
                            f64,
                            *self.second_orders.get(&op_inv).unwrap(),
                            *coef,
                            epsilon = 1.0e-12
                        )
                    {
                        continue;
                    } else {
                        return false;
                    }
                }
                _ => {
                    return false;
                }
            }
        }

        // fourth order terms
        for (op, coef) in &self.fourth_orders {
            let idx_0 = match op.0 {
                CrAnOperator::Annihilation(_) => None,
                CrAnOperator::Creation(idx) => Some(idx),
            };
            let idx_1 = match op.1 {
                CrAnOperator::Annihilation(_) => None,
                CrAnOperator::Creation(idx) => Some(idx),
            };
            let idx_2 = match op.2 {
                CrAnOperator::Annihilation(idx) => Some(idx),
                CrAnOperator::Creation(_) => None,
            };
            let idx_3 = match op.3 {
                CrAnOperator::Annihilation(idx) => Some(idx),
                CrAnOperator::Creation(_) => None,
            };
            match (idx_0, idx_1, idx_2, idx_3) {
                (Some(i), Some(j), Some(k), Some(l)) if i != j && k != l => {
                    let op_inv = (
                        CrAnOperator::Creation(l),
                        CrAnOperator::Creation(k),
                        CrAnOperator::Annihilation(j),
                        CrAnOperator::Annihilation(i),
                    );
                    if self.fourth_orders.get(&op_inv) != None
                        && approx_eq!(
                            f64,
                            *self.fourth_orders.get(&op_inv).unwrap(),
                            *coef,
                            epsilon = 1.0e-12
                        )
                    {
                        continue;
                    } else {
                        return false;
                    }
                }
                (Some(i), Some(j), Some(k), Some(l)) if i == j || k == l => {
                    continue;
                }
                _ => {
                    return false;
                }
            }
        }

        return true;
    }

    /// Return the Number Operators
    ///
    pub fn number_operators(&self) -> Result<NumberOperators, String> {
        let mut number_ops = NumberOperators::new();
        for (op, coef) in &self.second_orders {
            let idx_0 = match op.0 {
                CrAnOperator::Annihilation(_) => None,
                CrAnOperator::Creation(idx) => Some(idx),
            };
            let idx_1 = match op.1 {
                CrAnOperator::Annihilation(idx) => Some(idx),
                CrAnOperator::Creation(_) => None,
            };
            match (idx_0, idx_1) {
                (Some(i), Some(j)) if i == j => {
                    number_ops.insert(i, *coef);
                }
                _ => (),
            }
        }
	Ok(number_ops)
    }

    /// Return the Coulomb/Exchange Operators
    ///
    pub fn coulomb_exchange_operators(&self) -> Result<CoulombExchangeOperators, String> {
        let mut coulomb_exchange_ops = CoulombExchangeOperators::new();
        for (op, coef) in &self.fourth_orders {
            let idx_0 = match op.0 {
                CrAnOperator::Annihilation(_) => None,
                CrAnOperator::Creation(idx) => Some(idx),
            };
            let idx_1 = match op.1 {
                CrAnOperator::Annihilation(_) => None,
                CrAnOperator::Creation(idx) => Some(idx),
            };
            let idx_2 = match op.2 {
                CrAnOperator::Annihilation(idx) => Some(idx),
                CrAnOperator::Creation(_) => None,
            };
            let idx_3 = match op.3 {
                CrAnOperator::Annihilation(idx) => Some(idx),
                CrAnOperator::Creation(_) => None,
            };
            match (idx_0, idx_1, idx_2, idx_3) {
                (Some(i), Some(j), Some(k), Some(l)) if i == l && j == k && i != j => {
                    coulomb_exchange_ops.insert((i, j), *coef);
                }
                (Some(i), Some(j), Some(k), Some(l)) if i == k && j == l && i != j => {
                    coulomb_exchange_ops.insert((i, j), *coef * (-1.0));
                }
                _ => (),
            }
        }
	Ok(coulomb_exchange_ops)
    }

    /// Return the Excitation Operators
    ///
    pub fn excitation_operators(&self) -> Result<ExcitationOperators, String> {
        let mut excitation_ops = ExcitationOperators::new();
        for (op, coef) in &self.second_orders {
            let idx_0 = match op.0 {
                CrAnOperator::Annihilation(_) => None,
                CrAnOperator::Creation(idx) => Some(idx),
            };
            let idx_1 = match op.1 {
                CrAnOperator::Annihilation(idx) => Some(idx),
                CrAnOperator::Creation(_) => None,
            };
            match (idx_0, idx_1) {
                (Some(i), Some(j)) if i < j => {
                    excitation_ops.insert((i, j), *coef);
                }
                _ => (),
            }
        }
	Ok(excitation_ops)
    }

    /// Return the Number Excitation Operators
    ///
    pub fn number_excitation_operators(&self) -> Result<NumberExcitationOperators, String> {
        let mut number_excitation_ops = NumberExcitationOperators::new();
        for (op, coef) in &self.fourth_orders {
            let idx_0 = match op.0 {
                CrAnOperator::Annihilation(_) => None,
                CrAnOperator::Creation(idx) => Some(idx),
            };
            let idx_1 = match op.1 {
                CrAnOperator::Annihilation(_) => None,
                CrAnOperator::Creation(idx) => Some(idx),
            };
            let idx_2 = match op.2 {
                CrAnOperator::Annihilation(idx) => Some(idx),
                CrAnOperator::Creation(_) => None,
            };
            let idx_3 = match op.3 {
                CrAnOperator::Annihilation(idx) => Some(idx),
                CrAnOperator::Creation(_) => None,
            };
            match (idx_0, idx_1, idx_2, idx_3) {
                (Some(i), Some(j), Some(k), Some(l)) if j == k && i != j && k != l && i < l => {
	    	    *number_excitation_ops.entry((i, j, l)).or_insert(0.0) += *coef;
                }
                (Some(i), Some(j), Some(k), Some(l)) if i == l && i != j && k != l && j < k => {
	    	    *number_excitation_ops.entry((j, i, k)).or_insert(0.0) += *coef;
                }
                (Some(i), Some(j), Some(k), Some(l)) if i == k && i != j && k != l && j < l => {
	    	    *number_excitation_ops.entry((j, i, l)).or_insert(0.0) -= *coef;
                }
                (Some(i), Some(j), Some(k), Some(l)) if j == l && i != j && k != l && i < k => {
	    	    *number_excitation_ops.entry((i, j, k)).or_insert(0.0) -= *coef;
                }
                _ => (),
            }
        }
	Ok(number_excitation_ops)
    }

    /// Return the Double Excitation Operators
    ///
    pub fn double_excitation_operators(&self) -> Result<DoubleExcitationOperators, String> {
        let mut double_excitation_ops = DoubleExcitationOperators::new();
        for (op, coef) in &self.fourth_orders {
            let idx_0 = match op.0 {
                CrAnOperator::Annihilation(_) => None,
                CrAnOperator::Creation(idx) => Some(idx),
            };
            let idx_1 = match op.1 {
                CrAnOperator::Annihilation(_) => None,
                CrAnOperator::Creation(idx) => Some(idx),
            };
            let idx_2 = match op.2 {
                CrAnOperator::Annihilation(idx) => Some(idx),
                CrAnOperator::Creation(_) => None,
            };
            let idx_3 = match op.3 {
                CrAnOperator::Annihilation(idx) => Some(idx),
                CrAnOperator::Creation(_) => None,
            };
            match (idx_0, idx_1, idx_2, idx_3) {
                (Some(i), Some(j), Some(k), Some(l))
                    if i != j && i != k && i != l && j != k && j != l && k != l =>
                {
		    if i < j && k < l {
		    	*double_excitation_ops.entry((i,j,k,l)).or_insert(0.0) += *coef;
		    }
		    else if j < i && l < k {
		    	*double_excitation_ops.entry((j,i,l,k)).or_insert(0.0) += *coef;
		    }
		    else if i < j && l < k {
		    	*double_excitation_ops.entry((i,j,l,k)).or_insert(0.0) -= *coef;
		    }
		    else { // j < i && k < l
		    	*double_excitation_ops.entry((j,i,k,l)).or_insert(0.0) -= *coef;
		    }
                }
                _ => (),
            }
        }

        let mut double_excitation_ops_out = DoubleExcitationOperators::new();

        for (idx, coef) in double_excitation_ops {
            match (idx.0, idx.1, idx.2, idx.3) {
                (i, j, k, l)
		    if double_excitation_ops_out.get(&(i, j, k, l)) == None
		    && double_excitation_ops_out.get(&(k, l, i, j)) == None =>
		{
		    double_excitation_ops_out.insert((i,j,k,l), coef);
                }
                _ => {}
            }
        }
	
	Ok(double_excitation_ops_out)
    }
}

impl fmt::Display for CrAnOperator {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self {
            Self::Annihilation(idx) => {
                let _ = write!(f, "a[{}]", idx);
            }
            Self::Creation(idx) => {
                let _ = write!(f, "a^[{}]", idx);
            }
        }
        write!(f, "")
    }
}

impl fmt::Display for FermionOperator {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let _ = writeln!(f, "num_orbits: {}", self.num_orbits);
        let _ = writeln!(f, "constant: {}", self.constant);
        for (op, coef) in &self.second_orders {
            let _ = writeln!(f, "{}{}: {}", op.0, op.1, coef);
        }
        for (op, coef) in &self.fourth_orders {
            let _ = writeln!(f, "{}{}{}{}: {}", op.0, op.1, op.2, op.3, coef);
        }
        write!(f, "")
    }
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

    #[allow(dead_code)]
    pub fn exp_notation_str() -> &'static str {
        "\n
() 0.8e+10
((0, 1), (0, 0)) -1.3e-1
((0, 1), (1, 1), (1, 0), (0, 0)) 2.6e-10
"
    }
}

#[cfg(test)]
mod tests {
    use super::common_test_data::*;
    use super::*;

    #[test]
    fn new_success() {
        let fo = FermionOperator::new();
        assert_eq!(fo.num_orbits, 0);
        assert_eq!(fo.constant, 0.0);
        assert_eq!(fo.second_orders.len(), 0);
        assert_eq!(fo.fourth_orders.len(), 0);
    }

    #[test]
    fn from_str_exp_notation_success() {
        let fo = FermionOperator::from_str(exp_notation_str()).unwrap();
        assert!(approx_eq!(f64, fo.constant, 0.8e+10));
        assert!(approx_eq!(
            f64,
            *fo.second_orders
                .get(&(CrAnOperator::Creation(0), CrAnOperator::Annihilation(0)))
                .unwrap(),
            -1.3e-1
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(0),
                    CrAnOperator::Creation(1),
                    CrAnOperator::Annihilation(1),
                    CrAnOperator::Annihilation(0)
                ))
                .unwrap(),
            2.6e-10
        ));
    }
    
    #[test]
    fn from_str_h2_success() {
        let fo = FermionOperator::from_str(h2_str()).unwrap();
        let fo2 = FermionOperator::from_str(h2_str()).unwrap();
        assert_eq!(fo, fo2);
        assert_eq!(fo.num_orbits, 4);
        assert!(approx_eq!(f64, fo.constant, 0.8141187860307693));
        assert!(approx_eq!(
            f64,
            *fo.second_orders
                .get(&(CrAnOperator::Creation(0), CrAnOperator::Annihilation(0)))
                .unwrap(),
            -1.309509868464871
        ));
        assert!(approx_eq!(
            f64,
            *fo.second_orders
                .get(&(CrAnOperator::Creation(1), CrAnOperator::Annihilation(1)))
                .unwrap(),
            -1.309509868464871
        ));
        assert!(approx_eq!(
            f64,
            *fo.second_orders
                .get(&(CrAnOperator::Creation(2), CrAnOperator::Annihilation(2)))
                .unwrap(),
            -0.4100263808117848
        ));
        assert!(approx_eq!(
            f64,
            *fo.second_orders
                .get(&(CrAnOperator::Creation(3), CrAnOperator::Annihilation(3)))
                .unwrap(),
            -0.4100263808117848
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(0),
                    CrAnOperator::Creation(1),
                    CrAnOperator::Annihilation(1),
                    CrAnOperator::Annihilation(0)
                ))
                .unwrap(),
            0.3459522026149021 + 0.3459522026149021
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(0),
                    CrAnOperator::Creation(2),
                    CrAnOperator::Annihilation(2),
                    CrAnOperator::Annihilation(0)
                ))
                .unwrap(),
            0.33984195696523023 + 0.3398419569652301 - 0.08815922580510362 - 0.08815922580510362
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(0),
                    CrAnOperator::Creation(3),
                    CrAnOperator::Annihilation(3),
                    CrAnOperator::Annihilation(0)
                ))
                .unwrap(),
            0.33984195696523023 + 0.3398419569652301
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(1),
                    CrAnOperator::Creation(2),
                    CrAnOperator::Annihilation(2),
                    CrAnOperator::Annihilation(1)
                ))
                .unwrap(),
            0.33984195696523023 + 0.3398419569652301
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(1),
                    CrAnOperator::Creation(3),
                    CrAnOperator::Annihilation(3),
                    CrAnOperator::Annihilation(1)
                ))
                .unwrap(),
            0.33984195696523023 + 0.3398419569652301 - 0.08815922580510362 - 0.08815922580510362
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(2),
                    CrAnOperator::Creation(3),
                    CrAnOperator::Annihilation(3),
                    CrAnOperator::Annihilation(2)
                ))
                .unwrap(),
            0.35733555551906787 + 0.35733555551906787
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(1),
                    CrAnOperator::Creation(0),
                    CrAnOperator::Annihilation(2),
                    CrAnOperator::Annihilation(3)
                ))
                .unwrap(),
            0.08815922580510362
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(3),
                    CrAnOperator::Creation(2),
                    CrAnOperator::Annihilation(0),
                    CrAnOperator::Annihilation(1)
                ))
                .unwrap(),
            0.08815922580510362
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(0),
                    CrAnOperator::Creation(3),
                    CrAnOperator::Annihilation(1),
                    CrAnOperator::Annihilation(2)
                ))
                .unwrap(),
            0.08815922580510362
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(2),
                    CrAnOperator::Creation(1),
                    CrAnOperator::Annihilation(3),
                    CrAnOperator::Annihilation(0)
                ))
                .unwrap(),
            0.08815922580510362
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(0),
                    CrAnOperator::Creation(1),
                    CrAnOperator::Annihilation(3),
                    CrAnOperator::Annihilation(2)
                ))
                .unwrap(),
            0.08815922580510362
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(2),
                    CrAnOperator::Creation(3),
                    CrAnOperator::Annihilation(1),
                    CrAnOperator::Annihilation(0)
                ))
                .unwrap(),
            0.08815922580510362
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(1),
                    CrAnOperator::Creation(2),
                    CrAnOperator::Annihilation(0),
                    CrAnOperator::Annihilation(3)
                ))
                .unwrap(),
            0.08815922580510362
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(3),
                    CrAnOperator::Creation(0),
                    CrAnOperator::Annihilation(2),
                    CrAnOperator::Annihilation(1)
                ))
                .unwrap(),
            0.08815922580510362
        ));
    }

    #[test]
    fn from_str_failure() {
        let s = "\n
() 1.0
((0, 1), (0, 0)) 1.0
((0, 1), (1, 0)) 1.0
";
        assert!(FermionOperator::from_str(s).is_err());

        let s = "\n
() 1.0
((0, 1), (0, 0)) 1.0
((0, 1), (1, 1), (2, 0), (3, 0)) 1.0
";
        assert!(FermionOperator::from_str(s).is_err());
    }

    #[test]
    fn insert_constant_value_success() {
        let mut fo = FermionOperator::new();
        let _ = fo.insert_constant_value(1.0);
        assert!(approx_eq!(f64, fo.constant, 1.0));
        assert_eq!(fo.num_orbits, 0);
        assert_eq!(fo.second_orders.len(), 0);
        assert_eq!(fo.fourth_orders.len(), 0);
        assert_eq!(fo.is_hermitian(), true);

        let mut fo = FermionOperator::new();
        let _ = fo.insert_constant_value(1.0);
        let _ = fo.insert_constant_value(2.0);
        assert!(approx_eq!(f64, fo.constant, 3.0));
        assert_eq!(fo.num_orbits, 0);
        assert_eq!(fo.second_orders.len(), 0);
        assert_eq!(fo.fourth_orders.len(), 0);
        assert_eq!(fo.is_hermitian(), true);
    }

    #[test]
    fn insert_second_orders_success() {
        // Number Operator (a^[0] a[0])
        let mut fo = FermionOperator::new();
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        assert!(approx_eq!(
            f64,
            *fo.second_orders
                .get(&(CrAnOperator::Creation(0), CrAnOperator::Annihilation(0)))
                .unwrap(),
            1.0
        ));
        assert_eq!(fo.num_orbits, 1);
        assert_eq!(fo.constant, 0.0);
        assert_eq!(fo.second_orders.len(), 1);
        assert_eq!(fo.fourth_orders.len(), 0);
        assert_eq!(fo.is_hermitian(), true);

        // non-Hermitian (a^[0] a[1])
        let mut fo = FermionOperator::new();
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(1),
            1.0,
        );
        assert!(approx_eq!(
            f64,
            *fo.second_orders
                .get(&(CrAnOperator::Creation(0), CrAnOperator::Annihilation(1)))
                .unwrap(),
            1.0
        ));
        assert_eq!(fo.num_orbits, 2);
        assert_eq!(fo.constant, 0.0);
        assert_eq!(fo.second_orders.len(), 1);
        assert_eq!(fo.fourth_orders.len(), 0);
        assert_eq!(fo.is_hermitian(), false);

        // Excitation Operator (a^[0] a[1] + a^[1] a[0])
        let mut fo = FermionOperator::new();
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(1),
            1.0,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        assert!(approx_eq!(
            f64,
            *fo.second_orders
                .get(&(CrAnOperator::Creation(0), CrAnOperator::Annihilation(1)))
                .unwrap(),
            1.0
        ));
        assert!(approx_eq!(
            f64,
            *fo.second_orders
                .get(&(CrAnOperator::Creation(1), CrAnOperator::Annihilation(0)))
                .unwrap(),
            1.0
        ));
        assert_eq!(fo.num_orbits, 2);
        assert!(approx_eq!(f64, fo.constant, 0.0));
        assert_eq!(fo.second_orders.len(), 2);
        assert_eq!(fo.fourth_orders.len(), 0);
        assert_eq!(fo.is_hermitian(), true);

        // non-Hermitian (a^[0] a[1] + 2.0 a^[1] a[0])
        let mut fo = FermionOperator::new();
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(1),
            1.0,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(0),
            2.0,
        );
        assert!(approx_eq!(
            f64,
            *fo.second_orders
                .get(&(CrAnOperator::Creation(0), CrAnOperator::Annihilation(1)))
                .unwrap(),
            1.0
        ));
        assert!(approx_eq!(
            f64,
            *fo.second_orders
                .get(&(CrAnOperator::Creation(1), CrAnOperator::Annihilation(0)))
                .unwrap(),
            2.0
        ));
        assert_eq!(fo.num_orbits, 2);
        assert!(approx_eq!(f64, fo.constant, 0.0));
        assert_eq!(fo.second_orders.len(), 2);
        assert_eq!(fo.fourth_orders.len(), 0);
        assert_eq!(fo.is_hermitian(), false);
    }

    #[test]
    fn insert_second_orders_failure() {
        // a^[0] a^[1]
        let mut fo = FermionOperator::new();
        assert!(fo
            .insert_second_order(CrAnOperator::Creation(0), CrAnOperator::Creation(1), 1.0)
            .is_err());

        // a[0] a^[1]
        let mut fo = FermionOperator::new();
        assert!(fo
            .insert_second_order(
                CrAnOperator::Annihilation(0),
                CrAnOperator::Creation(1),
                1.0
            )
            .is_err());

        // a[0] a[1]
        let mut fo = FermionOperator::new();
        assert!(fo
            .insert_second_order(
                CrAnOperator::Annihilation(0),
                CrAnOperator::Annihilation(1),
                1.0
            )
            .is_err());
    }

    #[test]
    fn insert_fourth_orders_success() {
        // Coulomb/Exchange Operator (a^[0] a^[1] a[1] a[0])
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(0),
                    CrAnOperator::Creation(1),
                    CrAnOperator::Annihilation(1),
                    CrAnOperator::Annihilation(0)
                ))
                .unwrap(),
            1.0
        ));
        assert_eq!(fo.num_orbits, 2);
        assert!(approx_eq!(f64, fo.constant, 0.0));
        assert_eq!(fo.second_orders.len(), 0);
        assert_eq!(fo.fourth_orders.len(), 1);
        assert_eq!(fo.is_hermitian(), true);

        // Number Excitation Operator (a^[0] a^[1] a[1] a[2] + a^[2] a^[1] a[1] a[0])
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(2),
            1.0,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(0),
                    CrAnOperator::Creation(1),
                    CrAnOperator::Annihilation(1),
                    CrAnOperator::Annihilation(2)
                ))
                .unwrap(),
            1.0
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(2),
                    CrAnOperator::Creation(1),
                    CrAnOperator::Annihilation(1),
                    CrAnOperator::Annihilation(0)
                ))
                .unwrap(),
            1.0
        ));
        assert_eq!(fo.num_orbits, 3);
        assert!(approx_eq!(f64, fo.constant, 0.0));
        assert_eq!(fo.second_orders.len(), 0);
        assert_eq!(fo.fourth_orders.len(), 2);
        assert_eq!(fo.is_hermitian(), true);

        // Double Excitation Operator (a^[0] a^[1] a[2] a[3] + a^[3] a^[2] a[1] a[0])
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(2),
            CrAnOperator::Annihilation(3),
            1.0,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(3),
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(0),
                    CrAnOperator::Creation(1),
                    CrAnOperator::Annihilation(2),
                    CrAnOperator::Annihilation(3)
                ))
                .unwrap(),
            1.0
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(3),
                    CrAnOperator::Creation(2),
                    CrAnOperator::Annihilation(1),
                    CrAnOperator::Annihilation(0)
                ))
                .unwrap(),
            1.0
        ));
        assert_eq!(fo.num_orbits, 4);
        assert!(approx_eq!(f64, fo.constant, 0.0));
        assert_eq!(fo.second_orders.len(), 0);
        assert_eq!(fo.fourth_orders.len(), 2);
        assert_eq!(fo.is_hermitian(), true);

        // non-Hermitian (a^[0] a^[1] a[2] a[3] + 3.0 a^[3] a^[2] a[1] a[0])
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(2),
            CrAnOperator::Annihilation(3),
            1.0,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(3),
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            3.0,
        );
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(0),
                    CrAnOperator::Creation(1),
                    CrAnOperator::Annihilation(2),
                    CrAnOperator::Annihilation(3)
                ))
                .unwrap(),
            1.0
        ));
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(3),
                    CrAnOperator::Creation(2),
                    CrAnOperator::Annihilation(1),
                    CrAnOperator::Annihilation(0)
                ))
                .unwrap(),
            3.0
        ));
        assert_eq!(fo.num_orbits, 4);
        assert!(approx_eq!(f64, fo.constant, 0.0));
        assert_eq!(fo.second_orders.len(), 0);
        assert_eq!(fo.fourth_orders.len(), 2);
        assert_eq!(fo.is_hermitian(), false);

        // non-Hermitian (a^[0] a^[1] a[2] a[3])
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(2),
            CrAnOperator::Annihilation(3),
            1.0,
        );
        assert!(approx_eq!(
            f64,
            *fo.fourth_orders
                .get(&(
                    CrAnOperator::Creation(0),
                    CrAnOperator::Creation(1),
                    CrAnOperator::Annihilation(2),
                    CrAnOperator::Annihilation(3)
                ))
                .unwrap(),
            1.0
        ));
        assert_eq!(fo.num_orbits, 4);
        assert!(approx_eq!(f64, fo.constant, 0.0));
        assert_eq!(fo.second_orders.len(), 0);
        assert_eq!(fo.fourth_orders.len(), 1);
        assert_eq!(fo.is_hermitian(), false);
    }

    #[test]
    fn insert_fourth_orders_failure() {
        // a^[0] a^[1] a^[2] a^[3]
        let mut fo = FermionOperator::new();
        assert!(fo
            .insert_fourth_order(
                CrAnOperator::Creation(0),
                CrAnOperator::Creation(1),
                CrAnOperator::Creation(2),
                CrAnOperator::Creation(3),
                1.0
            )
            .is_err());

        // a^[0] a^[1] a[2] a^[3]
        let mut fo = FermionOperator::new();
        assert!(fo
            .insert_fourth_order(
                CrAnOperator::Creation(0),
                CrAnOperator::Creation(1),
                CrAnOperator::Annihilation(2),
                CrAnOperator::Creation(3),
                1.0
            )
            .is_err());

        // a^[0] a^[1] a^[2] a[3]
        let mut fo = FermionOperator::new();
        assert!(fo
            .insert_fourth_order(
                CrAnOperator::Creation(0),
                CrAnOperator::Creation(1),
                CrAnOperator::Creation(2),
                CrAnOperator::Annihilation(3),
                1.0
            )
            .is_err());

        // a[0] a^[1] a[2] a[3]
        let mut fo = FermionOperator::new();
        assert!(fo
            .insert_fourth_order(
                CrAnOperator::Annihilation(0),
                CrAnOperator::Creation(1),
                CrAnOperator::Annihilation(2),
                CrAnOperator::Annihilation(3),
                1.0
            )
            .is_err());

        // a^[0] a[1] a[2] a[3]
        let mut fo = FermionOperator::new();
        assert!(fo
            .insert_fourth_order(
                CrAnOperator::Creation(0),
                CrAnOperator::Annihilation(1),
                CrAnOperator::Annihilation(2),
                CrAnOperator::Annihilation(3),
                1.0
            )
            .is_err());

        // a[0] a[1] a[2] a[3]
        let mut fo = FermionOperator::new();
        assert!(fo
            .insert_fourth_order(
                CrAnOperator::Annihilation(0),
                CrAnOperator::Annihilation(1),
                CrAnOperator::Annihilation(2),
                CrAnOperator::Annihilation(3),
                1.0
            )
            .is_err());
    }

    #[test]
    fn number_operators_success() {
        // H2
        let fo = FermionOperator::from_str(h2_str()).unwrap();
        let number_ops = fo.number_operators().unwrap();
        assert_eq!(*number_ops.get(&0).unwrap(), -1.309509868464871);
        assert_eq!(*number_ops.get(&1).unwrap(), -1.309509868464871);
        assert_eq!(*number_ops.get(&2).unwrap(), -0.4100263808117848);
        assert_eq!(*number_ops.get(&3).unwrap(), -0.4100263808117848);

        // a^[0] a[0]
        let mut fo = FermionOperator::new();
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        let number_ops = fo.number_operators().unwrap();
        assert!(approx_eq!(f64, *number_ops.get(&(0)).unwrap(), 1.0));
        assert_eq!(number_ops.len(), 1);
    }

    #[test]
    fn coulomb_exchange_operators_success() {
        // H2
        let fo = FermionOperator::from_str(h2_str()).unwrap();
        let coulomb_exchange_ops = fo.coulomb_exchange_operators().unwrap();
        assert!(approx_eq!(
            f64,
            *coulomb_exchange_ops.get(&(0, 1)).unwrap(),
            0.3459522026149021 + 0.3459522026149021
        ));
        assert!(approx_eq!(
            f64,
            *coulomb_exchange_ops.get(&(0, 2)).unwrap(),
            0.33984195696523023 + 0.33984195696523023 - 0.08815922580510362 - 0.08815922580510362
        ));
        assert!(approx_eq!(
            f64,
            *coulomb_exchange_ops.get(&(0, 3)).unwrap(),
            0.33984195696523023 + 0.3398419569652301
        ));
        assert!(approx_eq!(
            f64,
            *coulomb_exchange_ops.get(&(1, 2)).unwrap(),
            0.33984195696523023 + 0.3398419569652301
        ));
        assert!(approx_eq!(
            f64,
            *coulomb_exchange_ops.get(&(1, 3)).unwrap(),
            0.33984195696523023 + 0.3398419569652301 - 0.08815922580510362 - 0.08815922580510362
        ));
        assert!(approx_eq!(
            f64,
            *coulomb_exchange_ops.get(&(2, 3)).unwrap(),
            0.35733555551906787 + 0.35733555551906787
        ));
        assert_eq!(coulomb_exchange_ops.get(&(0, 0)), None);
        assert_eq!(coulomb_exchange_ops.get(&(1, 1)), None);
        assert_eq!(coulomb_exchange_ops.get(&(2, 2)), None);
        assert_eq!(coulomb_exchange_ops.get(&(3, 3)), None);

        // a^[0] a^[1] a[1] a[0]
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        let coulomb_exchange_ops = fo.coulomb_exchange_operators().unwrap();
        assert!(approx_eq!(
            f64,
            *coulomb_exchange_ops.get(&(0, 1)).unwrap(),
            1.0
        ));
        assert_eq!(coulomb_exchange_ops.len(), 1);
    }

    #[test]
    fn excitation_operators_success() {
        // H2
        let fo = FermionOperator::from_str(h2_str()).unwrap();
        let excitation_ops = fo.excitation_operators().unwrap();
        assert_eq!(excitation_ops.get(&(0, 1)), None);
        assert_eq!(excitation_ops.get(&(1, 2)), None);
        assert_eq!(excitation_ops.get(&(2, 3)), None);
        assert_eq!(excitation_ops.get(&(0, 2)), None);

        // a^[0] a[1] + a^[1] a[0]
        let mut fo = FermionOperator::new();
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Annihilation(1),
            1.0,
        );
        let _ = fo.insert_second_order(
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        let excitation_ops = fo.excitation_operators().unwrap();
        assert!(approx_eq!(f64, *excitation_ops.get(&(0, 1)).unwrap(), 1.0));
        assert_eq!(excitation_ops.len(), 1);
    }

    #[test]
    fn number_excitation_operators_success() {
        // H2
        let fo = FermionOperator::from_str(h2_str()).unwrap();
        let number_excitation_ops = fo.number_excitation_operators().unwrap();
        assert_eq!(number_excitation_ops.get(&(0, 1, 2)), None);

        // a^[0] a^[1] a[1] a[2] + a^[2] a^[1] a[1] a[0]
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(2),
            1.0,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(2),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        let number_excitation_ops = fo.number_excitation_operators().unwrap();
        assert!(approx_eq!(
            f64,
            *number_excitation_ops.get(&(0, 1, 2)).unwrap(),
            1.0
        ));
        assert_eq!(number_excitation_ops.len(), 1);
    }

    #[test]
    fn double_excitation_operators_success() {
        // a^[0] a^[1] a[2] a[3] + a^[3] a^[2] a[1] a[0]
        let mut fo = FermionOperator::new();
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(0),
            CrAnOperator::Creation(1),
            CrAnOperator::Annihilation(2),
            CrAnOperator::Annihilation(3),
            1.0,
        );
        let _ = fo.insert_fourth_order(
            CrAnOperator::Creation(3),
            CrAnOperator::Creation(2),
            CrAnOperator::Annihilation(1),
            CrAnOperator::Annihilation(0),
            1.0,
        );
        let double_excitation_ops = fo.double_excitation_operators().unwrap();
        if double_excitation_ops.get(&(0, 1, 2, 3)) != None {
            assert!(approx_eq!(
                f64,
                *double_excitation_ops.get(&(0, 1, 2, 3)).unwrap(),
                1.0
            ));
        } else if double_excitation_ops.get(&(3, 2, 1, 0)) != None {
            assert!(approx_eq!(
                f64,
                *double_excitation_ops.get(&(3, 2, 1, 0)).unwrap(),
                1.0
            ));
        }
        assert_eq!(double_excitation_ops.len(), 1);
    }
}
