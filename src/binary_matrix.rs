//! Module for Binary Matrix

use pyo3::prelude::*;
use std::fmt;
use std::fs::File;
use std::io::prelude::*;

/// Struct of Binary Matrix
///
/// # Fields
///
/// * `elements` - elements of the matrix
///
#[pyclass]
#[derive(Debug, Clone, PartialEq)]
pub struct BinaryMatrix {
    pub elements: Vec<Vec<u8>>,
}

impl BinaryMatrix {
    /// Return a binary matrix that all elements are zero
    pub fn new() -> Self {
        Self { elements: vec![] }
    }

    /// Return a binary matrix that all elements are zero
    ///
    /// # Arguments
    ///
    /// * `size_row` - size of row of the matrix
    /// * `size_col` - size of column of the matrix
    ///
    pub fn zero(size_row: usize, size_col: usize) -> Result<Self, String> {
        Ok(Self {
            elements: vec![vec![0; size_col]; size_row],
        })
    }

    /// Return a binary matrix that all elements are one
    ///
    /// # Arguments
    ///
    /// * `size_row` - size of row of the matrix
    /// * `size_col` - size of column of the matrix
    ///
    pub fn one(size_row: usize, size_col: usize) -> Result<Self, String> {
        Ok(Self {
            elements: vec![vec![1; size_col]; size_row],
        })
    }

    /// Return an identity binary matrix
    ///
    /// # Arguments
    ///
    /// * `size` - size of the matrix
    ///
    pub fn identity(size: usize) -> Result<Self, String> {
        if size == 0 {
            return Err("The size must be greater than 0.".to_string());
        }
        let mut bmat = Self::zero(size, size)?;
        for i in 0..size {
            bmat.elements[i][i] = 1;
        }
        Ok(bmat)
    }

    /// Check if the binary matrix is valid or not
    ///
    pub fn is_valid(&self) -> bool {
        if self.elements.len() == 0 {
            return false;
        }

        let size_col = self.elements[0].len();
        for row_vec in &self.elements {
            if row_vec.len() != size_col {
                return false;
            }
        }

        for row_vec in &self.elements {
            for e in row_vec {
                if *e != 0 && *e != 1 {
                    return false;
                }
            }
        }

        true
    }

    /// Return a binary matrix with the specified matrix elements
    ///
    /// # Arguments
    ///
    /// * `elements` - binary matrix elements
    ///
    /// # Notes
    ///
    /// Each element of the matrix must be 0 or 1.
    ///
    pub fn any(elements: &Vec<Vec<u8>>) -> Result<Self, String> {
        let bmat = Self {
            elements: elements.clone(),
        };
        match bmat.is_valid() {
            true => Ok(bmat),
            false => Err("The binary matrix is invalid.".to_string()),
        }
    }

    /// Merge the multiple binary matrices in row direction
    ///
    /// # Arguments
    ///
    /// * `bmat_vec` - vector of the binary matrices
    ///
    /// # Notes
    ///
    /// Each binary matrix must have the same column size.
    ///
    pub fn merge_row(bmat_vec: &Vec<Self>) -> Result<Self, String> {
        let mut size_col = 0;
        for bm in bmat_vec {
            size_col = bm.size_of_col();
            if size_col == 0 {
                continue;
            } else {
                break;
            }
        }
        for bm in bmat_vec {
            if bm.size_of_col() != 0 && bm.size_of_col() != size_col {
                return Err("Each binary matrix must have the same column size.".to_string());
            }
        }

        let mut bmat = Self::new();
        for bm in bmat_vec {
            if bm.size_of_col() == 0 {
                continue;
            }
            for row_vec in bm.elements.iter() {
                bmat.elements.push((&row_vec).to_vec());
            }
        }
        Ok(bmat)
    }

    /// Merge the multiple binary matrices in column direction
    ///
    /// # Arguments
    ///
    /// * `bmat_vec` - vector of the binary matrices
    ///
    /// # Notes
    ///
    /// Each binary matrix must have the same row size.
    ///
    pub fn merge_col(bmat_vec: &Vec<Self>) -> Result<Self, String> {
        let mut size_row = 0;
        for bm in bmat_vec {
            size_row = bm.size_of_row();
            if size_row == 0 {
                continue;
            } else {
                break;
            }
        }
        for bm in bmat_vec {
            if bm.size_of_row() != 0 && bm.size_of_row() != size_row {
                return Err("Each binary matrix must have the same row size.".to_string());
            }
        }

        let mut bmat = Self::new();
        for bm in bmat_vec {
            if bm.size_of_row() == 0 {
                continue;
            } else if bmat.elements.len() == 0 {
                bmat = bm.clone();
                continue;
            }
            for i in 0..size_row {
                bmat.elements[i].extend(&bm.elements[i]);
            }
        }
        Ok(bmat)
    }

    /// Return the row size of the binary matrix
    pub fn size_of_row(&self) -> usize {
        self.elements.len()
    }

    /// Return the column size of the binary matrix
    pub fn size_of_col(&self) -> usize {
        if self.elements.len() == 0 {
            return 0;
        }
        self.elements[0].len()
    }

    /// Add the binary matrix
    ///
    /// # Arguments
    ///
    /// * `other` - binary matrix to add
    ///
    pub fn add(&self, other: &Self) -> Result<Self, String> {
        if !self.is_valid() || !other.is_valid() {
            return Err("The binary matrix is invalid.".to_string());
        }

        if (other.size_of_row() != self.size_of_row())
            || (other.size_of_row() != self.size_of_row())
        {
            return Err("Size of matrices are not match.".to_string());
        }

        let mut out = self.clone();
        for i in 0..self.size_of_row() {
            for j in 0..self.size_of_col() {
                out.elements[i][j] += other.elements[i][j];
                out.elements[i][j] %= 2;
            }
        }
        Ok(out)
    }

    /// Multiply the binary matrix
    ///
    /// # Arguments
    ///
    /// * `other` - binary matrix to multiply
    ///
    pub fn mul(&self, other: &Self) -> Result<Self, String> {
        if !self.is_valid() || !other.is_valid() {
            return Err("The binary matrix is invalid.".to_string());
        }

        if other.size_of_row() != self.size_of_col() {
            return Err("Size of matrices are not match.".to_string());
        }

        let mut out = Self::new();
        for i in 0..self.size_of_row() {
            let mut row: Vec<u8> = vec![];
            for j in 0..other.size_of_col() {
                let mut e: u8 = 0;
                for k in 0..self.size_of_col() {
                    e += self.elements[i][k] * other.elements[k][j];
                    e %= 2;
                }
                row.push(e);
            }
            out.elements.push(row);
        }
        Ok(out)
    }

    /// Transpose the binary matrix
    ///
    pub fn transpose(&self) -> Result<Self, String> {
        if !self.is_valid() {
            return Err("The binary matrix is invalid.".to_string());
        }

        let mut out = Self::new();
        for i in 0..self.size_of_col() {
            let mut row: Vec<u8> = vec![];
            for j in 0..self.size_of_row() {
                row.push(self.elements[j][i]);
            }
            out.elements.push(row);
        }
        Ok(out)
    }

    /// Return the inverse matrix of the binary matrix
    ///
    pub fn inverse(&self) -> Result<Self, String> {
        if !self.is_valid() {
            return Err("The binary matrix is invalid.".to_string());
        }

        let matrix = self.elements.clone();

        let n = matrix.len();
        if n == 0 {
            return Err("The input matrix must not be empty.".to_string());
        }

        // Check if it's a square matrix
        for row in matrix.iter() {
            if row.len() != n {
                return Err("The input matrix must be a square matrix.".to_string());
            }
        }

        // Create an enlarged matrix [A | I]
        let mut augmented_matrix: Vec<Vec<u8>> = Vec::with_capacity(n);
        for i in 0..n {
            let mut row = Vec::with_capacity(2 * n);
            for j in 0..n {
                row.push(matrix[i][j]);
            }
            for j in 0..n {
                row.push(if i == j { 1 } else { 0 }); // Add an identity matrix
            }
            augmented_matrix.push(row);
        }

        // Gauss-Jordan Method
        for i in 0..n {
            // Find a pivot
            let mut pivot_row = i;
            while pivot_row < n && augmented_matrix[pivot_row][i] == 0 {
                pivot_row += 1;
            }

            if pivot_row == n {
                // If no pivot is found, the singular matrix
                return Err(
                    "A matrix is a singular matrix, and there is no inverse matrix.".to_string(),
                );
            }

            // Swap pivot row to current row
            augmented_matrix.swap(i, pivot_row);

            // Clear other rows in the current column
            for j in 0..n {
                if i != j && augmented_matrix[j][i] == 1 {
                    for k in i..2 * n {
                        augmented_matrix[j][k] ^= augmented_matrix[i][k]; // Addition in GF(2) (XOR)
                    }
                }
            }
        }

        // Extract the inverse part on the right
        let mut inverse_matrix: Vec<Vec<u8>> = Vec::with_capacity(n);
        for i in 0..n {
            let mut row = Vec::with_capacity(n);
            for j in 0..n {
                row.push(augmented_matrix[i][j + n]);
            }
            inverse_matrix.push(row);
        }

        Self::any(&inverse_matrix)
    }

    /// Save the Binary Matrix to the file
    ///
    /// # Arguments
    ///
    /// * `bm_path` - Binary Matrix file path
    ///
    pub fn save(&self, bm_path: &str) -> Result<(), String> {
        let mut f = File::create(bm_path).expect("File not created.");
        f.write_all(self.to_string().as_bytes())
            .expect("Something went wrong writeing the file.");
        Ok(())
    }
}

impl fmt::Display for BinaryMatrix {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let size_row = self.size_of_row();
        let size_col = self.size_of_col();
        for i in 0..size_row {
            for j in 0..size_col {
                match j {
                    val if val == size_col - 1 => {
                        let _ = write!(f, "{}", self.elements[i][j]);
                    }
                    _ => {
                        let _ = write!(f, "{} ", self.elements[i][j]);
                    }
                }
            }
            let _ = write!(f, "\n");
        }
        write!(f, "")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_success() {
        let bmat = BinaryMatrix::new();
        assert_eq!(bmat.elements, Vec::<Vec<u8>>::new());
    }

    #[test]
    fn zero_success() {
        let bmat = BinaryMatrix::zero(2, 3).unwrap();
        assert_eq!(bmat.elements[0][0], 0);
        assert_eq!(bmat.elements[0][1], 0);
        assert_eq!(bmat.elements[0][2], 0);
        assert_eq!(bmat.elements[1][0], 0);
        assert_eq!(bmat.elements[1][1], 0);
        assert_eq!(bmat.elements[1][2], 0);
    }

    #[test]
    fn identity_failure() {
        assert!(BinaryMatrix::identity(0).is_err());
    }

    #[test]
    fn identity_success() {
        let bmat = BinaryMatrix::identity(2).unwrap();
        assert_eq!(bmat.elements[0][0], 1);
        assert_eq!(bmat.elements[0][1], 0);
        assert_eq!(bmat.elements[1][0], 0);
        assert_eq!(bmat.elements[1][1], 1);
    }

    #[test]
    fn any_success() {
        let bmat = BinaryMatrix::any(&vec![vec![0, 1], vec![1, 1]]).unwrap();
        assert_eq!(bmat.elements[0][0], 0);
        assert_eq!(bmat.elements[0][1], 1);
        assert_eq!(bmat.elements[1][0], 1);
        assert_eq!(bmat.elements[1][1], 1);
    }

    #[test]
    fn any_failure() {
        assert!(BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0]]).is_err());
    }

    #[test]
    fn fmt_success() {
        let bmat = BinaryMatrix::any(&vec![vec![0, 1], vec![1, 1]]).unwrap();
        assert_eq!(bmat.to_string(), "0 1\n1 1\n");
    }

    #[test]
    fn merge_row_success() {
        let bmat_0 = BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0]]).unwrap();
        let bmat_1 = BinaryMatrix::any(&vec![vec![1, 0, 1]]).unwrap();
        let bmat_2 = BinaryMatrix::any(&vec![vec![1, 1, 1], vec![0, 0, 0]]).unwrap();
        let bmat_nil = BinaryMatrix::new();

        let bmat = BinaryMatrix::merge_row(&vec![bmat_0.clone()]).unwrap();
        assert_eq!(bmat.to_string(), "1 0 0\n0 1 0\n");

        let bmat = BinaryMatrix::merge_row(&vec![bmat_0.clone(), bmat_1.clone()]).unwrap();
        assert_eq!(bmat.to_string(), "1 0 0\n0 1 0\n1 0 1\n");

        let bmat =
            BinaryMatrix::merge_row(&vec![bmat_0.clone(), bmat_1.clone(), bmat_2.clone()]).unwrap();
        assert_eq!(bmat.to_string(), "1 0 0\n0 1 0\n1 0 1\n1 1 1\n0 0 0\n");

        let bmat = BinaryMatrix::merge_row(&vec![
            bmat_nil.clone(),
            bmat_0.clone(),
            bmat_nil.clone(),
            bmat_1.clone(),
        ])
        .unwrap();
        assert_eq!(bmat.to_string(), "1 0 0\n0 1 0\n1 0 1\n");
    }

    #[test]
    fn merge_row_failure() {
        let bmat_0 = BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0]]).unwrap();
        let bmat_1 = BinaryMatrix::any(&vec![vec![1, 0, 1]]).unwrap();
        let bmat_2 = BinaryMatrix::any(&vec![vec![1, 1], vec![0, 0]]).unwrap();
        assert!(
            BinaryMatrix::merge_row(&vec![bmat_0.clone(), bmat_1.clone(), bmat_2.clone()]).is_err()
        );
    }

    #[test]
    fn merge_col_success() {
        let bmat_0 = BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0]]).unwrap();
        let bmat_1 = BinaryMatrix::any(&vec![vec![1, 0], vec![0, 1]]).unwrap();
        let bmat_2 = BinaryMatrix::any(&vec![vec![1, 1], vec![0, 0]]).unwrap();
        let bmat_nil = BinaryMatrix::new();

        let bmat = BinaryMatrix::merge_col(&vec![bmat_0.clone()]).unwrap();
        assert_eq!(bmat.to_string(), "1 0 0\n0 1 0\n");

        let bmat = BinaryMatrix::merge_col(&vec![bmat_0.clone(), bmat_1.clone()]).unwrap();
        assert_eq!(bmat.to_string(), "1 0 0 1 0\n0 1 0 0 1\n");

        let bmat =
            BinaryMatrix::merge_col(&vec![bmat_0.clone(), bmat_1.clone(), bmat_2.clone()]).unwrap();
        assert_eq!(bmat.to_string(), "1 0 0 1 0 1 1\n0 1 0 0 1 0 0\n");

        let bmat = BinaryMatrix::merge_col(&vec![
            bmat_nil.clone(),
            bmat_0.clone(),
            bmat_nil.clone(),
            bmat_1.clone(),
        ])
        .unwrap();
        assert_eq!(bmat.to_string(), "1 0 0 1 0\n0 1 0 0 1\n");
    }

    #[test]
    fn merge_col_failure() {
        let bmat_0 = BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0]]).unwrap();
        let bmat_1 = BinaryMatrix::any(&vec![vec![1, 0], vec![0, 1]]).unwrap();
        let bmat_2 = BinaryMatrix::any(&vec![vec![1, 1], vec![0, 0], vec![1, 0]]).unwrap();
        assert!(
            BinaryMatrix::merge_col(&vec![bmat_0.clone(), bmat_1.clone(), bmat_2.clone()]).is_err()
        );
    }

    #[test]
    fn add_success() {
        // example 1
        let bmat = BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]]).unwrap();
        let other = BinaryMatrix::any(&vec![vec![1, 0, 1], vec![0, 1, 0], vec![0, 1, 0]]).unwrap();
        let expect = BinaryMatrix::any(&vec![vec![0, 0, 1], vec![0, 0, 0], vec![0, 1, 1]]).unwrap();
        assert_eq!(bmat.add(&other).unwrap(), expect);

        // example 2
        let bmat = BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 0, 1]]).unwrap();
        let other = BinaryMatrix::any(&vec![vec![1, 0, 1], vec![0, 1, 0]]).unwrap();
        let expect = BinaryMatrix::any(&vec![vec![0, 0, 1], vec![0, 1, 1]]).unwrap();
        assert_eq!(bmat.add(&other).unwrap(), expect);

        // example 3
        let bmat = BinaryMatrix::any(&vec![vec![1, 0], vec![0, 1], vec![0, 0]]).unwrap();
        let other = BinaryMatrix::any(&vec![vec![1, 0], vec![0, 1], vec![0, 1]]).unwrap();
        let expect = BinaryMatrix::any(&vec![vec![0, 0], vec![0, 0], vec![0, 1]]).unwrap();
        assert_eq!(bmat.add(&other).unwrap(), expect);
    }

    #[test]
    fn add_failure() {
        // example 1
        let bmat = BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]]).unwrap();
        let other = BinaryMatrix::any(&vec![vec![1, 0, 1], vec![0, 1, 0]]).unwrap();
        assert!(bmat.add(&other).is_err());

        // example 2
        let bmat = BinaryMatrix::any(&vec![vec![1, 0], vec![0, 1], vec![0, 0]]).unwrap();
        let other = BinaryMatrix::any(&vec![vec![1, 0, 1], vec![0, 1, 0]]).unwrap();
        assert!(bmat.add(&other).is_err());
    }

    #[test]
    fn mul_success() {
        // example 1
        let bmat = BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]]).unwrap();
        let other = BinaryMatrix::any(&vec![vec![1, 0, 1], vec![0, 1, 0], vec![0, 1, 0]]).unwrap();
        let expect = BinaryMatrix::any(&vec![vec![1, 0, 1], vec![0, 1, 0], vec![0, 1, 0]]).unwrap();
        assert_eq!(bmat.mul(&other).unwrap(), expect);

        // example 2
        let bmat = BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0]]).unwrap();
        let other = BinaryMatrix::any(&vec![vec![1, 0], vec![0, 1], vec![0, 1]]).unwrap();
        let expect = BinaryMatrix::any(&vec![vec![1, 0], vec![0, 1]]).unwrap();
        assert_eq!(bmat.mul(&other).unwrap(), expect);

        // example 3
        let bmat = BinaryMatrix::any(&vec![vec![1, 0], vec![0, 1], vec![0, 0]]).unwrap();
        let other = BinaryMatrix::any(&vec![vec![1, 0, 1], vec![0, 1, 0]]).unwrap();
        let expect = BinaryMatrix::any(&vec![vec![1, 0, 1], vec![0, 1, 0], vec![0, 0, 0]]).unwrap();
        assert_eq!(bmat.mul(&other).unwrap(), expect);
    }

    #[test]
    fn mul_failure() {
        // example 1
        let bmat = BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]]).unwrap();
        let other = BinaryMatrix::any(&vec![vec![1, 0, 1], vec![0, 1, 0]]).unwrap();
        assert!(bmat.mul(&other).is_err());

        // example 2
        let bmat = BinaryMatrix::any(&vec![vec![1, 1], vec![0, 0], vec![0, 1]]).unwrap();
        let other = BinaryMatrix::any(&vec![vec![1, 0], vec![0, 1], vec![0, 0]]).unwrap();
        assert!(bmat.mul(&other).is_err());
    }

    #[test]
    fn transpose_success() {
        // example 1
        let bmat = BinaryMatrix::any(&vec![vec![1, 0, 1], vec![0, 1, 1], vec![0, 0, 1]]).unwrap();
        let expect = BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0], vec![1, 1, 1]]).unwrap();
        assert_eq!(bmat.transpose().unwrap(), expect);

        // example 2
        let bmat = BinaryMatrix::any(&vec![vec![1, 0], vec![0, 1], vec![1, 0]]).unwrap();
        let expect = BinaryMatrix::any(&vec![vec![1, 0, 1], vec![0, 1, 0]]).unwrap();
        assert_eq!(bmat.transpose().unwrap(), expect);
    }

    #[test]
    fn inverse_success() {
        // identity matrix
        let bmat = BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]]).unwrap();
        let bmat_inv =
            BinaryMatrix::any(&vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]]).unwrap();
        assert_eq!(bmat.inverse().unwrap(), bmat_inv);

        // bravyi-kitaev encoding matrix
        let bmat = BinaryMatrix::any(&vec![
            vec![1, 0, 0, 0],
            vec![1, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![1, 1, 1, 1],
        ])
        .unwrap();
        let bmat_inv = BinaryMatrix::any(&vec![
            vec![1, 0, 0, 0],
            vec![1, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![0, 1, 1, 1],
        ])
        .unwrap();
        assert_eq!(bmat.inverse().unwrap(), bmat_inv);

        // large invertible matrix
        let bmat = BinaryMatrix::any(&vec![
            vec![1, 0, 0, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![1, 1, 1, 1],
        ])
        .unwrap();
        let bmat_inv = BinaryMatrix::any(&vec![
            vec![1, 0, 0, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![1, 1, 1, 1],
        ])
        .unwrap();
        assert_eq!(bmat.inverse().unwrap(), bmat_inv);
    }

    #[test]
    fn inverse_failure() {
        // empty
        let bmat = BinaryMatrix::new();
        assert!(bmat.inverse().is_err());

        // non square matrix
        let bmat = BinaryMatrix::any(&vec![vec![1, 0, 1], vec![0, 1, 1]]).unwrap();
        assert!(bmat.inverse().is_err());
    }
}
