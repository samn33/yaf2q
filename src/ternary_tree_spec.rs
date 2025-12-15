//! Module for Ternary Tree Specification

use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::prelude::*;

use crate::binary_matrix::BinaryMatrix;
use crate::util;

/// Struct of Ternary Tree Specification
///
/// # Fields
///
/// * `indices` - indices vector
/// * `parent_map` - parent index map
/// * `mask` - mask of parent indices
///
#[derive(Debug, Clone, PartialEq)]
pub struct TernayTreeSpec {
    pub indices: Vec<usize>,
    pub parent_map: HashMap<usize, usize>,
    mask: Vec<u8>,
}

impl TernayTreeSpec {
    /// Return a Ternary Tree 
    ///
    pub fn new() -> Self {
        Self {
            indices: vec![],
            parent_map: HashMap::new(),
            mask: Vec::new(),
        }
    }

    /// Create a TernaryTreeSpec
    ///
    /// # Arguments
    ///
    /// * `indices` - indices vector
    /// * `parent_map` - parent index map
    ///
    pub fn create(
        indices: &Vec<usize>,
        parent_map: &HashMap<usize, usize>,
    ) -> Result<Self, String> {
        let num_indices = indices.len();
        let mut mask: Vec<u8> = vec![0; (num_indices - 1) * 3];
        for (_, v) in parent_map.iter() {
            mask[*v] = 1;
        }

        let ttspec = Self {
            indices: indices.clone(),
            parent_map: parent_map.clone(),
            mask: mask,
        };
        match ttspec.is_valid() {
            true => Ok(ttspec),
            false => Err("The Ternary Tree Generator is invalid.".to_string()),
        }
    }

    /// Create a TernaryTreeSpec randomly
    ///
    /// # Arguments
    ///
    /// * `num_indices` - number of indices
    /// * `seed_value` - seed value
    ///
    #[allow(deprecated)]
    pub fn random(num_indices: usize, seed_value: Option<u64>) -> Result<Self, String> {
        let mut rng = match seed_value {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_os_rng(),
        };

        let mut mask: Vec<u8> = vec![0; (num_indices - 1) * 3];

        let mut indices: Vec<usize> = vec![];
        for i in 0..num_indices {
            indices.push(i);
        }
        indices.shuffle(&mut rng);

        let mut parent_map = HashMap::<usize, usize>::new();
        for i in 1..num_indices {
            loop {
                let r = rng.gen_range(0..i * 3);
                if mask[r] == 1 {
                    continue;
                } else {
                    parent_map.insert(i, r);
                    mask[r] = 1;
                    break;
                }
            }
        }

        let ttspec = Self {
            indices: indices.clone(),
            parent_map: parent_map.clone(),
            mask: mask.clone(),
        };

        match ttspec.is_valid() {
            true => {}
            false => {
                return Err("The Ternary Tree Generator is invalid.".to_string());
            }
        };

        Ok(ttspec)
    }

    /// Check if the TernaryTreeSpec is valid or not
    ///
    pub fn is_valid(&self) -> bool {
        let num_indices = self.indices.len();

        // check duplications of the indices
        let set: HashSet<_> = self.indices.iter().collect();
        if num_indices != set.len() {
            return false;
        }

        // check length of the parent_map
        if self.parent_map.len() != num_indices - 1 {
            return false;
        }

        // validation of the parent_map
        let mut values = HashSet::<usize>::new();
        for i in 1..num_indices {
            match self.parent_map.get(&i) {
                Some(v) => {
                    if *v >= 3 * i {
                        return false;
                    } else if self.mask[*v] == 0 {
                        return false;
                    }
                    values.insert(*v);
                }
                None => {
                    return false;
                }
            }
        }
        if values.len() != num_indices - 1 {
            return false;
        }
        true
    }

    /// Swap two randomly selected elements in indices
    ///
    pub fn swap_indices(&self, rng: &mut StdRng) -> Self {
        let mut ttspec = self.clone();
        let num_indices = self.indices.len();
        let a = rng.random_range(0..num_indices);
        let mut b;
        loop {
            b = rng.random_range(0..num_indices);
            if b == a {
                continue;
            } else {
                break;
            }
        }

        (ttspec.indices[a], ttspec.indices[b]) = (ttspec.indices[b], ttspec.indices[a]);
        ttspec
    }

    /// Randomly replace randomly selected subnodes with another parent node
    ///
    pub fn change_parent(&self, rng: &mut StdRng) -> Self {
        let mut ttspec = self.clone();
        let num_indices = self.indices.len();
        let mut child;

        // select child randomly
        let mut cnt = 0;
        loop {
            child = rng.random_range(1..num_indices - 1);

            if cnt > 100 {
                println!("child = {}", child);
                println!("mask = {:?}", ttspec.mask);
                panic!("Infinete loop!");
            }
            cnt += 1;

            let mut mask_sum: usize = 0;
            for b in ttspec.mask[0..child * 3].iter() {
                mask_sum += *b as usize;
            }
            if mask_sum == child * 3 {
                continue;
            } else {
                break;
            }
        }

        // select parent randomly
        let parent_org = ttspec.parent_map[&child];
        let mut parent;
        let mut cnt = 0;
        loop {
            parent = rng.random_range(0..child * 3);

            if cnt > 100 {
                println!(
                    "child = {}, parent_org = {}, parent = {}",
                    child, parent_org, parent
                );
                println!("mask = {:?}", ttspec.mask);
                panic!("Infinete loop!");
            }
            cnt += 1;

            if parent == parent_org || ttspec.mask[parent] == 1 {
                continue;
            } else {
                ttspec.parent_map.insert(child, parent);
                ttspec.mask[parent_org] = 0;
                ttspec.mask[parent] = 1;
                break;
            }
        }
        ttspec
    }

    /// Get the TernaryTreeSpec from the text
    ///
    /// # Arguments
    ///
    /// * `contents` - text describing the TernaryTreeSpec
    ///
    pub fn from_str(contents: &str) -> Result<Self, String> {
        // indices, parent_map
        let mut indices = Vec::<usize>::new();
        let mut parent_map = HashMap::<usize, usize>::new();
        for (i, line) in contents
            .lines()
            .filter(|line| !line.trim().is_empty())
            .enumerate()
        {
            match i {
                0 => {
                    for x in line.trim().split_whitespace().map(|s| s.parse::<usize>()) {
                        indices.push(x.unwrap());
                    }
                }
                _ => {
                    let mut kv: Vec<usize> = vec![];
                    for x in line.trim().split_whitespace().map(|s| s.parse::<usize>()) {
                        kv.push(x.unwrap());
                    }
                    parent_map.insert(kv[0], kv[1]);
                }
            };
        }

        // mask
        let num_indices = indices.len();
        let mut mask: Vec<u8> = vec![0; (num_indices - 1) * 3];
        for (_, v) in parent_map.iter() {
            mask[*v] = 1;
        }

        let ttspec = Self {
            indices: indices,
            parent_map: parent_map,
            mask: mask,
        };
        Ok(ttspec)
    }

    /// Load a TernaryTreeSpec file
    ///
    /// # Arguments
    ///
    /// * `ttspec_path` - File path of the Ternary Tree Generator
    ///
    pub fn load(ttspec_path: &str) -> Result<Self, String> {
        let mut f = File::open(ttspec_path).expect("File not found.");
        let mut contents = String::new();
        f.read_to_string(&mut contents)
            .expect("Something went wrong reading the file.");
        Self::from_str(&contents)
    }
}

/// Enumeration of ternary tree
///
/// # Fields
///
/// * `Node` - node of the ternary tree
/// * `Nil` - termination node of the ternary tree
///
#[derive(Debug, Clone, PartialEq)]
pub enum TernaryTree {
    Node {
        index: usize,
        up: Box<TernaryTree>,
        mid: Box<TernaryTree>,
        down: Box<TernaryTree>,
    },
    Nil,
}

impl TernaryTree {
    /// Create a TernaryTree for the Jordan-Wigner transformation
    ///
    /// # Arguments
    ///
    /// * `num` - number of qubits (= molecular orbitals)
    ///
    pub fn create_jordan_wigner(num: usize) -> Result<Self, String> {
        let mut tree = Self::Nil;
        for i in (0..num).rev() {
            tree = Self::Node {
                index: i,
                up: Box::new(Self::Nil),
                mid: Box::new(Self::Nil),
                down: Box::new(tree.clone()),
            };
        }
        Ok(tree)
    }

    /// Create a TernaryTree for the Parity transformation
    ///
    /// # Arguments
    ///
    /// * `num` - Number of qubits (= molecular orbitals)
    ///
    pub fn create_parity(num: usize) -> Result<Self, String> {
        let mut tree = Self::Nil;
        for i in 0..num {
            tree = Self::Node {
                index: i,
                up: Box::new(tree.clone()),
                mid: Box::new(Self::Nil),
                down: Box::new(Self::Nil),
            };
        }
        Ok(tree)
    }

    /// Create a ternary sub-tree for the Bravyi-Kitaev transformation
    ///
    /// # Arguments
    ///
    /// * `index` - index of the node
    /// * `step` - the difference between the indexes (index_up and index_down)
    ///            of the two leaf nodes (up and down) that this node is connected to.
    ///            Calculated as index_up = index + step, index_down = index - step.
    ///
    fn subtree_bravyi_kitaev(index: usize, mut step: usize) -> Self {
        if index % 2 == 0 {
            Self::Node {
                index: index,
                up: Box::new(Self::Nil),
                mid: Box::new(Self::Nil),
                down: Box::new(Self::Nil),
            }
        } else {
            step /= 2;
            let tree_up = Self::subtree_bravyi_kitaev(index - step, step);
            let tree_down = Self::subtree_bravyi_kitaev(index + step, step);
            Self::Node {
                index: index,
                up: Box::new(tree_up.clone()),
                mid: Box::new(Self::Nil),
                down: Box::new(tree_down.clone()),
            }
        }
    }

    /// Create a TernaryTree for the Bravyi-Kitaev transformation
    ///
    /// # Arguments
    ///
    /// * `num` - Number of qubits
    ///
    /// # Notes
    ///
    /// The num must be a power of 2 in the case of "bravyi-kitaev"
    ///
    pub fn create_bravyi_kitaev(num: usize) -> Result<Self, String> {
        if num == 0 {
            return Ok(Self::Nil);
        } else if num == 1 {
            return Ok(Self::Node {
                index: num - 1,
                up: Box::new(Self::Nil),
                mid: Box::new(Self::Nil),
                down: Box::new(Self::Nil),
            });
        } else if util::is_power_of_two(num) == false {
            return Err("The num must be a power of two.".to_string());
        }

        let subtree = Self::subtree_bravyi_kitaev(num / 2 - 1, num / 2);
        let tree = Self::Node {
            index: num - 1,
            up: Box::new(subtree.clone()),
            mid: Box::new(Self::Nil),
            down: Box::new(Self::Nil),
        };
        Ok(tree)
    }

    /// Create a Ternary sub-tree by TernaryTreeSpec
    ///
    /// # Arguments
    ///
    /// * `index` - Index of the node
    /// * `ttspec` - ternary tree specification
    ///
    fn subtree_ternary_tree(index_root: usize, ttspec: &TernayTreeSpec) -> Self {
        let num_indices = ttspec.indices.len();
        if index_root == num_indices - 1 {
            Self::Node {
                index: ttspec.indices[num_indices - 1],
                up: Box::new(Self::Nil),
                mid: Box::new(Self::Nil),
                down: Box::new(Self::Nil),
            }
        } else {
            let mut tree_up = Self::Nil;
            let mut tree_mid = Self::Nil;
            let mut tree_down = Self::Nil;
            for i in index_root..num_indices - 1 {
                if ttspec.parent_map[&(i + 1)] / 3 != index_root {
                    continue;
                }
                match ttspec.parent_map[&(i + 1)] % 3 {
                    0 => {
                        tree_up = Self::subtree_ternary_tree(i + 1, ttspec);
                    }
                    1 => {
                        tree_mid = Self::subtree_ternary_tree(i + 1, ttspec);
                    }
                    _ => {
                        tree_down = Self::subtree_ternary_tree(i + 1, ttspec);
                    }
                }
            }

            Self::Node {
                index: ttspec.indices[index_root],
                up: Box::new(tree_up.clone()),
                mid: Box::new(tree_mid.clone()),
                down: Box::new(tree_down.clone()),
            }
        }
    }

    /// Create a TernaryTree by the TernaryTreeSpec
    ///
    /// # Arguments
    ///
    /// * `ttspec` - ternary tree specification
    ///
    pub fn create_ternary_tree(ttspec: &TernayTreeSpec) -> Result<Self, String> {
        let _ = match ttspec.is_valid() {
            true => Ok(()),
            false => Err("The Ternary Tree Generator is invalid.".to_string()),
        };
        Ok(TernaryTree::subtree_ternary_tree(0, ttspec))
    }

    /// Return the encoding matrix from the TernaryTree
    ///
    /// # Arguments
    ///
    /// * `ttspec` - ternary tree specification
    ///
    fn ternary_tree_to_encoding_matrix(tree: &Self) -> BinaryMatrix {
        match tree {
            Self::Nil => BinaryMatrix::new(),
            Self::Node {
                index: _,
                up,
                mid,
                down,
            } => {
                let bmat_up = Self::ternary_tree_to_encoding_matrix(up);
                let bmat_mid = Self::ternary_tree_to_encoding_matrix(mid);
                let bmat_down = Self::ternary_tree_to_encoding_matrix(down);

                // The row of Ex
                let zero_0_1 = BinaryMatrix::zero(bmat_up.size_of_row(), 1).unwrap();
                let zero_0_2 =
                    BinaryMatrix::zero(bmat_up.size_of_row(), bmat_mid.size_of_col()).unwrap();
                let zero_0_3 =
                    BinaryMatrix::zero(bmat_up.size_of_row(), bmat_down.size_of_col()).unwrap();

                // The row of 1
                let one_1_0 = BinaryMatrix::one(1, bmat_up.size_of_col()).unwrap();
                let one_1_1 = BinaryMatrix::one(1, 1).unwrap();
                let one_1_2 = BinaryMatrix::one(1, bmat_mid.size_of_col()).unwrap();
                let zero_1_3 = BinaryMatrix::zero(1, bmat_down.size_of_col()).unwrap();

                // The row of Ey
                let zero_2_0 =
                    BinaryMatrix::zero(bmat_mid.size_of_row(), bmat_up.size_of_col()).unwrap();
                let zero_2_1 = BinaryMatrix::zero(bmat_mid.size_of_row(), 1).unwrap();
                let zero_2_3 =
                    BinaryMatrix::zero(bmat_mid.size_of_row(), bmat_down.size_of_col()).unwrap();

                // The row of Ez
                let zero_3_0 =
                    BinaryMatrix::zero(bmat_down.size_of_row(), bmat_up.size_of_col()).unwrap();
                let zero_3_1 = BinaryMatrix::zero(bmat_down.size_of_row(), 1).unwrap();
                let zero_3_2 =
                    BinaryMatrix::zero(bmat_down.size_of_row(), bmat_mid.size_of_col()).unwrap();

                // Merge rows
                let row_0 =
                    BinaryMatrix::merge_col(&vec![bmat_up, zero_0_1, zero_0_2, zero_0_3]).unwrap();
                let row_1 =
                    BinaryMatrix::merge_col(&vec![one_1_0, one_1_1, one_1_2, zero_1_3]).unwrap();
                let row_2 =
                    BinaryMatrix::merge_col(&vec![zero_2_0, zero_2_1, bmat_mid, zero_2_3]).unwrap();
                let row_3 = BinaryMatrix::merge_col(&vec![zero_3_0, zero_3_1, zero_3_2, bmat_down])
                    .unwrap();

                BinaryMatrix::merge_row(&vec![row_0, row_1, row_2, row_3]).unwrap()
            }
        }
    }

    /// Return the encoding matrix
    ///
    pub fn encoding_matrix(&self) -> Result<BinaryMatrix, String> {
        let bmat = Self::ternary_tree_to_encoding_matrix(&self);
        Ok(bmat)
    }
}

impl fmt::Display for TernayTreeSpec {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, idx) in self.indices.iter().enumerate() {
            if i == self.indices.len() - 1 {
                let _ = write!(f, "{}\n", idx);
            } else {
                let _ = write!(f, "{} ", idx);
            }
        }
        for (k, v) in &self.parent_map {
            let _ = writeln!(f, "{} {}", k, v);
        }
        write!(f, "")
    }
}

impl fmt::Display for TernaryTree {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self {
            Self::Node {
                index: _i,
                up: _u,
                mid: _m,
                down: _d,
            } => {
                let _ = write!(f, "[{} ", _i);
                let _ = write!(f, "{} ", _u);
                let _ = write!(f, "{} ", _m);
                let _ = write!(f, "{}", _d);
                let _ = write!(f, "]");
            }
            Self::Nil => {
                let _ = write!(f, "nil");
            }
        };
        write!(f, "")
    }
}

mod common_test_data {

    #[allow(dead_code)]
    pub fn ttspec_str() -> &'static str {
        "\n
3 2 1 0
1 0
2 3
3 6
"
    }
}

#[cfg(test)]
mod tests {
    use super::common_test_data::*;
    use super::*;

    #[test]
    fn ttspec_from_str_success() {
        let ttspec = TernayTreeSpec::from_str(ttspec_str()).unwrap();
        assert_eq!(ttspec.indices, vec![3, 2, 1, 0]);
        assert_eq!(
            ttspec.parent_map,
            vec![(1, 0), (2, 3), (3, 6)].into_iter().collect()
        );
    }

    #[test]
    fn create_success() {
        // Jordan-Wigner
        let tree = TernaryTree::create_jordan_wigner(4).unwrap();
        assert_eq!(
            tree.to_string(),
            "[0 nil nil [1 nil nil [2 nil nil [3 nil nil nil]]]]"
        );

        // Parity
        let tree = TernaryTree::create_parity(4).unwrap();
        assert_eq!(
            tree.to_string(),
            "[3 [2 [1 [0 nil nil nil] nil nil] nil nil] nil nil]"
        );

        // Bravyi-Kitaev
        let tree = TernaryTree::create_bravyi_kitaev(4).unwrap();
        assert_eq!(
            tree.to_string(),
            "[3 [1 [0 nil nil nil] nil [2 nil nil nil]] nil nil]"
        );
    }

    #[test]
    fn create_failure() {
        assert!(TernaryTree::create_bravyi_kitaev(10).is_err());
    }

    #[test]
    fn create_ternary_tree_success() {
        let indices = vec![0, 1, 2, 3];
        let parent_map: HashMap<usize, usize> = vec![(1, 0), (2, 1), (3, 4)].into_iter().collect();
        let ttspec = TernayTreeSpec::create(&indices, &parent_map).unwrap();
        let ttree = TernaryTree::create_ternary_tree(&ttspec).unwrap();
        assert_eq!(
            ttree.to_string(),
            "[0 [1 nil [3 nil nil nil] nil] [2 nil nil nil] nil]"
        );

        let indices = vec![3, 2, 1, 0];
        let parent_map: HashMap<usize, usize> = vec![(1, 0), (2, 1), (3, 4)].into_iter().collect();
        let ttspec = TernayTreeSpec::create(&indices, &parent_map).unwrap();
        let ttree = TernaryTree::create_ternary_tree(&ttspec).unwrap();
        assert_eq!(
            ttree.to_string(),
            "[3 [2 nil [0 nil nil nil] nil] [1 nil nil nil] nil]"
        );

        let indices = vec![0, 1, 2, 3];
        let parent_map: HashMap<usize, usize> = vec![(1, 0), (2, 1), (3, 2)].into_iter().collect();
        let ttspec = TernayTreeSpec::create(&indices, &parent_map).unwrap();
        let ttree = TernaryTree::create_ternary_tree(&ttspec).unwrap();
        assert_eq!(
            ttree.to_string(),
            "[0 [1 nil nil nil] [2 nil nil nil] [3 nil nil nil]]"
        );

        let indices = vec![0, 1, 2, 3];
        let parent_map: HashMap<usize, usize> = vec![(1, 2), (2, 1), (3, 0)].into_iter().collect();
        let ttspec = TernayTreeSpec::create(&indices, &parent_map).unwrap();
        let ttree = TernaryTree::create_ternary_tree(&ttspec).unwrap();
        assert_eq!(
            ttree.to_string(),
            "[0 [3 nil nil nil] [2 nil nil nil] [1 nil nil nil]]"
        );

        let indices = vec![0, 1, 2, 3];
        let parent_map: HashMap<usize, usize> = vec![(1, 0), (2, 4), (3, 8)].into_iter().collect();
        let ttspec = TernayTreeSpec::create(&indices, &parent_map).unwrap();
        let ttree = TernaryTree::create_ternary_tree(&ttspec).unwrap();
        assert_eq!(
            ttree.to_string(),
            "[0 [1 nil [2 nil nil [3 nil nil nil]] nil] nil nil]"
        );

        let indices = vec![0, 1, 2];
        let parent_map: HashMap<usize, usize> = vec![(1, 0), (2, 4)].into_iter().collect();
        let ttspec = TernayTreeSpec::create(&indices, &parent_map).unwrap();
        let ttree = TernaryTree::create_ternary_tree(&ttspec).unwrap();
        assert_eq!(ttree.to_string(), "[0 [1 nil [2 nil nil nil] nil] nil nil]");
    }

    #[test]
    fn create_ternary_tree_failure() {
        // duplications of the indices
        let indices = vec![3, 0, 1, 0];
        let parent_map: HashMap<usize, usize> = vec![(1, 0), (2, 1), (3, 4)].into_iter().collect();
        assert!(TernayTreeSpec::create(&indices, &parent_map).is_err());

        // the size of the parent_map
        let indices = vec![3, 2, 1, 0];
        let parent_map: HashMap<usize, usize> =
            vec![(1, 0), (2, 1), (3, 3), (4, 6)].into_iter().collect();
        assert!(TernayTreeSpec::create(&indices, &parent_map).is_err());

        // the size of the parent_map
        let indices = vec![3, 2, 1, 0];
        let parent_map: HashMap<usize, usize> = vec![(1, 0), (3, 4)].into_iter().collect();
        assert!(TernayTreeSpec::create(&indices, &parent_map).is_err());

        // constraints of the adj_mat (sum in the row direction must be 0 or 1)
        let indices = vec![3, 2, 1, 0];
        let parent_map: HashMap<usize, usize> = vec![(1, 0), (2, 1), (3, 1)].into_iter().collect();
        assert!(TernayTreeSpec::create(&indices, &parent_map).is_err());
    }

    #[test]
    fn endoding_matrix_success() {
        // Jordan-Wigner
        let ttree = TernaryTree::create_jordan_wigner(4).unwrap();
        let bmat = ttree.encoding_matrix().unwrap();
        assert_eq!(bmat.to_string(), "1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1\n");

        // Parity
        let ttree = TernaryTree::create_parity(4).unwrap();
        let bmat = ttree.encoding_matrix().unwrap();
        assert_eq!(bmat.to_string(), "1 0 0 0\n1 1 0 0\n1 1 1 0\n1 1 1 1\n");

        // Bravyi-Kitaev
        let ttree = TernaryTree::create_bravyi_kitaev(4).unwrap();
        let bmat = ttree.encoding_matrix().unwrap();
        assert_eq!(bmat.to_string(), "1 0 0 0\n1 1 0 0\n0 0 1 0\n1 1 1 1\n");

        // Jordan-Wigner from the Ternary Tree Generator
        let indices = vec![0, 1, 2, 3];
        let parent_map: HashMap<usize, usize> = vec![(1, 2), (2, 5), (3, 8)].into_iter().collect();
        let ttspec = TernayTreeSpec::create(&indices, &parent_map).unwrap();
        let ttree = TernaryTree::create_ternary_tree(&ttspec).unwrap();
        assert_eq!(
            ttree.to_string(),
            "[0 nil nil [1 nil nil [2 nil nil [3 nil nil nil]]]]"
        );
        let bmat = ttree.encoding_matrix().unwrap();
        assert_eq!(bmat.to_string(), "1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1\n");

        // Parity from the Ternary Tree Generator
        let indices = vec![3, 2, 1, 0];
        let parent_map: HashMap<usize, usize> = vec![(1, 0), (2, 3), (3, 6)].into_iter().collect();
        let ttspec = TernayTreeSpec::create(&indices, &parent_map).unwrap();
        let ttree = TernaryTree::create_ternary_tree(&ttspec).unwrap();
        assert_eq!(
            ttree.to_string(),
            "[3 [2 [1 [0 nil nil nil] nil nil] nil nil] nil nil]"
        );
        let bmat = ttree.encoding_matrix().unwrap();
        assert_eq!(bmat.to_string(), "1 0 0 0\n1 1 0 0\n1 1 1 0\n1 1 1 1\n");

        // Bravyi-Kitaev from the Ternary Tree Generator
        let indices = vec![3, 1, 0, 2];
        let parent_map: HashMap<usize, usize> = vec![(1, 0), (2, 3), (3, 5)].into_iter().collect();
        let ttspec = TernayTreeSpec::create(&indices, &parent_map).unwrap();
        let ttree = TernaryTree::create_ternary_tree(&ttspec).unwrap();
        assert_eq!(
            ttree.to_string(),
            "[3 [1 [0 nil nil nil] nil [2 nil nil nil]] nil nil]"
        );
        let bmat = ttree.encoding_matrix().unwrap();
        assert_eq!(bmat.to_string(), "1 0 0 0\n1 1 0 0\n0 0 1 0\n1 1 1 1\n");
    }

    #[test]
    fn random_success() {
        let ttspec = TernayTreeSpec::random(4, None).unwrap();
        assert!(ttspec.is_valid());
    }
}
