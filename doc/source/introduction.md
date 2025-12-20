Introduction
============

## Feature

yaf2q is a Python library for performing fermion-to-qubit mapping used in quantum chemical calculations. It provides the following two functions.

* **Performing fermion-to-qubit mapping**
  - It supports well-known fermion-qubit mapping methods such as the Jordan-Wigner, Parity, and Bravyi-Kitaev transformations, as well as a more generalized approach using ternary trees.

* **Optimization of fermion-to-qubit mapping**
  - It is possible to search for a ternary tree structure that minimizes a chosen index, such as the Pauli weight in the mapped Hamiltonian, or the depth of the quantum circuit required to execute a quantum algorithm like Quantum Phase Estimation (QPE), when converting a fermionic operator to a qubit operator.

## Install

```bash
$ git clone https://github.com/samn33/yaf2q.git
$ cd yaf2q
$ make install
```

## Uninstall

```bash
$ make uninstall
```

## Test

```bash
$ make test
```

## References

Papers about fermion-to-qubit mapping using ternary tree.

1. Haytham McDowall-Rose, Razin A. Shaikh, Lia Yeh, "From fermions to Qubits: A ZX-Calculus Perspective",
[arXiv:2505.06212](https://arxiv.org/abs/2505.06212)

## Requirements

- Linux (Ubuntu 24.04 LTS)
- Python 3.12
- Rust 1.88.0

## Licence

MIT

## Author

[Sam.N](http://github.com/samn33)
