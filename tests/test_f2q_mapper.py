import pytest
import random

from yaf2q.ternary_tree_spec import TernaryTreeSpec
from yaf2q.f2q_mapper import F2QMapper


def test_constructor():

    # named mapper: jordan-wigner
    f2q_mapper = F2QMapper(
        kind = "jordan-wigner",
        num_qubits = 4,
    )
    assert f2q_mapper.kind == "jordan-wigner"
    assert f2q_mapper.num_qubits == 4
    assert f2q_mapper.ttspec is None
    assert f2q_mapper.to_string() == "named mapper - jordan-wigner (num_qubits:4)"

    # maned mapper: parity
    f2q_mapper = F2QMapper(
        kind = "parity",
        num_qubits = 4,
    )
    assert f2q_mapper.kind == "parity"
    assert f2q_mapper.num_qubits == 4
    assert f2q_mapper.ttspec is None
    assert f2q_mapper.to_string() == "named mapper - parity (num_qubits:4)"

    # named mapper: bravyi-kitaev
    f2q_mapper = F2QMapper(
        kind = "bravyi-kitaev",
        num_qubits = 4,
    )
    assert f2q_mapper.kind == "bravyi-kitaev"
    assert f2q_mapper.num_qubits == 4
    assert f2q_mapper.ttspec is None
    assert f2q_mapper.to_string() == "named mapper - bravyi-kitaev (num_qubits:4)"

    # ternary tree mapper
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec(
            indices = [1, 0, 3, 2],
            edges = {1: (0, 'X'), 2: (1, 'X'), 3: (0, 'Y')},
        )
    )
    assert f2q_mapper.kind is None
    assert f2q_mapper.num_qubits == 4
    assert f2q_mapper.ttspec == TernaryTreeSpec(
        indices = [1, 0, 3, 2],
        edges = {1: (0, 'X'), 2: (1, 'X'), 3: (0, 'Y')},
    )
    assert f2q_mapper.to_string() == "ternary tree mapper - indices:[1, 0, 3, 2], edges:{1: (0, 'X'), 2: (1, 'X'), 3: (0, 'Y')"


def test_eigenvalue_H2(fermion_hamiltonian_H2):

    num_qubits = fermion_hamiltonian_H2.one_body_tensor.shape[0]

    # jordan-wigner
    f2q_mapper = F2QMapper(kind="jordan-wigner", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2).eigenvalues(num=1)[0]
    expect = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2, method="of").eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)

    # parity
    f2q_mapper = F2QMapper(kind="parity", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)

    # bravyi-kiaev
    f2q_mapper = F2QMapper(kind="bravyi-kitaev", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)

    # ternary-tree
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec(
            indices = [1, 0, 3, 2],
            edges = {1: (0, 'X'), 2: (1, 'X'), 3: (0, 'Y')},
        )
    )
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)


def test_eigenvalue_H3(fermion_hamiltonian_H3):

    num_qubits = fermion_hamiltonian_H3.one_body_tensor.shape[0]

    # jordan-wigner
    f2q_mapper = F2QMapper(kind="jordan-wigner", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H3).eigenvalues(num=1)[0]
    expect = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H3, method="of").eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)

    # parity
    f2q_mapper = F2QMapper(kind="parity", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H3).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)
    
    # ternary-tree
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec(
            indices = [5, 3, 4, 2, 0, 1],
            edges = {1: (0, 'X'), 2: (1, 'Y'), 3: (0, 'Z'), 4: (1, 'X'), 5: (3, 'Y')}
        )
    )
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H3).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)


def test_eigenvalue_H4(fermion_hamiltonian_H4):

    num_qubits = fermion_hamiltonian_H4.one_body_tensor.shape[0]

    # jordan-wigner
    f2q_mapper = F2QMapper(kind="jordan-wigner", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H4).eigenvalues(num=1)[0]
    expect = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H4, method="of").eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)

    # parity
    f2q_mapper = F2QMapper(kind="parity", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H4).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)

    # bravyi-kiaev
    f2q_mapper = F2QMapper(kind="bravyi-kitaev", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H4).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)

    # ternary-tree
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec(
            indices = [2, 4, 0, 6, 7, 3, 5, 1],
            edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'X'), 4: (2, 'Z'), 5: (0, 'Y'), 6: (5, 'Y'), 7: (6, 'Y')}
        )
    )
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H4).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)
    

def test_eigenvalue_LiH(fermion_hamiltonian_LiH):

    num_qubits = fermion_hamiltonian_LiH.one_body_tensor.shape[0]

    # jordan-wigner
    f2q_mapper = F2QMapper(kind="jordan-wigner", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_LiH).eigenvalues(num=1)[0]
    expect = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_LiH, method="of").eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)

    # parity
    f2q_mapper = F2QMapper(kind="parity", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_LiH).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)

    # ternary-tree
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec(
            indices = [8, 5, 3, 9, 4, 10, 2, 7, 0, 11, 6, 1],
            edges = {1: (0, 'X'), 2: (1, 'Y'), 3: (0, 'Z'), 4: (1, 'X'), 5: (3, 'Y'), 6: (4, 'X'),
                     7: (0, 'Y'), 8: (7, 'Z'), 9: (5, 'Y'), 10: (7, 'X'), 11: (10, 'Z'),}
        )
    )
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_LiH).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)


def test_eigenvalue_BeH2(fermion_hamiltonian_BeH2):

    num_qubits = fermion_hamiltonian_BeH2.one_body_tensor.shape[0]

    # jordan-wigner
    f2q_mapper = F2QMapper(kind="jordan-wigner", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_BeH2).eigenvalues(num=1)[0]
    expect = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_BeH2, method="of").eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)

    # parity
    f2q_mapper = F2QMapper(kind="parity", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_BeH2).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)

    # ternary-tree.random (seed = 1)
    random.seed(1)
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec.random(num_qubits)
    )
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_BeH2).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)


def test_eigenvalue_H2O(fermion_hamiltonian_H2O):

    num_qubits = fermion_hamiltonian_H2O.one_body_tensor.shape[0]

    # jordan-wigner
    f2q_mapper = F2QMapper(kind="jordan-wigner", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2O).eigenvalues(num=1)[0]
    expect = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2O, method="of").eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)
    
    # parity
    f2q_mapper = F2QMapper(kind="parity", num_qubits=num_qubits)
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2O).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)
    
    # ternary-tree.random (seed = 1)
    random.seed(1)
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec.random(num_qubits)
    )
    actual = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2O).eigenvalues(num=1)[0]
    assert actual == pytest.approx(expect, abs=0.0001)


def test_fock_to_qubit_state_H2(fermion_hamiltonian_H2):

    num_qubits = fermion_hamiltonian_H2.one_body_tensor.shape[0]

    # jordan-wigner
    f2q_mapper = F2QMapper(kind="jordan-wigner", num_qubits=num_qubits)
    actual = f2q_mapper.fock_to_qubit_state([1, 1, 0, 0])
    expect = [1, 1, 0, 0]
    assert actual == expect

    # parity
    f2q_mapper = F2QMapper(kind="parity", num_qubits=num_qubits)
    actual = f2q_mapper.fock_to_qubit_state([1, 1, 0, 0])
    expect = [1, 0, 0, 0]
    assert actual == expect

    # bravyi-kitaev
    f2q_mapper = F2QMapper(kind="bravyi-kitaev", num_qubits=num_qubits)
    actual = f2q_mapper.fock_to_qubit_state([1, 1, 0, 0])
    expect = [1, 0, 0, 0]
    assert actual == expect
    
    # ternary-tree
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec(
            indices = [1, 0, 3, 2],
            edges = {1: (0, 'X'), 2: (1, 'X'), 3: (0, 'Y')},
        )
    )
    actual = f2q_mapper.fock_to_qubit_state([1, 1, 0, 0])
    expect = [1, 0, 0, 0]
    assert actual == expect


def test_fock_to_qubit_state_H3(fermion_hamiltonian_H3):

    num_qubits = fermion_hamiltonian_H3.one_body_tensor.shape[0]

    # jordan-wigner
    f2q_mapper = F2QMapper(kind="jordan-wigner", num_qubits=num_qubits)
    actual = f2q_mapper.fock_to_qubit_state([1, 1, 1, 0, 0, 0])
    expect = [1, 1, 1, 0, 0, 0]
    assert actual == expect

    # parity
    f2q_mapper = F2QMapper(kind="parity", num_qubits=num_qubits)
    actual = f2q_mapper.fock_to_qubit_state([1, 1, 1, 0, 0, 0])
    expect = [1, 0, 1, 1, 1, 1]
    assert actual == expect

    # ternary-tree
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec(
            indices = [5, 3, 4, 2, 0, 1],
            edges = {1: (0, 'X'), 2: (1, 'Y'), 3: (0, 'Z'), 4: (1, 'X'), 5: (3, 'Y')}
        )
    )
    actual = f2q_mapper.fock_to_qubit_state([1, 1, 1, 0, 0, 0])
    expect = [1, 1, 1, 1, 0, 0]
    assert actual == expect


def test_fock_to_qubit_state_H4(fermion_hamiltonian_H4):

    num_qubits = fermion_hamiltonian_H4.one_body_tensor.shape[0]

    # jordan-wigner
    f2q_mapper = F2QMapper(kind="jordan-wigner", num_qubits=num_qubits)
    actual = f2q_mapper.fock_to_qubit_state([1, 1, 1, 1, 0, 0, 0, 0])
    expect = [1, 1, 1, 1, 0, 0, 0, 0]
    assert actual == expect

    # parity
    f2q_mapper = F2QMapper(kind="parity", num_qubits=num_qubits)
    actual = f2q_mapper.fock_to_qubit_state([1, 1, 1, 1, 0, 0, 0, 0])
    expect = [1, 0, 1, 0, 0, 0, 0, 0]
    assert actual == expect

    # bravyi-kitaev
    f2q_mapper = F2QMapper(kind="bravyi-kitaev", num_qubits=num_qubits)
    actual = f2q_mapper.fock_to_qubit_state([1, 1, 1, 1, 0, 0, 0, 0])
    expect = [1, 0, 1, 0, 0, 0, 0, 0]    
    assert actual == expect

    # ternary-tree
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec(
            indices = [2, 4, 0, 6, 7, 3, 5, 1],
            edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'X'), 4: (2, 'Z'), 5: (0, 'Y'), 6: (5, 'Y'), 7: (6, 'Y')}
        )
    )
    actual = f2q_mapper.fock_to_qubit_state([1, 1, 1, 1, 0, 0, 0, 0])
    expect = [1, 0, 0, 1, 0, 0, 0, 0]
    assert actual == expect


def test_qubit_to_fock_state_H2(fermion_hamiltonian_H2):

    num_qubits = fermion_hamiltonian_H2.one_body_tensor.shape[0]
    fock_state = [1, 1, 0, 0]

    # jordan-wigner
    f2q_mapper = F2QMapper(kind="jordan-wigner", num_qubits=num_qubits)
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    actual = f2q_mapper.qubit_to_fock_state(qubit_state)
    assert actual == fock_state

    # parity
    f2q_mapper = F2QMapper(kind="parity", num_qubits=num_qubits)
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    actual = f2q_mapper.qubit_to_fock_state(qubit_state)
    assert actual == fock_state

    # bravyi-kitaev
    f2q_mapper = F2QMapper(kind="bravyi-kitaev", num_qubits=num_qubits)
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    actual = f2q_mapper.qubit_to_fock_state(qubit_state)
    assert actual == fock_state
    
    # ternary-tree
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec(
            indices = [1, 0, 3, 2],
            edges = {1: (0, 'X'), 2: (1, 'X'), 3: (0, 'Y')},
        )
    )
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    actual = f2q_mapper.qubit_to_fock_state(qubit_state)
    assert actual == fock_state


def test_qubit_to_fock_state_H3(fermion_hamiltonian_H3):

    num_qubits = fermion_hamiltonian_H3.one_body_tensor.shape[0]
    fock_state = [1, 1, 1, 0, 0, 0]

    # jordan-wigner
    f2q_mapper = F2QMapper(kind="jordan-wigner", num_qubits=num_qubits)
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    actual = f2q_mapper.qubit_to_fock_state(qubit_state)
    assert actual == fock_state

    # parity
    f2q_mapper = F2QMapper(kind="parity", num_qubits=num_qubits)
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    actual = f2q_mapper.qubit_to_fock_state(qubit_state)
    assert actual == fock_state

    # ternary-tree
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec(
            indices = [5, 3, 4, 2, 0, 1],
            edges = {1: (0, 'X'), 2: (1, 'Y'), 3: (0, 'Z'), 4: (1, 'X'), 5: (3, 'Y')}
        )
    )
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    actual = f2q_mapper.qubit_to_fock_state(qubit_state)
    assert actual == fock_state


def test_qubit_to_fock_state_H4(fermion_hamiltonian_H4):

    num_qubits = fermion_hamiltonian_H4.one_body_tensor.shape[0]
    fock_state = [1, 1, 1, 1, 0, 0, 0, 0]

    # jordan-wigner
    f2q_mapper = F2QMapper(kind="jordan-wigner", num_qubits=num_qubits)
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    actual = f2q_mapper.qubit_to_fock_state(qubit_state)
    assert actual == fock_state

    # parity
    f2q_mapper = F2QMapper(kind="parity", num_qubits=num_qubits)
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    actual = f2q_mapper.qubit_to_fock_state(qubit_state)
    assert actual == fock_state

    # bravyi-kitaev
    f2q_mapper = F2QMapper(kind="bravyi-kitaev", num_qubits=num_qubits)
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    actual = f2q_mapper.qubit_to_fock_state(qubit_state)
    assert actual == fock_state
    
    # ternary-tree
    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec(
            indices = [2, 4, 0, 6, 7, 3, 5, 1],
            edges = {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'X'), 4: (2, 'Z'), 5: (0, 'Y'), 6: (5, 'Y'), 7: (6, 'Y')}
        )
    )
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    actual = f2q_mapper.qubit_to_fock_state(qubit_state)
    assert actual == fock_state
