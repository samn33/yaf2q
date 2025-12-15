import pytest

from yaf2q.ternary_tree_spec import TernaryTreeSpec
from yaf2q.f2q_mapper import F2QMapper
from yaf2q.optimizer.sa import SAParams, SA


def test_sa_H2(fermion_hamiltonian_H2):

    num_qubits = fermion_hamiltonian_H2.one_body_tensor.shape[0]

    def objective_func(ttspec: TernaryTreeSpec) -> float:
        f2q_mapper = F2QMapper(ttspec=ttspec)
        qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2)
        pauli_weights = qubit_hamiltonian.pauli_weights()
        return sum(pauli_weights) / len(pauli_weights)

    ttspec_opt = SA(
        num_qubits = num_qubits,
        objective_func = objective_func,
        params = SAParams(
            init_sampling = 10,
            num_steps = 3,
            cooling_factor = 1.0,
            seed = 123,
        ),
    ).optimize()
    
    f2q_mapper = F2QMapper(ttspec=ttspec_opt)
    qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2)

    pauli_weights = qubit_hamiltonian.pauli_weights()
    actual = sum(pauli_weights)/len(pauli_weights)
    expect = 2.1333333333333333
    assert actual == pytest.approx(expect, abs=1.0e-8)


def test_sa_H3(fermion_hamiltonian_H3):

    num_qubits = fermion_hamiltonian_H3.one_body_tensor.shape[0]

    def objective_func(ttspec: TernaryTreeSpec) -> float:
        f2q_mapper = F2QMapper(ttspec=ttspec)
        qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H3)
        pauli_weights = qubit_hamiltonian.pauli_weights()
        return sum(pauli_weights) / len(pauli_weights)

    ttspec_opt = SA(
        num_qubits = num_qubits,
        objective_func = objective_func,
        params = SAParams(
            init_sampling = 10,
            num_steps = 3,
            cooling_factor = 1.0,
            seed = 123,
        ),
    ).optimize()
    
    f2q_mapper = F2QMapper(ttspec=ttspec_opt)
    qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H3)

    pauli_weights = qubit_hamiltonian.pauli_weights()
    actual = sum(pauli_weights)/len(pauli_weights)
    expect = 3.2580645161290325
    assert actual == pytest.approx(expect, abs=1.0e-8)


def test_sa_H4(fermion_hamiltonian_H4):

    num_qubits = fermion_hamiltonian_H4.one_body_tensor.shape[0]

    def objective_func(ttspec: TernaryTreeSpec) -> float:
        f2q_mapper = F2QMapper(ttspec=ttspec)
        qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H4)
        pauli_weights = qubit_hamiltonian.pauli_weights()
        return sum(pauli_weights) / len(pauli_weights)

    ttspec_opt = SA(
        num_qubits = num_qubits,
        objective_func = objective_func,
        params = SAParams(
            init_sampling = 10,
            num_steps = 3,
            cooling_factor = 1.0,
            seed = 123,
        ),
    ).optimize()
    
    f2q_mapper = F2QMapper(ttspec=ttspec_opt)
    qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H4)

    pauli_weights = qubit_hamiltonian.pauli_weights()
    actual = sum(pauli_weights)/len(pauli_weights)
    expect = 4.454054054054054
    assert actual == pytest.approx(expect, abs=1.0e-8)
