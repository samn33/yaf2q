import pytest

from yaf2q.ternary_tree_spec import TernaryTreeSpec
from yaf2q.f2q_mapper import F2QMapper


def test_qiskit_form_jordan_wigner(fermion_hamiltonian_H2):

    num_qubits = fermion_hamiltonian_H2.one_body_tensor.shape[0]
    f2q_mapper = F2QMapper(num_qubits=num_qubits, kind="jordan-wigner")
    qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2)
    qiskit_hamiltonian = qubit_hamiltonian.qiskit_form
    
    assert (qiskit_hamiltonian.paulis.to_labels()
            == ['IIII', 'IIIZ', 'IIZI', 'IIZZ', 'IZII', 'IZIZ', 'IZZI', 'XXYY',
                'XYYX', 'YXXY', 'YYXX', 'ZIII', 'ZIIZ', 'ZIZI', 'ZZII'])
    assert (qiskit_hamiltonian.coeffs.tolist()
            == pytest.approx([(0.03775110394645531+0j), (0.18601648886230604+0j),
                              (0.18601648886230604+0j), (0.17297610130745106+0j),
                              (-0.26941693141631995+0j), (0.12584136558006329+0j),
                              (0.1699209784826151+0j), (-0.04407961290255181+0j),
                              (0.04407961290255181+0j), (0.04407961290255181+0j),
                              (-0.04407961290255181+0j), (-0.26941693141631995+0j),
                              (0.1699209784826151+0j), (0.12584136558006329+0j),
                              (0.17866777775953394+0j)]))


def test_qiskit_form_parity(fermion_hamiltonian_H2):

    num_qubits = fermion_hamiltonian_H2.one_body_tensor.shape[0]
    f2q_mapper = F2QMapper(num_qubits=num_qubits, kind="parity")
    qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2)
    qiskit_hamiltonian = qubit_hamiltonian.qiskit_form
    
    assert (qiskit_hamiltonian.paulis.to_labels()
            == ['IIII', 'IIIZ', 'IIZI', 'IIZZ', 'IXZX', 'IYIY', 'IZIZ', 'IZZI',
                'IZZZ', 'ZIZI', 'ZXZX', 'ZYIY', 'ZZII', 'ZZIZ', 'ZZZZ'])
    assert (qiskit_hamiltonian.coeffs.tolist()
            == pytest.approx([(0.03775110394645531+0j), (0.18601648886230604+0j),
                              (0.17297610130745106+0j), (0.18601648886230604+0j),
                              (0.04407961290255181+0j), (0.04407961290255181+0j),
                              (0.1699209784826151+0j), (-0.26941693141631995+0j),
                              (0.12584136558006329+0j), (0.17866777775953394+0j),
                              (0.04407961290255181+0j), (0.04407961290255181+0j),
                              (-0.26941693141631995+0j), (0.1699209784826151+0j),
                              (0.12584136558006329+0j)]))


def test_qiskit_form_bravyi_kitaev(fermion_hamiltonian_H2):

    num_qubits = fermion_hamiltonian_H2.one_body_tensor.shape[0]
    f2q_mapper = F2QMapper(num_qubits=num_qubits, kind="bravyi-kitaev")
    qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2)
    qiskit_hamiltonian = qubit_hamiltonian.qiskit_form
    
    assert (qiskit_hamiltonian.paulis.to_labels()
            == ['IIII', 'IIIZ', 'IIZI', 'IIZZ', 'IXZX', 'IYZY', 'IZII', 'IZIZ',
                'IZZZ', 'ZIZI', 'ZXZX', 'ZYZY', 'ZZIZ', 'ZZZI', 'ZZZZ'])
    assert (qiskit_hamiltonian.coeffs.tolist()
            == pytest.approx([(0.03775110394645509+0j), (0.18601648886230604+0j),
                              (0.17297610130745106+0j), (0.18601648886230604+0j),
                              (0.04407961290255181+0j), (0.04407961290255181+0j),
                              (-0.26941693141631995+0j), (0.12584136558006329+0j),
                              (0.1699209784826151+0j), (0.17866777775953394+0j),
                              (0.04407961290255181+0j), (0.04407961290255181+0j),
                              (0.12584136558006329+0j), (-0.26941693141631995+0j),
                              (0.1699209784826151+0j)]))


def test_qiskit_form_ternary_tree(fermion_hamiltonian_H2):

    f2q_mapper = F2QMapper(
        ttspec = TernaryTreeSpec(
            indices = [1, 0, 3, 2],
            edges = {1: (0, 'X'), 2: (1, 'X'), 3: (0, 'Y')},
        )
    )
    qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian_H2)
    qiskit_hamiltonian = qubit_hamiltonian.qiskit_form
    
    assert (qiskit_hamiltonian.paulis.to_labels()
            == ['IIII', 'IIIZ', 'IIZI', 'IIZZ', 'IZZI', 'XIZX', 'XZZX', 'YIIY',
                'YZIY', 'ZIII', 'ZIIZ', 'ZIZZ', 'ZZIZ', 'ZZZI', 'ZZZZ'])
    assert (qiskit_hamiltonian.coeffs.tolist()
            == pytest.approx([(0.03775110394645542+0j), (0.18601648886230604+0j),
                              (0.17297610130745106+0j), (0.18601648886230604+0j),
                              (0.17866777775953394+0j), (0.04407961290255181+0j),
                              (0.04407961290255181+0j), (0.04407961290255181+0j),
                              (0.04407961290255181+0j), (-0.26941693141631995+0j),
                              (0.1699209784826151+0j), (0.12584136558006329+0j),
                              (0.1699209784826151+0j), (-0.26941693141631995+0j),
                              (0.12584136558006329+0j)]))
