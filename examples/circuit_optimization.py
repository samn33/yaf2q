""" QPE circut optimizer using yaf2q """
import numpy as np
import argparse
import warnings

from openfermion import MolecularData
from openfermionpyscf import run_pyscf

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import PhaseEstimation
from qiskit.circuit.library import PauliEvolutionGate
from qiskit_aer import AerSimulator
from qiskit.synthesis import SuzukiTrotter

from yaf2q.ternary_tree_spec import TernaryTreeSpec
from yaf2q.f2q_mapper import F2QMapper
from yaf2q.optimizer.sa import SAParams, SA

warnings.simplefilter("ignore", category=DeprecationWarning) # for PhaseEstimation


NUM_ANCILLA_QUBITS = 8
TIME_INTERVAL = 0.5


def make_qcircuit(qubit_state, qubit_hamiltonian, num_ancilla_qubits, qiskit_optlevel: int = 0):
    """ make QPE circuit """

    # Unitary operator: U = exp(-i * H * t)
    num_qubits = qubit_hamiltonian.num_qubits
    evolution_gate = PauliEvolutionGate(qubit_hamiltonian.qiskit_form, time=TIME_INTERVAL, synthesis=SuzukiTrotter())
    unitary_circuit = QuantumCircuit(num_qubits)
    unitary_circuit.append(evolution_gate, range(num_qubits))

    # QPE circuit
    qpe_logic = PhaseEstimation(num_ancilla_qubits, unitary_circuit)
    init_circuit = QuantumCircuit(num_ancilla_qubits + num_qubits)
    for i, q in enumerate(qubit_state): # state preparation
        if q == 1:
            init_circuit.x(num_ancilla_qubits + i)
    init_circuit.compose(qpe_logic, inplace=True)
    init_circuit.measure_all()

    # Transpile
    backend = AerSimulator()
    qc = transpile(init_circuit, backend)

    return qc


def get_energy(qc: QuantumCircuit, num_ancilla_qubits: int) -> float:

    backend = AerSimulator()
    counts = backend.run(qc, shots=2048).result().get_counts()
    highest_bitstring = max(counts, key=counts.get).split()[0][::-1][:num_ancilla_qubits]
    phase = int(highest_bitstring, 2) / (2**num_ancilla_qubits)

    # e^{i 2 pi phi} = e^{-i E t} => E = - (2π * phi) / t 
    if phase > TIME_INTERVAL:
        phase -= 1.0
    energy = -(phase * 2 * np.pi) / TIME_INTERVAL
    return energy

    
def main():

    parser = argparse.ArgumentParser(description="Experimental program for optimizing QPE circuit using yaf2q")
    parser.add_argument("--kind", "-k", help="Kind of Fermion to Qubit Operator mapping(jordan-wigner/parity/bravyi-kitaev/optimize-sa)",
                        default="jordan-wigner", type=str)
    parser.add_argument("--qiskit-optlevel", "-o", help="Qiskit circuit optimization level (0:none,1,2,3:best)",
                        default=0, type=int)
    parser.add_argument("--num-ancilla", "-a", help="Number of ancilla qubits",
                        default=NUM_ANCILLA_QUBITS, type=int)
    parser.add_argument("--calc-energy", "-e", help="Specify when calculating ground state energy",
                        action="store_true")
    args = parser.parse_args()
    kind = args.kind
    qiskit_optlevel = args.qiskit_optlevel
    num_ancilla = args.num_ancilla
    calc_energy = args.calc_energy

    basis = 'sto-3g'
    distance = 0.65
    geometry = [('H', (0, 0, 0)), ('H', (0, 0, distance))]
    multiplicity = 1
    charge = 0
        
    molecule = MolecularData(
        geometry=geometry,
        basis=basis,
        multiplicity=multiplicity,
        charge=charge,
    )
    molecule = run_pyscf(molecule, run_scf=1, run_fci=1)
    fermion_hamiltonian = molecule.get_molecular_hamiltonian()

    num_qubits = fermion_hamiltonian.one_body_tensor.shape[0]
    
    if qiskit_optlevel < 0 or qiskit_optlevel > 3:
        raise ValueError("qiskit_optlevel must be 0,1,2 or 3.")

    if kind not in ("jordan-wigner", "bravyi-kitaev", "parity", "optimize-sa"):
        raise ValueError("kind must be jordan-wigner, bravyi-kitaev, or optimize-sa.")

    def objective_func(ttspec: TernaryTreeSpec):
        num_qubits = ttspec.num_qubits
        f2q_mapper = F2QMapper(ttspec=ttspec)
        qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian)
        fock_state = [1 if x < num_qubits / 2 else 0 for x in range(num_qubits)]
        qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
        qc = make_qcircuit(
            qubit_state = qubit_state,
            qubit_hamiltonian = qubit_hamiltonian,
            num_ancilla_qubits = num_ancilla,
            qiskit_optlevel = qiskit_optlevel,
        )
        obj = qc.depth() # depth
        #obj = qc.count_ops().get('cx', 0) # cnot-count
        #obj = sum(list(qc.count_ops().values())) # gate-count
        return obj
    
    # get the quantum circuit
    if kind in ("jordan-wigner","parity","bravyi-kitaev"):
        f2q_mapper = F2QMapper(num_qubits=num_qubits, kind=kind)
        qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian)
        fock_state = [1 if x < num_qubits / 2 else 0 for x in range(num_qubits)]
        qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
        qc = make_qcircuit(
            qubit_state = qubit_state,
            qubit_hamiltonian = qubit_hamiltonian,
            num_ancilla_qubits = num_ancilla,
            qiskit_optlevel = qiskit_optlevel,
        )
    elif kind == "optimize-sa":
        ttspec_opt = SA(
            num_qubits = num_qubits,
            objective_func = objective_func,
            params = SAParams(
                init_sampling = 10,
                num_steps = 5,
                cooling_factor = 1.0,
            ),
            verbose = True,
        ).optimize()
        f2q_mapper = F2QMapper(ttspec=ttspec_opt)
        qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian)
        fock_state = [1 if x < num_qubits / 2 else 0 for x in range(num_qubits)]
        qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
        qc = make_qcircuit(
            qubit_state = qubit_state,
            qubit_hamiltonian = qubit_hamiltonian,
            num_ancilla_qubits = num_ancilla,
            qiskit_optlevel = qiskit_optlevel,
        )
    else:
        raise ValueError("kind must be jordan-wigner,bravyi-kitaev,or optimize-sa.")
    
    # result
    print(f"kind = {kind}")
    weights = qubit_hamiltonian.pauli_weights()
    print("[pauli weight]")
    print(f"ave = {sum(weights) / len(weights)}")
    print(f"min = {min(weights)}")
    print(f"max = {max(weights)}")

    depth = qc.depth()
    count_ops = qc.count_ops()
    cnot_count = count_ops.get('cx', 0)
    gate_count = sum(list(count_ops.values()))
    print("[circuit stats]")
    print(f"depth      = {depth}")
    print(f"cnot count = {cnot_count}")
    print(f"gate count = {gate_count}")

    if calc_energy is True:
        print("[ground state energy]")
        print(f"diagonalization(correct) = {qubit_hamiltonian.eigenvalues()[0]}")
        print(f"quantum phase estimation = {get_energy(qc, num_ancilla)}")


if __name__ == "__main__":
    main()
