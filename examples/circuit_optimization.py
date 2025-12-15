""" QPE circut optimizer using yaf2q """
import numpy as np
import argparse

from openfermion import MolecularData
from openfermionpyscf import run_pyscf

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import QFTGate
from qiskit_aer import AerSimulator
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.synthesis import SuzukiTrotter

from yaf2q.ternary_tree_spec import TernaryTreeSpec
from yaf2q.f2q_mapper import F2QMapper
from yaf2q.optimizer.sa import SAParams, SA

REPS = 10
NUM_ANCILLA_QUBITS = 4

def make_qcircuit(qubit_state, qubit_hamiltonian,
                  num_ancilla_qubits, qiskit_optlevel: int = 0):
    """ make QPE circuit """

    num_qubits = qubit_hamiltonian.num_qubits
    qiskit_hamiltonian = qubit_hamiltonian.qiskit_form
    estimated_energy = abs(qubit_hamiltonian.eigenvalues(num=1)[0])
    evolution_time = np.pi / estimated_energy

    qc = QuantumCircuit(num_qubits + num_ancilla_qubits, num_ancilla_qubits)

    # hadamard gates for ancilla qubits
    for i in range(num_ancilla_qubits):
        qc.h(i)

    # state preparation
    for i, q in enumerate(qubit_state):
        if q == 1:
            qc.x(num_ancilla_qubits + i)
    
    # controlled time evolution operator
    trotter_gate = SuzukiTrotter(reps=REPS)
    evolution_op = PauliEvolutionGate(qiskit_hamiltonian, time=evolution_time, synthesis=trotter_gate)
    controlled_evolution_op = evolution_op.control()
    for i in range(num_ancilla_qubits):
        t = 2 ** i
        for _ in range(t):
            qc.append(
                controlled_evolution_op,
                [i] + list(range(num_ancilla_qubits, num_qubits + num_ancilla_qubits)),
            )

    # inverse quantum fourier transform
    qft_gate = QFTGate(num_ancilla_qubits)
    iqft_gate = qft_gate.inverse()
    qc.append(iqft_gate, range(num_ancilla_qubits))

    # measurement
    qc.measure(range(num_ancilla_qubits), range(num_ancilla_qubits))

    # decompose circuit
    qc = qc.decompose()

    # backend
    simulator = AerSimulator()

    # circuit optimization by qiskit
    qc = transpile(qc, simulator, optimization_level=qiskit_optlevel)

    return qc


def get_energy(qc: QuantumCircuit, num_ancilla_qubits: int, estimated_energy: float) -> float:

    estimated_energy = estimated_energy
    evolution_time = np.pi / estimated_energy

    # run by qiskit-aer simulator
    simulator = AerSimulator()
    job = simulator.run(qc, shots=1000)
    result = job.result()

    # evaluate the energy
    counts = result.get_counts(qc)
    most_common_bitstring = max(counts, key=counts.get)
    phase = int(most_common_bitstring, 2) / (2 ** num_ancilla_qubits)
    energy = - 2.0 * np.pi * phase / evolution_time

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
    geometry = [('H', (0, 0, 0)),
                ('H', (0, 0, distance))]
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
        depth = qc.depth()
        #count_ops = qc.count_ops()
        #cnot_count = count_ops.get('cx', 0)
        #gate_count = sum(list(count_ops.values()))

        return depth
    
    # get the quantum circuit
    if kind in ("jordan-wigner","bravyi-kitaev","parity"):
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
                num_steps = 10,
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
        estimated_energy = abs(qubit_hamiltonian.eigenvalues(num=1)[0])
        print("[ground state energy]")
        print(f"diagonalization(correct) = {qubit_hamiltonian.eigenvalues()[0]}")
        print(f"quantum phase estimation = {get_energy(qc, num_ancilla, estimated_energy)}")


if __name__ == "__main__":
    main()
