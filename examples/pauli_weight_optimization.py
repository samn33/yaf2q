""" pauli weight optimization using yaf2q """
from openfermion import MolecularData
from openfermionpyscf import run_pyscf

from yaf2q.ternary_tree_spec import TernaryTreeSpec
from yaf2q.f2q_mapper import F2QMapper
from yaf2q.optimizer.sa import SAParams, SA

    
def main():

    distance = 0.65
    molecule = MolecularData(
        geometry = [('H', (0, 0, 0)), ('H', (0, 0, distance)), ('H', (0, 0, 2*distance)), ('H', (0, 0, 3*distance))],
        basis = 'sto-3g',
        multiplicity = 1,
        charge = 0,
    )
    molecule = run_pyscf(molecule, run_scf=1, run_fci=1)
    fermion_hamiltonian = molecule.get_molecular_hamiltonian()
    num_qubits = fermion_hamiltonian.one_body_tensor.shape[0]

    for kind in ("jordan-wigner", "parity", "bravyi-kitaev"):
        print(f"[{kind}]")
        f2q_mapper = F2QMapper(num_qubits=num_qubits, kind=kind)
        qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian)
        weights = qubit_hamiltonian.pauli_weights()
        weight_ave = sum(weights) / len(weights)
        print(f"* pauli weight (ave): {weight_ave}")

    print("[ternary tree optimization]")
    def objective_func(ttspec: TernaryTreeSpec):
        f2q_mapper = F2QMapper(ttspec=ttspec)
        qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian)
        weights = qubit_hamiltonian.pauli_weights()
        weight_ave = sum(weights) / len(weights)
        return weight_ave
    
    ttspec_opt = SA(
            num_qubits = num_qubits,
            objective_func = objective_func,
            params = SAParams(init_sampling=10, num_steps=30, cooling_factor=1.0),
            verbose = False,
    ).optimize()

    print(f"* ternary tree:\n{ttspec_opt}")
    f2q_mapper = F2QMapper(ttspec=ttspec_opt)
    qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian)
    weights = qubit_hamiltonian.pauli_weights()
    weight_ave = sum(weights) / len(weights)

    print(f"* pauli weight (ave): {weight_ave}")

if __name__ == "__main__":
    main()
