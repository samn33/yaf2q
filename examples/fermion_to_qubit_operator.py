""" Simple fermion to qubit mapper using yaf2q """
import argparse
import random
from openfermion import MolecularData
from openfermionpyscf import run_pyscf

from yaf2q.ternary_tree_spec import TernaryTreeSpec
from yaf2q.f2q_mapper import F2QMapper


def main():

    parser = argparse.ArgumentParser(description="Simple fermion to qubit mapper")
    parser.add_argument("--molecule", "-m", help="molecule name (H2,H3,H4)",
                        default="H2", type=str)
    parser.add_argument("--kind", "-k", help="kind of fermion to qubit mapping (jordan-wigner,parity,bravyi-kitaev,random-tt)",
                        default="random-tt", type=str)
    args = parser.parse_args()

    if args.kind not in ("jordan-wigner", "bravyi-kitaev", "parity", "random-tt"):
        raise ValueError("kind must be jordan-wigner, bravyi-kitaev, or random-tt")

    if args.molecule == "H2":
        basis = 'sto-3g'
        distance = 0.65
        geometry = [('H', (0, 0, 0)), ('H', (0, 0, distance))]
        multiplicity = 1
        charge = 0
    elif args.molecule == "H3":
        basis = 'sto-3g'
        distance = 0.65
        geometry = [('H', (0, 0, 0)), ('H', (0, 0, distance)), ('H', (0, 0, 2*distance))]
        multiplicity = 2
        charge = 0
    elif args.molecule == "H4":
        basis = 'sto-3g'
        distance = 0.65
        geometry = [('H', (0, 0, 0)), ('H', (0, 0, distance)), ('H', (0, 0, 2*distance)), ('H', (0, 0, 3*distance))]
        multiplicity = 1
        charge = 0
    else:
        raise ValueError("molecule must be H2, H3, or H4.")

    # fermion hamiltonian
    molecule = MolecularData(
        geometry=geometry,
        basis=basis,
        multiplicity=multiplicity,
        charge=charge,
    )
    molecule = run_pyscf(molecule, run_scf=1, run_fci=1)
    fermion_hamiltonian = molecule.get_molecular_hamiltonian()
    num_qubits = fermion_hamiltonian.one_body_tensor.shape[0]

    # fermion to qubit mapper
    if args.kind in ("jordan-wigner", "bravyi-kitaev", "parity"):
        f2q_mapper = F2QMapper(kind=args.kind, num_qubits=num_qubits)
        print(f"* f2q mapper:\n{f2q_mapper}")
    elif args.kind == "random-tt":
        ttspec = TernaryTreeSpec.random(num_qubits)
        f2q_mapper = F2QMapper(ttspec=ttspec)
        print(f"* f2q mapper:\n{f2q_mapper}")
    else:
        raise ValueError("unknown kind is specified.")

    # qubit hamiltonian
    qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian)

    # fock state to qubit state
    fock_state = [1 if x < num_qubits / 2 else 0 for x in range(num_qubits)]
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    print("* fock to qubit state:")
    print(f"{fock_state} -> {qubit_state}")
    for _ in range(3):
        fock_state_tmp = random.sample(fock_state, len(fock_state))
        qubit_state_tmp = f2q_mapper.fock_to_qubit_state(fock_state_tmp)
        print(f"{fock_state_tmp} -> {qubit_state_tmp}")
    
    # average length of pauli product
    weights = qubit_hamiltonian.pauli_weights()
    print(f"* average length of pauli product:\n{sum(weights) / len(weights)}")

    # ground state energy
    print(f"* ground state energy:\n{qubit_hamiltonian.eigenvalues()[0]}")


if __name__ == "__main__":
    main()
