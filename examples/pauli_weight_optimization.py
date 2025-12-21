""" pauli weight optimization using yaf2q """
import numpy as np
from openfermion import MolecularData
from openfermionpyscf import run_pyscf

from yaf2q.ternary_tree_spec import TernaryTreeSpec
from yaf2q.f2q_mapper import F2QMapper
from yaf2q.optimizer.sa import SAParams, SA

def molecule_H4():

    distance = 0.65
    molecule = MolecularData(
        geometry = [('H', (0, 0, 0)), ('H', (0, 0, distance)), ('H', (0, 0, 2*distance)), ('H', (0, 0, 3*distance))],
        basis = 'sto-3g',
        multiplicity = 1,
        charge = 0,
    )
    return run_pyscf(molecule, run_scf=1, run_fci=1)


def molecule_BeH2():

    distance = 1.33
    molecule = MolecularData(
        geometry = [('Be',(0.,0.,0.)), ('H',(0.,0.,-distance)), ('H',(0.,0.,distance))],
        basis = 'sto-3g',
        multiplicity=1,
        charge=0,
    )
    return run_pyscf(molecule, run_scf=1, run_fci=1)


def molecule_H2O():

    distance = 0.96
    angle_deg = 104.5
    angle_half_rad = np.deg2rad(angle_deg / 2)
    x = distance * np.sin(angle_half_rad)
    z = distance * np.cos(angle_half_rad)
    molecule = MolecularData(
        geometry = [('O',(0., 0., 0.)), ('H',(x, 0., z)), ('H',(-x, 0., z))],
        basis = 'sto-3g',
        multiplicity = 1,
        charge = 0,
    )
    return run_pyscf(molecule, run_scf=1, run_fci=1)


def molecule_CH4():

    geometry = [
        ('C', (0.0, 0.0, 0.0)),
        ('H', (0.6291, 0.6291, 0.6291)),
        ('H', (-0.6291, -0.6291, 0.6291)),
        ('H', (0.6291, -0.6291, -0.6291)),
        ('H', (-0.6291, 0.6291, -0.6291))
    ]
    molecule = MolecularData(
        geometry=geometry,
        basis = 'sto-3g',
        multiplicity = 1,
        charge = 0,
    )
    return run_pyscf(molecule, run_scf=1, run_fci=1)


def main():

    #molecule = molecule_H4()
    #molecule = molecule_BeH2()
    molecule = molecule_CH4()
    #molecule = molecule_H2O()
    fermion_hamiltonian = molecule.get_molecular_hamiltonian()
    num_qubits = fermion_hamiltonian.one_body_tensor.shape[0]
    print(f"num_qubits = {num_qubits}")

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
            verbose = True,
    ).optimize()

    print(f"* ternary tree:\n{ttspec_opt}")
    f2q_mapper = F2QMapper(ttspec=ttspec_opt)
    qubit_hamiltonian = f2q_mapper.fermion_to_qubit_operator(fermion_hamiltonian)
    weights = qubit_hamiltonian.pauli_weights()
    weight_ave = sum(weights) / len(weights)

    print(f"* pauli weight (ave): {weight_ave}")


if __name__ == "__main__":
    main()
