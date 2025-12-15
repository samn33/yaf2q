import pytest
import numpy as np
from openfermion import MolecularData
from openfermionpyscf import run_pyscf


@pytest.fixture
def fermion_hamiltonian_H2():

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

    return fermion_hamiltonian
    

@pytest.fixture
def fermion_hamiltonian_H3():

    basis = 'sto-3g'
    distance = 0.65
    geometry = [('H', (0, 0, 0)), ('H', (0, 0, distance)), ('H', (0, 0, 2*distance))]
    multiplicity = 2
    charge = 0

    molecule = MolecularData(
        geometry=geometry,
        basis=basis,
        multiplicity=multiplicity,
        charge=charge,
    )
    molecule = run_pyscf(molecule, run_scf=1, run_fci=1)
    fermion_hamiltonian = molecule.get_molecular_hamiltonian()

    return fermion_hamiltonian
    

@pytest.fixture
def fermion_hamiltonian_H4():

    basis = 'sto-3g'
    distance = 0.65
    geometry = [('H', (0, 0, 0)), ('H', (0, 0, distance)), ('H', (0, 0, 2*distance)), ('H', (0, 0, 3*distance))]
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

    return fermion_hamiltonian


@pytest.fixture
def fermion_hamiltonian_LiH():

    basis = 'sto-3g'
    distance = 1.546
    geometry = [('Li',(0.,0.,0.)),('H',(0.,0.,distance))]
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

    return fermion_hamiltonian


@pytest.fixture
def fermion_hamiltonian_BeH2():

    basis = 'sto-3g'
    distance = 1.33
    geometry = [
        ('Be',(0.,0.,0.)),
        ('H',(0.,0.,-distance)),
        ('H',(0.,0.,distance)),
    ]
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

    return fermion_hamiltonian


@pytest.fixture
def fermion_hamiltonian_H2O():

    basis = 'sto-3g'
    distance = 0.96
    angle_deg = 104.5
    angle_half_rad = np.deg2rad(angle_deg / 2)
    x = distance * np.sin(angle_half_rad)
    z = distance * np.cos(angle_half_rad)
    geometry = [
        ('O',(0., 0., 0.)),
        ('H',(x, 0., z)),
        ('H',(-x, 0., z)),
    ]
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

    return fermion_hamiltonian
