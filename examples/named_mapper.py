""" Usage example of named mapper """
import random

from yaf2q.f2q_mapper import F2QMapper


def main():

    #random.seed(123)
    num_qubits = 4
    #kind = "jordan-wigner"
    #kind = "parity"
    kind = "bravyi-kitaev"

    # conventional named mapper
    f2q_mapper = F2QMapper(kind=kind, num_qubits=num_qubits)
    print(f"* f2q mapper:\n{f2q_mapper}")

    # encoding matrix
    encoding_matrix = f2q_mapper.encoding_matrix
    print(f"* encoding_matrix:\n{encoding_matrix}")

    # fock state to qubit state
    fock_state = [1 if x < num_qubits / 2 else 0 for x in range(num_qubits)]
    qubit_state = f2q_mapper.fock_to_qubit_state(fock_state)
    print("* fock to qubit state (some examples):")
    print(f"{fock_state} -> {qubit_state}")
    for _ in range(3):
        fock_state_tmp = random.sample(fock_state, len(fock_state))
        qubit_state_tmp = f2q_mapper.fock_to_qubit_state(fock_state_tmp)
        print(f"{fock_state_tmp} -> {qubit_state_tmp}")


if __name__ == "__main__":
    main()
