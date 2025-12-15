""" Usage example of ternary tree mapper """
import random

from yaf2q.ternary_tree_spec import TernaryTreeSpec
from yaf2q.f2q_mapper import F2QMapper


def main():

    # ternary tree mapper
    ttspec = TernaryTreeSpec(
        indices = [1, 0, 3, 2],
        edges = {1: (0, 'X'), 2: (1, 'Y'), 3: (1, 'Z')},
    )
    ttspec.show()
    f2q_mapper = F2QMapper(ttspec=ttspec)
    print(f"* f2q mapper:\n{f2q_mapper}")

    # endoding matrix
    encoding_matrix = f2q_mapper.encoding_matrix
    print(f"* encoding_matrix:\n{encoding_matrix}")

    # fock state to qubit state
    num_qubits = f2q_mapper.num_qubits
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
