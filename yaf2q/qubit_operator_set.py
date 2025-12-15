import numpy as np
from dataclasses import dataclass, field
from scipy.sparse.linalg import eigsh
from openfermion.linalg import qubit_operator_sparse
from openfermion.ops.operators.qubit_operator import QubitOperator
from qiskit.quantum_info import SparsePauliOp


@dataclass
class QubitOperatorSet:
    """
    Set of qubit operators

    Attributes
    ----------
    num_qubits : int
        number of qubits
    openfermion_form : QubitOperator
        qubit operator in openfermion format
    qiskit_form :
        qubit operator in qiskit format

    """
    num_qubits: int                 = field(default=0, init=True)
    openfermion_form: QubitOperator = field(default=None, init=True)
    _qiskit_form: SparsePauliOp     = field(default=None, init=False)

    @property
    def qiskit_form(self):
        """ getter of the qiskit_form """

        if self._qiskit_form is not None:
            return self._qiskit_form
        
        # qubit operator in qiskit format
        qo_list = []
        for pp, coef in self.openfermion_form.terms.items():
            pp_dict = {}
            for i,p in pp:
                pp_dict[i] = p
        
            pp_list = []
            for i in range(self.num_qubits):
                if i in pp_dict:
                    pp_list.append((i, pp_dict[i]))
                else:
                    pp_list.append((i, 'I'))
        
            pp_str = ""
            for i, p in pp_list:
                pp_str += p
        
            qo_list.append((pp_str, coef))
        
        qo_list_sorted = sorted(qo_list)
        
        self._qiskit_form = SparsePauliOp.from_list(qo_list_sorted)

        return self._qiskit_form


    def __str__(self) -> str:
        return self.to_string()


    def to_string(self) -> str:
        """
        Get the string of QubitOperator
        
        Parameters
        ----------
        None
        
        Returns
        -------
        str
            string of the QubitOperator

        """
        pp_list = []
        for pp, coef in self.openfermion_form.terms.items():
            pp_str = ""
            if len(pp) == 0:
                pp_str += "I"
            for i,p in pp:
                pp_str += f"{p}[{i}]"
            pp_str += f" {coef}"
            pp_list.append(pp_str)
        pp_list = sorted(pp_list)

        s = ""
        for pp_str in pp_list:
            s += (pp_str + "\n")

        return s[:-1]


    def eigsh(self, num:int = 1) -> tuple[np.ndarray,np.ndarray[np.ndarray]]:
        """
        Get the eigen-values and eigen-vectors of QubitOperator
        
        Parameters
        ----------
        num : int
            number of the eigen-values,eigen-vectors in order from smallest to largest
        
        Returns
        -------
        eigenvalues, eigenvectors : tuple[list,list]
            list of the eigenvalues

        """
        operator_matrix = qubit_operator_sparse(self.openfermion_form, n_qubits=self.num_qubits)
        eigenvalues, eigenvectors = eigsh(operator_matrix, k=num, which="SA")

        return eigenvalues, eigenvectors


    def eigenvalues(self, num:int = 1) -> np.ndarray[float]:
        """
        Get the eigen-values of QubitOperator
        
        Parameters
        ----------
        num : int
            number of the eigen-values in order from smallest to largest
        
        Returns
        -------
        eigenvalues : list[float]
            list of eigenvalues

        """
        eigenvalues, _ = self.eigsh(num=num)

        return eigenvalues


    def eigenvectors(self, num:int = 1) -> np.ndarray[np.ndarray[complex]]:
        """
        Get the eigen-vectors of QubitOperator
        
        Parameters
        ----------
        num : int
            number of the eigen-values in order from smallest to largest
        
        Returns
        -------
        eigenvectors : np.ndarray
            list of eigenvectors

        """
        _, eigenvectors = self.eigsh(num=num)

        return eigenvectors


    def pauli_weights(self) -> list[int]:
        """
        Get the pauli weights of QubitOperator
        
        Parameters
        ----------
        None
        
        Returns
        -------
        weights : list[int]
            pauli weights of QubitOperator

        Notes
        -----
        If the qubit operator is 0.1*X[0]*Z[1]+0.2*X[1]*Y[2]*Z[3], then returned list is [2, 3].

        """
        weights = []
        for pp, coef in self.openfermion_form.terms.items():
            weights.append(len(pp))

        return weights
