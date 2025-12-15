import copy
import numpy as np

import openfermion as of
from openfermion.ops import InteractionOperator
from openfermion.transforms import jordan_wigner, bravyi_kitaev

from yaf2q.ternary_tree_spec import TernaryTreeSpec
from yaf2q.yaf2q import yaf2q_operator, yaf2q_encoding_matrix, yaf2q_encoding_matrix_inv
from yaf2q.qubit_operator_set import QubitOperatorSet


def _is_power_of_two(n: int) -> bool:
    return n > 0 and (n & (n - 1) == 0)


class F2QMapper:
    """
    Fermion to qubit mapper

    Attributes
    ----------
    kind : str
        kind of the fermion to qubit mapper (jordan-wigner,parity,bravyi-kitaev)
    num_qubits : int
        number of qubits
    ttspec : TernaryTreeSpec
        ternary tree specification
    encoding_matrix : numpy.ndarray
        endoding matrix for the fermion to qubit mapper
    encoding_matrix_inv : numpy.ndarray
        inverse of the endoding matrix for the fermion to qubit mapper

    """
    def __init__(self, kind: str | None = None,
                 num_qubits: int = 0,
                 ttspec: TernaryTreeSpec | None = None):
        """
        Constructor of F2QMapper

        Parameters
        ----------
        kind : str
            kind of fermion to qubit mapper (jordan-wigner,parity,bravyi-kitaev)
        num_qubits : int
            number of qubits
        ttspec : TernaryTreeSpec
            ternary tree specification

        """
        if kind is None and ttspec is None:
            raise ValueError("kind and ttspec are both not specified. need to specify either.")
        elif kind is not None and ttspec is not None:
            raise ValueError("kind and ttspec are both specified. need to specify either.")
        elif kind is not None and num_qubits <= 0:
            raise ValueError("num_qubits is not specified.")
        elif ttspec is not None and num_qubits > 0:
            raise ValueError("num_qubits must not be specified when ttspec is specified.")

        if (kind is not None and (not isinstance(kind, str) or
                                  kind not in ("jordan-wigner","parity","bravyi-kitaev"))):
            raise TypeError("kind must be a string jordan-wigner, parity, or bravyi-kitaev.")
        elif ttspec is not None and not isinstance(ttspec, TernaryTreeSpec):
            raise TypeError("ttspec must be a TernaryTreeSpec.")

        if kind in ("jordan-wigner", "parity"):
            self._num_qubits = num_qubits
            self._encoding_matrix = np.array(yaf2q_encoding_matrix(kind, None, num_qubits))
            self._encoding_matrix_inv = np.array(yaf2q_encoding_matrix_inv(kind, None, num_qubits))
        elif kind == "bravyi-kitaev":
            self._num_qubits = num_qubits
            self._encoding_matrix = None
            if _is_power_of_two(num_qubits):
                self._encoding_matrix = np.array(yaf2q_encoding_matrix(kind, None, num_qubits))
                self._encoding_matrix_inv = np.array(yaf2q_encoding_matrix_inv(kind, None, num_qubits))
            else:
                raise ValueError("if kind is brabyi-kitaev, num_qubit must be a power of 2.")
        else:
            self._num_qubits = len(ttspec.indices)
            self._encoding_matrix = np.array(yaf2q_encoding_matrix(None, ttspec._to_string_rust(), num_qubits))
            self._encoding_matrix_inv = np.array(yaf2q_encoding_matrix_inv(None, ttspec._to_string_rust(), num_qubits))
        
        self._kind = kind
        self._ttspec = ttspec


    @property
    def kind(self) -> str:
        """ getter of the kind """
        return self._kind

    
    @property
    def num_qubits(self) -> int:
        """ getter of the num_qubits """
        return self._num_qubits

    
    @property
    def ttspec(self) -> int:
        """ getter of the ttspec """
        return self._ttspec

    
    @property
    def encoding_matrix(self) -> np.ndarray:
        """ getter of the encoding_matrix """
        return self._encoding_matrix

     
    @property
    def encoding_matrix_inv(self) -> np.ndarray:
        """ getter of the encoding_matrix_inv """
        return self._encoding_matrix_inv

    
    def __str__(self) -> str:
        return self.to_string()


    def to_string(self) -> str:
        """
        Get the string of the F2QMapper

        Parameters
        ----------
        None
        
        Returns
        -------
        str
            string of the F2QMapper

        """
        s = ""
        if self._kind is not None:
            s += f"named mapper - {self._kind} (num_qubits:{self._num_qubits})"
        else:
            s += f"ternary tree mapper - {self._ttspec}"
        return s


    def __eq__(self, other) -> bool:
        return self._kind == other.kind and self._num_qubits == other.num_qubits and self._ttspec == other.ttspec

        
    def _fermion_to_qubit_operator_tt(self, fermion_operator: InteractionOperator) -> QubitOperatorSet:
        """
        Get the qubit operator from the fermion operator using ternary tree
        
        Parameters
        ----------
        fermion_operator : InteractionOperator
            fermion operator

        Returns
        -------
        QubitOperatorSet
            set of the qubit operators

        """
        num_qubits = 0
        if hasattr(fermion_operator, "one_body_tensor"):
            num_qubits = fermion_operator.one_body_tensor.shape[0]
        elif hasattr(fermion_operator, "two_body_tensor"):
            num_qubits = fermion_operator.two_body_tensor.shape[0]

        if self._num_qubits != num_qubits:
            print(self._num_qubits, num_qubits, type(fermion_operator))
            raise ValueError("the number of qubits of F2QMapper and fermion operator do not match.")

        # qubit operator in openfermion format (openfermion_form)
        if self._kind is None:
            qo_dict = yaf2q_operator(None, self._ttspec._to_string_rust(), str(fermion_operator))
        elif self._kind == "jordan-wigner":
            qo_dict = yaf2q_operator(self._kind, None, str(fermion_operator))
        elif self._kind == "parity":
            qo_dict = yaf2q_operator(self._kind, None, str(fermion_operator))
        elif self._kind == "bravyi-kitaev":
            if _is_power_of_two(self._num_qubits):
                qo_dict = yaf2q_operator(self._kind, None, str(fermion_operator))
            else:
                raise ValueError("num_qubits must be power of two.")
        else:
            raise ValueError("invalid f2q_mapper is specified.")
        
        openfermion_form = of.ops.QubitOperator()
        for k,v in qo_dict.items():
            openfermion_form += of.ops.QubitOperator(k, v.real)
        
        return QubitOperatorSet(
            num_qubits = self._num_qubits,
            openfermion_form = openfermion_form,
        )


    def _fermion_to_qubit_operator_of(self, fermion_operator: InteractionOperator) -> QubitOperatorSet:
        """
        Get the qubit operator from the fermion operator by openfermion
        
        Parameters
        ----------
        fermion_operator : InteractionOperator
            fermion operator

        Returns
        -------
        QubitOperatorSet
            set of qubit operators

        """
        if self._num_qubits != fermion_operator.one_body_tensor.shape[0]:
            raise ValueError("the number of qubits of F2QMapper and fermion operator do not match.")
        
        # qubit operator in openfermion format (openfermion_form)
        if self._kind is None:
            raise ValueError("kind is not specified.")
        elif self._kind == "jordan-wigner":
            openfermion_form = jordan_wigner(fermion_operator)
            #encoding_matrix = None
        elif self._kind == "parity":
            raise ValueError("parity kind is not supported.")
        elif self._kind == "bravyi-kitaev":
            openfermion_form = bravyi_kitaev(fermion_operator)
            #encoding_matrix = None
        else:
            raise ValueError("invalid f2q_mapper is specified.")
        
        return QubitOperatorSet(
            num_qubits = self._num_qubits,
            openfermion_form = openfermion_form,
        )


    def fermion_to_qubit_operator(self, fermion_operator: InteractionOperator, method: str = "tt") -> QubitOperatorSet:
        """
        Get the qubit operator from the fermion operator
        
        Parameters
        ----------
        fermion_operator : InteractionOperator
            fermion operator

        Returns
        -------
        QubitOperatorSet
            set of the qubit operators

        """
        if method == "tt":
            return self._fermion_to_qubit_operator_tt(fermion_operator)
        elif method == "of":
            return self._fermion_to_qubit_operator_of(fermion_operator)
        else:
            raise ValueError("unknown method string is specified.")

        
    def fock_to_qubit_state(self, fock_state: list[int]) -> list[int]:
        """
        Get the qubit state from the fock state
        
        Parameters
        ----------
        fock_state : list[int]
            fock state

        Returns
        -------
        qubit_state : list[int]
            qubit state

        """
        if len(fock_state) != self._num_qubits:
            raise ValueError("dimension of the fock state is not match the num_qubits.")

        if self._kind == "jordan-wigner":
            qubit_state = copy.deepcopy(fock_state)
        else:
            qubit_state = list(map(int, (self._encoding_matrix @ fock_state) % 2))
        return qubit_state


    def qubit_to_fock_state(self, qubit_state: list[int]) -> list[int]:
        """
        Get the fock state from the qubit state
        
        Parameters
        ----------
        qubit_state : list[int]
            qubit state

        Returns
        -------
        fock_state : list[int]
            fock state

        """
        if len(qubit_state) != self._num_qubits:
            raise ValueError("dimension of the qubit state is not match the num_qubits.")

        if self._kind == "jordan-wigner":
            fock_state = copy.deepcopy(qubit_state)
        else:
            fock_state = list(map(int, (self._encoding_matrix_inv @ qubit_state) % 2))
        return fock_state
