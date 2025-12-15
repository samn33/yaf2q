import copy
import random
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from io import BytesIO
from typing import Self


class TernaryTreeSpec:
    """
    Specification of ternary tree
        
    Attributes
    ----------
    num_qubits : int
        number of qubits
    indices : list[int]
        node indices list
    edges : dict[int,tuple[int,str]]
        dictionary data with child node indices as keys, parent node indices and Pauli strings as values
    _parent_dict : dict[int,int]
        dictionary data with child indices as keys, parent indices as values
    _mask : list[int]
        list with a size three times the number of nodes in the graph,
        that the position corresponding to the value of the parent dict is 1 and the other positions are 0.

    Notes
    -----
    Child index is equal to the graph node index, parent index is 3 times the graph node index
    plus the Pauli operator index. Pauli operator index is defined as X=>0,Y=>1,Z=>2.

    """
    def __init__(self, indices: list[int] = [], edges: dict[int,int] = {}):
        """
        Constructor of TernaryTreeSpec
        
        Parameters
        ----------
        indices : list[int]
            node indices list
        edges : dict[int,tuple[int,str]]
            dictionary data with child node indices as keys, parent node indices and Pauli strings as values

        """
        self._indices = indices

        pauli = {'X': 0, 'Y': 1, 'Z': 2}
        self._parent_dict = {}
        for k,v in edges.items():
            self._parent_dict[k] = v[0] * 3 + pauli[v[1]]

        self._mask = [0] * ((len(indices)-1) * 3)
        if self._mask == []:
            return
        for v in self._parent_dict.values():
            self._mask[v] = 1

        if not self.is_valid():
            raise ValueError("the ternary tree spec is invalid.")


    @property
    def num_qubits(self) -> int:
        """ Get the number of qubits """
        return len(self._indices)

    
    @property
    def indices(self) -> int:
        """ Get the indices """
        return self._indices

    
    @property
    def edges(self) -> dict[int,int]:
        """
        Get the edges of the TernaryTreeSpec

        Parameters
        ----------
        None

        Reterns
        -------
        edges : dict[int,tuple[int,str]]
            dictionary data with child node numbers as keys, parent node numbers and Pauli strings as values

        """
        pauli_inv = {0: 'X', 1: 'Y', 2: 'Z'}
        edges = {}
        for k,v in self._parent_dict.items():
            edges[k] = (v // 3, pauli_inv[v % 3])
        return edges
        

    @classmethod
    def random(cls, num_indices: int) -> Self:
        """
        Get a TernaryTreeSpec randomly
        
        Parameters
        ----------
        num_indices : int
            number of indices

        Returns
        -------
        ttspec : TernaryTreeSpec
            ternary tree specification

        """
        ttspec = TernaryTreeSpec()

        ttspec._indices = list(range(num_indices))
        random.shuffle(ttspec._indices)

        ttspec._mask = [0] * ((num_indices-1) * 3)
        for i in range(1,num_indices):
            while True:
                r = random.randint(0, i*3-1)
                if ttspec._mask[r] == 1:
                    continue
                else:
                    ttspec._parent_dict[i] = r
                    ttspec._mask[r] = 1
                    break

        return ttspec

            
    def _to_string_rust(self) -> str:
        """
        Get the string of the TernaryTreeSpec aimed at interfacing with rust

        Parameters
        ----------
        None

        Returns
        -------
        s : str
            string of the TernaryTreeSpec

        Examples
        --------
        If the ternary tree specification is indices: [1, 0, 3, 2] and edges: {1: (0, 'Z'), 2: (1, 'X'), 3: (0, 'Y')},
        then it returns following string, "1 0 3 2\n1 2\n2 3\n3 1"

        """
        s = ""
        for i in self._indices:
            s += (str(i) + " ")
        s = s.strip()
        for k,v in self._parent_dict.items():
            s += ("\n" + str(k) + " " + str(v))
        return s

    
    def __str__(self) -> str:
        return self.to_string()

    def to_string(self) -> str:
        """
        Get the string of TernaryTreeSpec
        
        Parameters
        ----------
        None

        Returns
        -------
        s : str
            string of the TernaryTreeSpec

        """
        s = "indices:"
        s += (str(self._indices))
        s += ", edges:"
        s += (str(self.edges))
        return s[:-1]


    def __eq__(self, other) -> bool:
        return self._indices == other._indices and self._parent_dict == other._parent_dict
        

    def swap_indices(self) -> Self:
        """
        get the TernaryTreeSpec that swapped two indices randomly

        Parameters
        ----------
        None

        Returns
        -------
        ttspec : TernaryTreeSpec
            TernaryTreeSpec that swaped two indices randomly

        """
        ttspec = copy.deepcopy(self)
        num_indices = len(ttspec._indices)
        a = random.randint(0, num_indices-1)
        while True:
            b = random.randint(0, num_indices-1)
            if b == a:
                continue
            else:
                break
        ttspec._indices[a], ttspec._indices[b] = ttspec._indices[b], ttspec._indices[a]

        return ttspec


    def change_parent(self) -> Self:
        """
        get the TernaryTreeSpec that parent node corresponding to the child node has been changed randomly
        
        Parameters
        ----------
        None

        Returns
        -------
        ttspec : TernaryTreeSpec
            TernaryTreeSpec that parent node corresponding to the child node has been changed randomly

        """
        ttspec = copy.deepcopy(self)
        num_indices = len(ttspec._indices)

        while True:
            child = random.randint(1, num_indices-1)
            mask_sum = sum(ttspec._mask[:child*3])
            if mask_sum == child * 3:
                continue
            else:
                break

        parent_org = ttspec._parent_dict[child]
        while True:
            parent = random.randint(0, child*3-1)
            if parent == parent_org or ttspec._mask[parent] == 1:
                continue
            else:
                ttspec._parent_dict[child] = parent
                ttspec._mask[parent_org] = 0
                ttspec._mask[parent] = 1
                break

        return ttspec


    def is_valid(self) -> bool:
        """
        check validity

        Parameters
        ----------
        None

        Returns
        -------
        bool
            the TernaryTreeSpec is valid or not

        """
        if len(self._indices) != len(set(self._indices)):
            return False

        edges = self.edges
        for k,v in edges.items():
            if k <= v[0]:
                return False

        child_list = list(edges.keys())
        parent_list = list(edges.values())

        if len(child_list) != len(set(self._indices)) - 1:
            return False
        elif min(child_list) != 1 or max(child_list) != len(self._indices) - 1:
            return False
        elif len(parent_list) != len(set(parent_list)):
            return False

        return True


    def _make_graph(self) -> nx.Graph:
        """
        make the graph of ternary tree

        Parameters
        ----------
        None

        Returns
        -------
        G : nx.Graph
            graph object of the ternary tree

        """
        indices = self._indices
        edges = self.edges
        
        G = nx.Graph()
        e_dict = {}
        for i in range(len(indices)):
            e_dict[i] = {}
        for k,v in edges.items():
            e_dict[v[0]][v[1]] = k

        cnt = len(indices)
        for k,v in e_dict.items():
            for p in ('X','Y','Z'):
                if p in v:
                    G.add_edge(indices[k], indices[v[p]], label=p)
                else:
                    G.add_edge(indices[k], cnt, label=p)
                    G.nodes[cnt]["style"] = "invisible"
                    G.nodes[cnt]["width"] = 0
                    G.nodes[cnt]["height"] = 0
                    cnt += 1
        return G

    
    def show(self):
        """
        show the ternary tree
        
        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        G = self._make_graph()
        P = nx.nx_pydot.to_pydot(G)
        P.set_rankdir("LR")
        png_data = P.create_png()
        png_stream = BytesIO(png_data)

        img = mpimg.imread(png_stream, format='svg')
        plt.tight_layout()
        plt.imshow(img)
        plt.axis('off')
        plt.show()
