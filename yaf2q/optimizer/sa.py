import copy
import math
import random
import statistics
from dataclasses import dataclass, field

from yaf2q.ternary_tree_spec import TernaryTreeSpec

DEF_INIT_SAMPLING = 10
DEF_NUM_STEPS = 10
DEF_COOLING_FACTOR = 1.0

@dataclass
class SAParams:
    """
    Parameters of simulated annealing

    Attributes
    ----------
    init_sampling : int, default 10
        sampling number for determining the initial temperature
    num_steps : int, default 10
        anealing steps
    cooling_factor : float, default 1.0
        coooling factor for annealing
    seed : int
        seed for random generation

    Notes
    -----
    Annealing scheduling curve is T = (x + 1)^(-a), where T is temperature, x is annealing step (int), and a is cooling factor.

    """
    init_sampling: int    = field(default=DEF_INIT_SAMPLING, init=True)
    num_steps: int        = field(default=DEF_NUM_STEPS, init=True)
    cooling_factor: float = field(default=DEF_COOLING_FACTOR, init=True)
    seed: int             = field(default=None, init=True)


class SA:
    """
    Oprtimizer using simulated annealing

    Attributes
    ----------
    num_qubits : int
        number of qubits
    objective_func : func
        objective function
    params : SAParams
        parameters for simulated annealing
    verbose : bool
        displey detailed information

    """
    def __init__(self, num_qubits: int, objective_func, params: SAParams, verbose: bool = False):
        """
        Constructor of SA

        Parameters
        ----------
        num_qubits : int
            number of qubits
        objective_func : func
            objective function
        params : SAParams
            parameters for simulated annealing
        verbose : bool
            displey detailed information

        """
        self.num_qubits = num_qubits
        self.objective_func = objective_func
        self.params = params
        self.verbose = verbose

        if self.params.seed is not None:
            random.seed(self.params.seed)


    def optimize(self):
        """
        optimize a ternary tree specification

        Parameters
        ----------
        None

        Returns
        -------
        ttspec_opt : TernaryTreeSpec
            optimized ternary tree specification

        """
        num_qubits = self.num_qubits
        
        ttspec_init = TernaryTreeSpec.random(num_qubits)
        ttspec = copy.deepcopy(ttspec_init)
        ttspec_opt = copy.deepcopy(ttspec_init)

        obj = self.objective_func(ttspec = ttspec_init)
        obj_opt = obj

        obj_list = []
        for _ in range(self.params.init_sampling):
            obj = self.objective_func(ttspec = TernaryTreeSpec.random(num_qubits))
            obj_list.append(obj)
        temper_init = statistics.stdev(obj_list)
        
        num_steps = self.params.num_steps
        cooling_factor = self.params.cooling_factor
        for i in range(num_steps):
            temper = temper_init * pow(((i + 1) * num_steps) / num_steps, -cooling_factor)
            if self.verbose:
                print("\r", f"step = {i + 1}/{num_steps}, objective value = {obj_opt}", end="")
                if i == num_steps - 1:
                    print()

            for j in range(num_qubits * 2):
                if j // num_qubits == 0:
                    ttspec_new = ttspec.swap_indices()
                else:
                    ttspec_new = ttspec.change_parent()

                obj_new = self.objective_func(ttspec = ttspec_new)
                if obj_new < obj or random.random() < math.exp((obj - obj_new) / temper):
                    if obj_new < obj_opt:
                        obj_opt = obj_new
                        ttspec_opt = copy.deepcopy(ttspec_new)
                    obj = obj_new
                    ttspec = copy.deepcopy(ttspec_new)

        return ttspec_opt
