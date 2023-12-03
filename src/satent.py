
from typing import List, Tuple, Dict
import random

from abc import ABC, abstractmethod

from visibility import get_gsp, get_visibility, GSPair
from gurobipy import Model, GRB, quicksum

from networkx import DiGraph
import networkx as nx

import matplotlib.pyplot as plt


class SatAssign(ABC):
    def __init__(self,
        sat: List[int],
        gs: List[int],
        gsp: List[GSPair], 
        fidelity: Dict[Tuple[int, List[GSPair]], float],
        # Fth: Dict[Tuple[int, GSPair], float],
        weight: Dict[Tuple[int, GSPair], float],
        ) -> None:

        self.S = sat
        self.G = gs
        self.P = gsp
        self.F = fidelity
        self.w = weight

        self.params: Tuple = None

    def set_params(self,
            L: List[int],
            R: List[int],
            T: List[int],
            Fth: Dict[Tuple[int, GSPair], float],
        ) -> None:

        self.params = (L, R, T, Fth)

    @abstractmethod
    def solve(self, ):
        pass

class SatAssignILP(SatAssign):
    """
    Integer linear programming solution, using Gurobi
    """
    def __init__(self, 
        sat: List[int], 
        gs: List[int], 
        gsp: List[GSPair], 
        fidelity: Dict[Tuple[int, List[GSPair]], float], 
        weight: Dict[Tuple[int, GSPair], float]
        ) -> None:
        super().__init__(sat, gs, gsp, fidelity, weight)

        self.model = Model('SA-ILP')
        self.model_built = False

        # add decision variables
        # x[sat][gsp] = 1 if sat is assigned to gsp
        self.x = {}
        for s in self.S:
            for p in self.P:
                self.x[s, p] = self.model.addVar(
                    vtype=GRB.BINARY, 
                    name='x_{}_{}'.format(s, p)
                    )
        assert len(self.w) == len(self.F) == len(self.x)

    def solve(self, ): 
        assert self.model_built

        obj = 0
        for s in self.S:
            for p in self.P:
                # if self.F[s, p] > self.Fth[p]:
                    obj += self.w[s, p] * self.x[s, p]    
        self.model.setObjective(obj, GRB.MAXIMIZE)
        
        self.model.update()
        self.model.optimize()

    def add_constrs(self, ):
        assert self.params is not None

        L, R, T, Fth = self.params
        self._add_constr_1a(L)
        self._add_constr_1b(R)
        self._add_constr_1c(T)
        self._add_constr_1d(Fth)

        self.model_built = True

    def _add_constr_1a(self, L: List[int]):
        """
        each ground station pair requires at least L[j] entanglments
        """
        assert len(L) == len(self.P)

        for j, p in enumerate(self.P):
            constr = 0
            for s in self.S:
                constr += self.x[s, p]

            self.model.addConstr(constr <= L[j])

    def _add_constr_1b(self, R: List[int]):
        """
        each ground station can accept at most R[k] entanglments
        """
        assert len(R) == len(self.G)

        for k, g in enumerate(self.G):
            constr = 0
            for s in self.S:
                for p in self.P:
                    if g in p:
                        constr += self.x[s, p]

            self.model.addConstr(constr <= R[k])

    def _add_constr_1c(self, T: List[int]):
        """
        each satellite has T[i] transmitters
        """
        assert len(T) == len(self.S)

        for i, s in enumerate(self.S):
            constr = 0
            for p in self.P:
                constr += self.x[s, p]

            self.model.addConstr(constr <= T[i])

    def _add_constr_1d(self, Fth: Dict[Tuple[int, GSPair], float]):
        """
        Entanglements should satisfy the fidelity requirement
        """
        # assert len(Fth.items()) == len(self.F.items())

        for s in self.S:
            for p in self.P:
                # when p is visible to s
                if self.F[s, p] > 0:
                    self.model.addConstr(
                        self.F[s, p] >= self.x[s, p] * Fth[p]
                    )


class SatAssignMF(SatAssign):
    """
    Max flow solution, using networkx
    """

    def __init__(self, 
        sat: List[int], 
        gs: List[int], 
        gsp: List[GSPair], 
        fidelity: Dict[Tuple[int, List[GSPair]], float], 
        # Fth: Dict[Tuple[int, GSPair], float], 
        weight: Dict[Tuple[int, GSPair], float]
        ) -> None:
        super().__init__(sat, gs, gsp, fidelity, weight)

        self.graph = DiGraph()
        self.graph_built = False

    def build_graph(self, ):
        assert self.params is not None

        L, R, T, Fth = self.params
        
        # add source and sink
        self.graph.add_node('s', subset=0)
        self.graph.add_node('t', subset=4)

        # add nodes for satellites
        for s in self.S:
            self.graph.add_node(f's_{str(s)}', subset=1)
        # add edges from source to satellites
        for i, s in enumerate(self.S):
            self.graph.add_edge('s', f's_{str(s)}', capacity=T[i])

        # add nodes for ground stations
        for g in self.G:
            self.graph.add_node(f'g_{str(g)}', subset=3)

        # add nodes for ground station pairs
        # add edges from ground station pairs to ground stations
        for j, p in enumerate(self.P):
            node = f'p_{str(p[0])}_{str(p[1])}'
            for num in range(L[j]):
                node_n = node + f'_{str(num)}'
                self.graph.add_node(node_n, subset=2)
                self.graph.add_edge(node_n, f'g_{str(p[0])}', capacity=1)
                self.graph.add_edge(node_n, f'g_{str(p[1])}', capacity=1)

        # add edges from satellites to ground station pairs
        for i, s in enumerate(self.S):
            for j, p in enumerate(self.P):
                if self.w[s, p] > 0:
                    for num in range(L[j]):
                        node_n = f'p_{str(p[0])}_{str(p[1])}_{str(num)}'
                        self.graph.add_edge(f's_{str(s)}', node_n, capacity=2)
        
        # add edges from ground stations to sink
        for k, g in enumerate(self.G):
            self.graph.add_edge(f'g_{str(g)}', 't', capacity=R[k])

        # fidelity constraints
        for s in self.S:
            for p in self.P:
                for num in range(L[j]):
                    node_n = f'p_{str(p[0])}_{str(p[1])}_{str(num)}'
                    if self.F[s, p] < Fth[p] and self.graph.has_edge(f's_{str(s)}', node_n):
                        edge = self.graph.edges[f's_{str(s)}', node_n]
                        edge['capacity'] = 0


        self.graph_built = True

    def solve(self, ):
        assert self.graph_built

        # find max flow from source to sink
        flow_value, flow_dict = nx.maximum_flow(self.graph, 's', 't')

        return flow_value, flow_dict


def test_sa(sat, gs, gsp, fideliy, weight, L, R, T, Fth):

    sa = SatAssignILP(sat, gs, gsp, fideliy, weight)
    sa.set_params(L, R, T, Fth)
    sa.add_constrs()
    sa.solve()

def test_mf(sat, gs, gsp, fideliy, weight, L, R, T, Fth):
    mf = SatAssignMF(sat, gs, gsp, fideliy, weight)
    mf.set_params(L, R, T, Fth)
    mf.build_graph()
    flow_value, flow_dict = mf.solve()

    print(flow_value)
    # save graph with flow
    pos = nx.multipartite_layout(mf.graph, subset_key='subset', align='vertical')
    nx.draw(mf.graph, pos=pos, with_labels=True)
    # convert flow_dict to edge_labels
    edge_labels = {}
    for start, end_dict in flow_dict.items():
        for end, flow in end_dict.items():
            edge_labels[start, end] = flow

    nx.draw_networkx_edge_labels(mf.graph, pos=pos, edge_labels=edge_labels)
    plt.savefig('flow.png')


def test_example():
    sat = list(range(0, 3))
    gs = list(range(0, 3))
    gsp = [(0, 1), (1, 2)]

    fideliy = {
        (0, (0, 1)): 0.9,
        (0, (1, 2)): 0.9,
        (1, (0, 1)): 0.9,
        (1, (1, 2)): 0.9,
        (2, (0, 1)): 0.9,
        (2, (1, 2)): 0.9,
    }

    weight = {
        (0, (0, 1)): 1,
        (0, (1, 2)): 1,
        (1, (0, 1)): 0,
        (1, (1, 2)): 1,
        (2, (0, 1)): 0,
        (2, (1, 2)): 1,
    }

    L = [2, 1]
    R = [1, 2, 1]
    T = [4, 2, 2]
    Fth = {
        (0, 1): 0.8,
        (1, 2): 0.8,
    }

    test_sa(sat, gs, gsp, fideliy, weight, L, R, T, Fth)
    test_mf(sat, gs, gsp, fideliy, weight, L, R, T, Fth)

def test_rand():
    random.seed(2)
    SAT_NUM = 3
    GS_NUM = 3
    GSP_NUM = 2
    sat = list(range(0, SAT_NUM))
    gs = list(range(0, GS_NUM))
    gsp = get_gsp(gs, GSP_NUM)

    visibility = get_visibility(sat, gsp)

    fideliy = {}
    base_fidelity = 0.9
    for s in sat:
        for p in gsp:
            if p in visibility[s]:
                fideliy[s, p] = base_fidelity + random.random() * (1 - base_fidelity)
            else:
                fideliy[s, p] = 0

    weight = {}
    for s in sat:
        for p in gsp:
            if p in visibility[s]:
                weight[s, p] = 1
            else:
                weight[s, p] = 0

    L = [3] * len(gsp)
    # randomly assign L[j] = 1 for 1/10 of gsp
    # for i in random.sample(range(0, len(gsp)), len(gsp) // 10):
    #     L[i] = 1
    R = [5] * len(gs)
    T = [5] * len(sat)
    Fth = {}
    for p in gsp:
        Fth[p] = 0.8


    test_sa(sat, gs, gsp, fideliy, weight, L, R, T, Fth)
    test_mf(sat, gs, gsp, fideliy, weight, L, R, T, Fth)


if __name__ == '__main__':
    # test_example()
    test_rand()