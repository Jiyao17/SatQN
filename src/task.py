
from typing import List, Tuple, Dict
import random

from visibility import get_gsp, get_visibility, GSPair
from gurobipy import Model, GRB, quicksum


class SatAssign:
    def __init__(self,
        sat: List[int],
        gs: List[int],
        gsp: List[GSPair], 
        fidelity: Dict[Tuple[int, List[GSPair]], float],
        Fth: Dict[Tuple[int, GSPair], float],
        weight: Dict[Tuple[int, GSPair], float],
        ) -> None:

        self.S = sat
        self.G = gs
        self.P = gsp
        self.F = fidelity
        self.Fth = Fth
        self.w = weight

        self.model = Model('SA')
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
        obj = 0
        for s in self.S:
            for p in self.P:
                # if self.F[s, p] > self.Fth[p]:
                    obj += self.w[s, p] * self.x[s, p]    
        self.model.setObjective(obj, GRB.MAXIMIZE)
        
        self.model.update()
        self.model.optimize()

    def add_constrs(self, 
        L: List[int],
        R: List[int],
        T: List[int],
        Fth: Dict[Tuple[int, GSPair], float],
        ):
        
        self.add_constr_1a(L)
        self.add_constr_1b(R)
        self.add_constr_1c(T)
        self.add_constr_1d(Fth)

    def add_constr_1a(self, L: List[int]):
        """
        each ground station pair requires at least L[j] entanglments
        """
        assert len(L) == len(self.P)

        for j, p in enumerate(self.P):
            constr = 0
            for s in self.S:
                constr += self.x[s, p]

            self.model.addConstr(constr <= L[j])

    def add_constr_1b(self, R: List[int]):
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

    def add_constr_1c(self, T: List[int]):
        """
        each satellite has T[i] transmitters
        """
        assert len(T) == len(self.S)

        for i, s in enumerate(self.S):
            constr = 0
            for p in self.P:
                constr += self.x[s, p]

            self.model.addConstr(constr <= T[i])

    def add_constr_1d(self, Fth: Dict[Tuple[int, GSPair], float]):
        """
        Entanglements should satisfy the fidelity requirement
        """
        # assert len(Fth.items()) == len(self.F.items())

        for s in self.S:
            for p in self.P:
                # when p is visible to s
                if self.F[s, p] > 0:
                    self.model.addConstr(
                        self.F[s, p] >= Fth[p] * self.x[s, p]
                    )


class SatEntDist:
    def __init__(self,
        sat: List[int],
        gs: List[int],
        gsp: List[GSPair], 
        fidelity: Dict[Tuple[int, List[GSPair]], float],
        Fth: Dict[Tuple[int, GSPair], float],
        weight: Dict[Tuple[int, GSPair], float],
        ) -> None:

        self.S = sat
        self.G = gs
        self.P = gsp
        self.F = fidelity
        self.Fth = Fth
        self.w = weight

        self.model = Model('SA')
        # add decision variables
        # x[sat][gsp] = 1 if sat is assigned to gsp
        self.x = {}
        self.y = {}
        for s in self.S:
            for p in self.P:
                self.x[s, p] = self.model.addVar(
                    vtype=GRB.BINARY, 
                    name='x_{}_{}'.format(s, p)
                    )
                self.y[s, p] = self.model.addVar(
                    vtype=GRB.INTEGER, 
                    name='y_{}_{}'.format(s, p)
                    )
        
        assert len(self.w) == len(self.F) == len(self.x) == len(self.y)
   
    def solve(self, ):        
        obj = 0
        for s in self.S:
            for p in self.P:
                # if self.F[s, p] > self.Fth[p]:
                    obj += self.w[s, p] * self.x[s, p] * self.y[s, p]  
        self.model.setObjective(obj, GRB.MAXIMIZE)
        
        self.model.update()
        self.model.optimize()

    def add_constrs(self, 
        L: List[int],
        R: List[int],
        T: List[int],
        Fth: Dict[Tuple[int, GSPair], float],
        GM: List[int],
        SM: List[int],
        ):
        
        self.add_constr_2a(L)
        self.add_constr_2b(R)
        self.add_constr_2c(T)
        self.add_constr_2d(Fth)
        self.add_constr_2e(GM)
        self.add_constr_2f(SM)


    def add_constr_2a(self, L: List[int]):
        """
        each ground station pair requires at least L[j] entanglments
        """
        assert len(L) == len(self.P)

        for j, p in enumerate(self.P):
            constr = 0
            for s in self.S:
                constr += self.x[s, p]

            self.model.addConstr(constr <= L[j])

    def add_constr_2b(self, R: List[int]):
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

    def add_constr_2c(self, T: List[int]):
        """
        each satellite has T[i] transmitters
        """
        assert len(T) == len(self.S)

        for i, s in enumerate(self.S):
            constr = 0
            for p in self.P:
                constr += self.x[s, p]

            self.model.addConstr(constr <= T[i])

    def add_constr_2d(self, Fth: Dict[Tuple[int, GSPair], float]):
        """
        Entanglements should satisfy the fidelity requirement
        """
        # assert len(Fth.items()) == len(self.F.items())

        for s in self.S:
            for p in self.P:
                # when p is visible to s
                if self.F[s, p] > 0:
                    self.model.addConstr(
                        self.F[s, p] >= Fth[p] * self.x[s, p]
                    )

    def add_constr_2e(self, GM: List[int]):
        """
        Ground station memory constraint
        """
        assert len(GM) == len(self.G)

        for k, g in enumerate(self.G):
            constr = 0
            for s in self.S:
                for p in self.P:
                    if g in p:
                        constr += self.x[s, p] * self.y[s, p]

            self.model.addConstr(constr <= GM[k])

    def add_constr_2f(self, SM: List[int]):
        """
        Satellite memory constraint
        """
        assert len(SM) == len(self.S)

        for i, s in enumerate(self.S):
            constr = 0
            for p in self.P:
                constr += self.x[s, p] * self.y[s, p]

            self.model.addConstr(constr <= SM[i])

    def add_constr_2g(self, SM: List[int]):
        """
        Value constraint of y
        """
        for i, s in enumerate(self.S):
            for p in self.P:
                self.model.addConstr(self.y[s, p] >= 0)
                self.model.addConstr(self.y[s, p] <= SM[i])

def test_sa():
    random.seed(2)
    SAT_NUM = 20
    GS_NUM = 5
    GSP_NUM = 5
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

    L = [SAT_NUM] * len(gsp)
    # randomly assign L[j] = 1 for 1/10 of gsp
    # for i in random.sample(range(0, len(gsp)), len(gsp) // 10):
    #     L[i] = 1
    R = [SAT_NUM] * len(gs)
    T = [GSP_NUM] * len(sat)
    Fth = {}
    for p in gsp:
        Fth[p] = 0.9

    sa = SatAssign(sat, gs, gsp, fideliy, Fth, weight)
    sa.add_constrs(L, R, T, Fth)
    sa.solve()

def test_sed():
    # random.seed(1)
    SAT_NUM = 5
    GS_NUM = 10
    GSP_NUM = GS_NUM * 2
    sat = list(range(0, SAT_NUM))
    gs = list(range(0, GS_NUM))
    gsp = get_gsp(gs, GSP_NUM)
    SM = [random.randint(50, 51) for _ in range(0, SAT_NUM)]
    GM = [random.randint(50, 51) for _ in range(0, GS_NUM)]

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

    L = [SAT_NUM] * len(gsp)
    # randomly assign L[j] = 1 for 1/10 of gsp
    # for i in random.sample(range(0, len(gsp)), len(gsp) // 10):
    #     L[i] = 1
    R = [SAT_NUM] * len(gs)
    T = [GSP_NUM] * len(sat)
    Fth = {}
    for p in gsp:
        Fth[p] = 0.9

    sa = SatEntDist(sat, gs, gsp, fideliy, Fth, weight)
    sa.add_constrs(L, R, T, Fth, GM, SM)
    sa.solve()


if __name__ == '__main__':
    # test_sa()
    test_sed()