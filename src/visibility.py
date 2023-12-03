
from typing import List, Tuple, Dict, NewType
import random
import math

GSPair = NewType('GSPair', Tuple[int, int])


def get_gsp(gs: List[int], num: int) -> List[GSPair]:
    """
    generate groud station pairs
    """

    gsp = []
    # all possible pairs
    for i in range(0, len(gs)):
        for j in range(i + 1, len(gs)):
            gsp.append((gs[i], gs[j]))
    # randomly select num pairs
    gsp = random.sample(gsp, num)
    
    return gsp


def get_visibility(sat: List[int], gsp: List[GSPair]) \
    -> Dict[int, GSPair]:
    
    visibility = {}
    for s in sat:
        visible_num = random.randint(0, len(gsp) // 2)
        visible_gs = random.sample(gsp, visible_num)
        visibility[s] = visible_gs

    return visibility



if __name__ == "__main__":
    SAT_NUM = 5
    GS_NUM = 10
    GSP_NUM = GS_NUM * (GS_NUM - 1) // 2
    sat = list(range(0, SAT_NUM))
    gs = list(range(0, GS_NUM))
    gsp = get_gsp(gs, GSP_NUM)

    visibility = get_visibility(sat, gsp)
    print(visibility)
