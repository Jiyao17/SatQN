
from typing import List, Tuple, Dict

def get_gsp(gs):
    pass


def get_visibility(
    sat: List[int], 
    gs: List[int], 
    gsp: List[Tuple[int]]) -> Dict[int, Tuple[int]]:
    pass



if __name__ == "__main__":
    sat = list(range(0, 5))
    gs = list(range(0, 10))
    gsp = [(0, 1), (1, 2), (1, 3)]
    get_visibility(sat, gs, gsp)