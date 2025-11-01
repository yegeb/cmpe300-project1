from itertools import combinations, permutations
from typing import List, Sequence


def hamiltonian_check(H: List[List[int]], perm: Sequence[int]) -> bool:
    """
    Algorithm 5: Hamiltonian*Check(H, perm)
    Returns True iff every consecutive pair in perm is adjacent in H.
    """
    n = len(H)  # number of rows/cols of H
    for i in range(n - 1):  # i = 0 to n-2
        if H[perm[i]][perm[i + 1]] == 0:
            return False  # early exit on first missing edge
    return True  # all consecutive pairs are adjacent


def all_permutations(H: List[List[int]], s: int, t: int) -> bool:
    """
    Algorithm 4: AllPermutations(H, s, t)
    Enumerate all permutations of {0..n-1} with perm(0)=s and perm(n-1)=t,
    returning True if any permutation forms a Hamiltonian* path in H.
    """
    n = len(H)  # number of rows/cols of H  
    middle = [i for i in range(n) if (i != s and i != t)]
    for mid_perm in permutations(middle):
        perm = (s, *mid_perm, t)
        if hamiltonian_check(H, perm):
            return True
    return False


def hamiltonian_naive(graph: List[List[int]], start: int, end: int) -> bool:
    """
    Algorithm 3: Hamiltonian*Naive(3n×3n Graph, start, end)
    Tries every size-n subset S ⊆ V that contains {start, end}, builds the induced
    submatrix H in a fixed order L, and checks all permutations with fixed endpoints.
    """
    N = len(graph)
    n = N // 3  # N is 3n by problem definition

    # V ← {0, 1, ..., 3n-1}
    vertices = list(range(N))

    # for each subset S of size n that contains start and end
    others = [v for v in vertices if v not in (start, end)]
    for middle in combinations(others, n - 2):
        # S and fixed-order list L (any fixed order is allowed; use sorted)
        S = set(middle)
        S.add(start)
        S.add(end)
        L = sorted(S)  # "any fixed order"; sorting gives a deterministic one

        # Build H[a][b] ← graph[L[a]][L[b]] for all a, b
        H = [[graph[L[a]][L[b]] for b in range(n)] for a in range(n)]

        # s, t are indices of start and end within L
        s = L.index(start)
        t = L.index(end)

        # If any valid permutation yields a Hamiltonian* path, return True
        if all_permutations(H, s, t):
            return True

    # No subset produced a valid path
    return False

def _component_nodes(graph, src):
    """
    Return the (sorted) list of vertices in the connected component of src.
    Uses Depth-First Search (DFS) on the undirected adjacency matrix.
    """
    N = len(graph)
    visited = [False] * N
    comp = []

    def dfs(u):
        visited[u] = True
        comp.append(u)
        for v, edge in enumerate(graph[u]):
            # If there's an edge and it's unvisited, continue DFS
            if edge and not visited[v]:
                dfs(v)

    dfs(src)
    return sorted(comp)  # ensure deterministic vertex order 

def hamiltonian_optimized(graph: List[List[int]], start: int, end: int) -> bool:
    """
    Optimized Hamiltonian* algorithm that exploits the graph structure.

    Since the 3n×3n graph consists of three disjoint n-node components,
    we only need to consider the component containing `start`.
    If `end` is not in the same component, no Hamiltonian* path can exist.

    Steps:
      1. Find the component C containing `start`.
      2. If `end` ∉ C, return False.
      3. Build the induced adjacency matrix H for C.
      4. Call AllPermutations(H, s, t).
    """
    N = len(graph)
    n = N // 3  # guaranteed by generator

    # Find connected component of start
    L = _component_nodes(graph, start)

    # Must be in the same component
    if end not in L:
        return False

    # Build H[a][b] = graph[L[a]][L[b]]
    H = [[graph[L[a]][L[b]] for b in range(n)] for a in range(n)]

    # Map start and end to indices in L
    s = L.index(start)
    t = L.index(end)

    # Use Algorithm 4 as subroutine
    return all_permutations(H, s, t)

