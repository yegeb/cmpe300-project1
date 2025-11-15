from itertools import combinations, permutations
from typing import List, Sequence


def hamiltonian_check(H: List[List[int]], perm: Sequence[int]) -> bool:
    """
    Algorithm 5: Hamiltonian*Check(H, perm)
    Returns True iff every consecutive pair in perm is adjacent in H.
    """
    n = len(H)  # number of rows/cols of H
    for i in range(n - 1):  
        if H[perm[i]][perm[i + 1]] == 0:
            return False 
    return True  


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


def hamiltonian_naive(graph, start, end) -> bool:
    """
    Algorithm 3: Hamiltonian*Naive
    Input: 3nx3n adjacency matrix 'graph', start, end
    Output: True if a Hamiltonian* path exists from start to end
    """
    N = len(graph)
    n = N // 3  
    V = list(range(N))  

    # Iterate over all S element of V 
    for S in combinations(V, n):
        if start not in S or end not in S:
            continue

        # L = list of vertices in S 
        L = list(S)   

        # Build H[a][b] = graph[L[a]][L[b]] 
        H = [[graph[L[a]][L[b]] for b in range(n)] for a in range(n)]

        # map start and end to indices in L 
        s = L.index(start)
        t = L.index(end)

        # check permutations 
        if all_permutations(H, s, t):
            return True

    # If no subset worked
    return False

def dfs(graph, src) -> List:
    """
    Depth-First Search on an undirected adjacency matrix.
    Returns the list of all vertices reachable from `src`.
    """
    N = len(graph)
    visited = [False] * N
    stack = [src]
    comp = []

    while stack:
        u = stack.pop()
        if not visited[u]:
            visited[u] = True
            comp.append(u)
            # push neighbors
            for v, edge in enumerate(graph[u]):
                if edge and not visited[v]:
                    stack.append(v)

    return comp


def hamiltonian_optimized(graph: List[List[int]], start: int, end: int) -> bool:
    """
    Algorithm 6: Optimized Hamiltonian* algorithm that exploits the graph structure.
    Since the 3nx3n graph consists of three disjoint n-node components,
    we only need to consider the component containing `start`.
    If `end` is not in the same component, no Hamiltonian* path can exist.
    """

    N = len(graph)
    n = N // 3  # guaranteed by generator

    # Find connected component of start
    L = dfs(graph, start)

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