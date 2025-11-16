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
    n = N // 3  

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



def compress_graph(adj, cluster_nodes):
    """
    Build a compact adjacency matrix for the induced subgraph on cluster_nodes.
    Returns (compressed_adj_matrix, index_map).
    """
    idx = {node: i for i, node in enumerate(cluster_nodes)}
    n = len(cluster_nodes)

    comp_adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if adj[cluster_nodes[i]][cluster_nodes[j]] == 1:
                comp_adj[i][j] = 1

    return comp_adj, idx


def hamiltonian_bonus(graph, start, end):
    """
    Return True if a Hamiltonian* path exists from start to end using bitmask DP,
    searching only inside the connected component of start; otherwise return False.
    """
    # Extract connected component containing `start`
    cluster_nodes = dfs(graph, start)

    # If end is not reachable from start, no Hamiltonian* path is possible
    if end not in cluster_nodes:
        return False

    # Compress graph to 0..n-1 indexing
    comp_adj, idx = compress_graph(graph, cluster_nodes)
    n = len(cluster_nodes)

    start_l = idx[start]
    end_l = idx[end]

    # dp[mask][v] = path uses exactly vertices in mask and ends at v
    dp = [[False]*n for _ in range(1 << n)]
    dp[1 << start_l][start_l] = True

    # Fill DP table
    for mask in range(1 << n):
        if not (mask & (1 << start_l)):  # must include start
            continue

        for v in range(n):
            if not dp[mask][v]:
                continue

            for u in range(n):
                if not comp_adj[v][u]:   # no edge
                    continue
                if mask & (1 << u):      # already visited
                    continue

                next_mask = mask | (1 << u)

                # End can only be visited last
                if u == end_l and next_mask != (1 << n) - 1:
                    continue

                dp[next_mask][u] = True

    # Full mask means all nodes in cluster visited exactly once
    full_mask = (1 << n) - 1

    # Check whether Hamiltonian path exists
    return dp[full_mask][end_l]

