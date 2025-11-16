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
    This function takes the original full adjacency matrix `adj` and the list 
    of vertices `cluster_nodes` found by DFS, and builds a smaller adjacency 
    matrix that contains ONLY the vertices inside this cluster.

    Why we need it:
      - The bitmask DP requires vertices to be labeled from 0 to n-1.
      - `cluster_nodes` contains arbitrary original labels (e.g., [7, 12, 3]).
      - We must remap:  original_label → compressed_label.
      - Then we create an n×n adjacency matrix for the cluster.

    Output:
      comp_adj : compressed adjacency matrix of the cluster
      idx      : dictionary mapping original node → compressed index
    """
    # map original node → 0..n-1
    idx = {node: i for i, node in enumerate(cluster_nodes)}
    n = len(cluster_nodes)

    # compressed adjacency
    comp_adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if adj[cluster_nodes[i]][cluster_nodes[j]] == 1:
                comp_adj[i][j] = 1

    return comp_adj, idx


def hamiltonian_path_bitmask_cluster(adj, start, end, cluster_nodes):
    """
    This function finds a Hamiltonian path from 'start' to 'end', but ONLY inside 
    the connected component given by `cluster_nodes`.

    Algorithm steps:
      1. Compress the graph:
         - Convert cluster nodes to 0..n-1 indexing
         - Extract the adjacency matrix of only the cluster

      2. Dynamic Programming (Held–Karp style):
         - dp[mask][v] = True if there is a path using exactly the nodes in `mask`
                         and ending at vertex `v`.
         - mask is a bitmask representing which nodes have been visited.

      3. Transition:
         - From dp[mask][v], try to extend to a neighbor `u`
         - Skip u if:
             (a) it is not adjacent to v   <-- BASIC OPERATION
             (b) it is already in mask
             (c) u is the end node but it's not the last step
       4. Final check:
         - If dp[(1<<n)-1][end] is True, a Hamiltonian path exists.

       5. Reconstruct the path by backtracking through the DP table.
    """
    # Step 1: compress graph to local (0..n-1)
    comp_adj, idx = compress_graph(adj, cluster_nodes)
    n = len(cluster_nodes)

    start_l = idx[start]
    end_l   = idx[end]

    # dp[mask][v] = whether we can reach v using nodes in mask
    dp = [[False]*n for _ in range(1 << n)]
    dp[1 << start_l][start_l] = True

    # Fill DP
    for mask in range(1 << n):
        if not (mask & (1 << start_l)): 
            continue  # path must include start

        for v in range(n):
            if not dp[mask][v]:
                continue

            # Try extending to neighbors
            for u in range(n):
                if not comp_adj[v][u]:
                    continue
                if mask & (1 << u):
                    continue  # already used

                # Cannot visit end unless this is the last step
                next_mask = mask | (1 << u)
                if u == end_l and next_mask != (1 << n) - 1:
                    continue

                dp[next_mask][u] = True

    full_mask = (1 << n) - 1

    # Check if solution exists
    if not dp[full_mask][end_l]:
        return None

    # Reconstruct path
    path = []
    mask = full_mask
    cur = end_l

    for _ in range(n):
        path.append(cluster_nodes[cur])  # convert back to original numbers
        prev_mask = mask ^ (1 << cur)

        if prev_mask == 0:
            break

        # find predecessor
        for u in range(n):
            if dp[prev_mask][u] and comp_adj[u][cur]:
                # end must be last
                if u == end_l and prev_mask != (1 << end_l):
                    continue
                cur = u
                mask = prev_mask
                break

    return path[::-1]

