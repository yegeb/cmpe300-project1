def dfs_from(adj_matrix, start_idx):
    n = len(adj_matrix)
    visited = [False] * n
    order = []
    stack = [start_idx]

    while stack:
        v = stack.pop()
        if visited[v]:
            continue
        visited[v] = True
        order.append(v)

        # Collect neighbors of v (where there's a 1 in the row)
        neighbors = [u for u, edge in enumerate(adj_matrix[v]) if edge]
        # Push in reverse so we pop the smallest index first
        for u in reversed(sorted(neighbors)):
            if not visited[u]:
                stack.append(u)

    return order


def compress_graph(adj, cluster_nodes):
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

graph, start, end = generate_tricky_graph(n)
cluster_nodes = dfs_from(graph, start)
graph, start, end = generate_tricky_graph(n)
if end not in cluster_nodes:
    print("No Hamiltonian path found (end not reachable from start)")
else:
    path = hamiltonian_path_bitmask_cluster(graph, start, i, cluster_nodes)
    if path is None:
        print("No Hamiltonian path found")
    else:
        print("Hamiltonian path found:", path)