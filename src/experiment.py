import time
import matplotlib.pyplot as plt
import random
from graph_construction import generate_tricky_graph
from solution import hamiltonian_optimized, hamiltonian_naive, hamiltonian_bonus

random.seed(31)

def measure_runtime(func, n, rounds=10):
    """
    Measure average runtime across different random graphs.
    Graph construction is excluded from timing.
    """
    times = []
    for _ in range(rounds):
        # generate a new graph (excluded from timing)
        graph, start, end = generate_tricky_graph(n)

        # measure only Hamiltonian* call
        t0 = time.perf_counter()
        func(graph, start, end)
        t1 = time.perf_counter()

        times.append(t1 - t0)
    return sum(times) / len(times)


def main():
    ns = [4, 5, 6, 7, 8, 9]   # ≥ 4, at least 5 different values
    avg_times = []

    for n in ns:
        print(f"\nRunning experiments for n = {n}...")

        avg_time = measure_runtime(hamiltonian_naive, n, rounds=10)
        avg_times.append(avg_time)

        print(f"Average time for n={n}: {avg_time:.7f} seconds")

    plt.figure(figsize=(8,5))
    plt.plot(ns, avg_times, marker='o', linestyle='-')
    plt.xlabel("n (size of each component)")
    plt.ylabel("Average Execution Time (seconds)")
    plt.title("Hamiltonian* Naive: n vs Execution Time")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("naive_n_vs_time.png")
    plt.show()


if __name__ == "__main__":
    main()
