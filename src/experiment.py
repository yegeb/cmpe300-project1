import time
import matplotlib.pyplot as plt
import random
from src.graph_construction import generate_tricky_graph  # adjust name if needed
from src.solution import hamiltonian_naive, hamiltonian_optimized

random.seed(31)

def measure_runtime(func, *args, rounds=10):
    """Return average execution time over multiple rounds."""
    times = []
    for _ in range(rounds):
        start_t = time.perf_counter()
        func(*args)
        end_t = time.perf_counter()
        times.append(end_t - start_t)
    return sum(times) / len(times)

def main():
    ns = [4, 5, 6, 7, 8]   # ≥4, at least five different values
    avg_times = []

    for n in ns:
        print(f"\nRunning experiments for n = {n}...")
        N = 3 * n

        # generate graph once per round group (excluded from timing)
        graph, start, end = generate_tricky_graph(n)

        # measure runtime for 10 runs (excluding construction)
        avg_time = measure_runtime(hamiltonian_optimized, graph, start, end, rounds=10)
        avg_times.append(avg_time)

        print(f"Average time for n={n}: {avg_time:.7f} seconds")

    # plot results
    plt.figure(figsize=(8,5))
    plt.plot(ns, avg_times, marker='o', linestyle='-', color='b')
    plt.xlabel("n (size of each component)")
    plt.ylabel("Average Execution Time (seconds)")
    plt.title("Hamiltonian* Optimized: n vs Execution Time")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("optimized_n_vs_time.png")
    plt.show()

if __name__ == "__main__":
    main()
