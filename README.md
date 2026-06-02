# conservation-spectral-topology-c

You're monitoring a distributed system. CPU is at 40%. Memory at 70%. Disk at 90%. Are you healthy?

The answer depends on the **topology** of your dependencies, not just the numbers.

Disk at 90% means nothing if it's an isolated log server. Disk at 90% on the node that three critical services route through? That's a fuse burning toward the powder keg. The numbers didn't change. The **topology** did.

That's the starting point for Conservation-Spectral Topology: every system is a conserved quantity flowing over a topological space, and you can analyze it spectrally. This library is a clean C11 implementation of that idea — no dependencies beyond libm.

## The Insight

Every monitoring system checks numbers. But numbers without topology are just noise. CST asks a different question: **what flows where, and is the flow pattern stable?**

Three ideas, one framework:

1. **Conservation** — Energy, mass, charge, request rate, money. Every system has quantities that should be accounted for. If you're putting 100 requests/sec in and getting 93 out, something's wrong. `dQ/dt = sources − sinks`. This is just bookkeeping, but it's the kind of bookkeeping that catches leaks.

2. **Topology** — Your system has a shape. Connected components (isolated clusters), loops (circular dependencies), voids (enclosed failure domains). That shape constrains what can happen. A cycle in your dependency graph isn't just a code smell — it's a topological feature (Betti number β₁ = 1) that affects how failures propagate.

3. **Spectral analysis** — The graph Laplacian `L = D − A` encodes your system's connectivity in its eigenvalues. The second-smallest eigenvalue (the Fiedler value) tells you how well-connected the graph is. When it shrinks toward zero, your system is losing cohesion — nodes are becoming isolated, bottlenecks are forming. **The spectral gap predicts cascading failures before they happen.**

Put them together and you get something powerful: verify conservation on a topology, then monitor the spectral gap to know if the topology itself is healthy. Each layer catches what the others miss. Conservation violations detect leaks. Betti numbers detect structural changes (a new dependency cycle appeared, a component got isolated). The spectral gap detects *incipient* structural decay — the topology hasn't changed yet, but it's about to.

## Why the Spectral Gap Matters for YOUR Systems

The spectral gap is the distance between the smallest and second-smallest eigenvalue of the graph Laplacian.

- **Large gap** → Well-connected graph. Robust. Information (or requests, or energy) flows freely. A single node failure doesn't partition the system.
- **Small gap** → Fragile. There's a bottleneck or near-partition. The graph is close to splitting into disconnected components.
- **Zero gap** → Already disconnected. Game over.

Here's the thing: the spectral gap shrinks *continuously* as the system degrades. It doesn't jump from "fine" to "broken." That means you can watch it trend downward and **raise an alarm before the failure cascade starts.**

This library computes it in ~50 lines of sparse power iteration. No GPU, no Python, no NumPy. Pure C. It runs in microseconds on graphs with thousands of nodes. You can poll it every second in a monitoring loop and watch the gap evolve in real time.

## 10 Lines That Do Real Work

```c
#include "cst_core.h"

// Build your dependency graph as a simplicial complex
SimplicialComplex sc = {0};
sc.num_vertices = 4;
int edges[] = {0,1, 0,2, 1,2, 1,3, 2,3};
sc.edges = edges; sc.num_edges = 5;

// Compute the Laplacian and its eigenvalues
SparseMatrix L = laplacian_from_complex(&sc);
double evals[4], evecs[16];
power_iteration(&L, 4, evals, evecs, 1000, 1e-8);

// The spectral gap tells you how robust the topology is
double gap = spectral_gap(evals, 4);
printf("Spectral gap: %.4f — %s\n", gap,
       gap > 0.5 ? "healthy" : gap > 0.1 ? "degrading" : "critical");
```

That's it. Build a graph, compute its Laplacian, extract eigenvalues, read the spectral gap. If you've ever used a graph library, this should feel familiar — it's just that the math is doing more work for you than you might expect.

## The Viral Insight

Every distributed system IS a conserved quantity on a graph.

- Your load balancer routes requests (conserved quantity: request count) across backends (topology: a bipartite graph).
- Your power grid moves electricity (conserved: energy) across transmission lines (topology: a weighted graph).
- Your neural network propagates activations (approximately conserved: gradient norm) across layers (topology: a DAG).
- Your social network carries influence (approximately conserved: attention) across follows (topology: a directed graph).

You were already doing CST. You just didn't know it.

The "conservation" part is just: does what goes in come out? The "topology" part is just: what connects to what? The "spectral" part is just: how healthy are those connections? You've been checking these things separately. CST says they're the same problem.

## Building

```bash
make          # builds libcstcore.a static library
make test     # builds and runs all 37 tests
make clean    # remove build artifacts
```

Requires C11 + libm. No other dependencies. Compiles with `-Wall -Wextra -std=c11 -O2`. Link with `-lcstcore -lm`.

## What's Inside

```
include/cst_core.h     — Public API (everything below)
src/conservation.c     — Conservation law checking, discrete Noether verification
src/topology.c         — Simplicial complexes, Betti numbers, persistent homology
src/spectral.c         — Sparse Laplacian, power iteration eigenvalues, Cheeger constant
src/cst_unified.c      — CSTSystem: unified conservation + topology + spectral health
tests/test_cst.c       — 37 tests, no framework, just assertions and a runner
```

### Conservation

Check that flux equals sources minus sinks. Detect leaks. Verify charge conservation under diffusion on a graph (discrete Noether).

```c
ConservedQuantity q = {CST_QTYPE_ENERGY, .value=100, .flux=5, .source=10, .sink=5};
bool ok = conservation_check(&q, 1e-9);  // true: 5 == 10 - 5

// Check that diffusion on a graph preserves total charge
bool noether_ok = discrete_noether_check(adjacency, charges, n, dt, tol);
```

### Topology

Build simplicial complexes (vertices, edges, triangles, tetrahedra). Compute Betti numbers — β₀ is connected components, β₁ is loops, β₂ is enclosed voids. Run persistent homology to see which topological features survive across scales.

```c
BettiNumbers bn = compute_betti_numbers(&sc);
// bn.b0 = components, bn.b1 = loops, bn.b2 = voids

PersistenceDiagram pd = persistence_via_filtration(&sc, vertex_values, n);
// pd.pairs[i] = {birth, death, dimension}
```

An unfilled triangle (three edges, no face) has β₁ = 1 — one loop. Fill the face and β₁ drops to 0. A sphere (octahedron surface) has β₀ = 1, β₁ = 0, β₂ = 1 — one connected component, no tunnels, one enclosed void. The topology code computes all of this via mod-2 homology with Gaussian elimination.

### Spectral

Sparse graph Laplacian, power iteration with deflation for top-k eigenvalues, spectral gap computation, and Cheeger constant estimation from the Fiedler vector.

```c
SparseMatrix L = laplacian_from_complex(&sc);
double evals[k], evecs[k*n];
power_iteration(&L, k, evals, evecs, 1000, 1e-8);
double gap = spectral_gap(evals, k);
double h = cheeger_constant(fiedler_vector, n);
```

The Cheeger constant measures the "bottleneck-ness" of your graph — the minimum cut that separates the graph into two roughly equal pieces. It's bounded by the Fiedler value (Cheeger's inequality), so the spectral gap gives you a quick estimate of how bad the worst bottleneck is.

### Unified CST

Wire everything together into a `CSTSystem`. Get a single health score.

```c
CSTSystem sys = { ... };
bool conserved    = cst_verify(&sys);            // is conservation holding?
double health     = cst_spectral_health(&sys);    // spectral gap
CstPhase phase    = cst_phase_detect(&sys);       // Normal / Warning / Critical
double distance   = cst_predict(&sys);            // 0 = transition imminent, 1 = safe
```

Phase detection combines conservation and spectral health:
- **Normal:** Conservation holds, spectral gap is large.
- **Warning:** Conservation violated OR spectral gap is small.
- **Critical:** Conservation violated AND spectral gap is tiny. Something is very wrong.

`cst_predict()` returns a normalized distance-to-transition. Near 0 means the system is about to undergo a phase change (graph partition, cascade failure). Near 1 means it's far from trouble. This is your early warning system.

## Test Coverage

37 tests, zero dependencies beyond the library itself:

- **Conservation (6):** trivial conserved, source-sink balance, leak detection, discrete Noether on triangle/path/single-vertex graphs
- **Topology (9):** filled/unfilled triangle, tetrahedron, octahedral sphere (β₂=1), cycle graph (β₁=1), disconnected components, persistent homology
- **Spectral (7):** Laplacian construction, eigenvalue verification for path and complete graphs, spectral gap, Cheeger constant
- **Unified (8):** verification, spectral health, phase detection across all three phases, transition prediction (far and near)
- **Edge cases (7):** empty system, single vertex, null inputs, degenerate filtration, zero tolerance

Every test uses known mathematical results — the eigenvalues of a path-3 graph are {0, 1, 3}, K₃ has eigenvalues {0, 3, 3}, an octahedron surface has β₂ = 1. This isn't testing that the code runs. It's testing that the math is correct. If you port this to Rust, Go, or Python, these tests are your oracle — get the same numbers, and your port is correct.

## Reference

### Conservation

$$\frac{dQ}{dt} = \text{sources} - \text{sinks} \quad (\text{within tolerance } \varepsilon)$$

Discrete Noether: total charge is conserved under adjacency-weighted diffusion.

### Topology (mod-2 homology)

$$\beta_k = \dim(C_k) - \text{rank}(d_k) - \text{rank}(d_{k+1})$$

where $d_k$ are boundary operators over $\mathbb{Z}/2\mathbb{Z}$ with $d^2 = 0$.

### Spectral

$$L = D - A$$

$\lambda_2$ (Fiedler value) controls connectivity. The spectral gap $\lambda_2 - \lambda_1$ signals phase transitions.

**Cheeger's inequality:**

$$\frac{\lambda_2}{2} \leq h \leq \sqrt{2\lambda_2}$$

where $h$ is the Cheeger constant — the minimum edge cut that separates the graph into two roughly equal parts.

### Phase Detection

- **Normal:** conservation holds ∧ gap > warning threshold
- **Warning:** ¬conservation ∨ gap < warning threshold
- **Critical:** ¬conservation ∧ gap < critical threshold

### Transition Prediction

$$\text{distance} = \frac{\lambda_2}{\lambda_2 + \text{warning\_gap}}$$

0 = transition imminent (spectral gap collapsed). 1 = far from transition.

## License

MIT
