# conservation-spectral-topology-c

A C11 implementation of the **Conservation-Spectral-Topology (CST)** framework — a unified mathematical thesis that every system can be understood as a **conserved quantity flowing over a topological space, analyzed spectrally**.

No dependencies beyond standard C11 + libm. Builds with `-Wall -Wextra -std=c11`.

## The CST Thesis

The CST framework rests on three pillars:

1. **Conservation Laws** — Every physical (and many abstract) system has conserved quantities governed by symmetries (Noether's theorem). Energy, mass, charge, momentum — the bookkeeping of what flows where.

2. **Topological Substrate** — These quantities flow over a topological space with structure: connected components, loops, voids. The shape of the space constrains what's possible. Betti numbers characterize this shape.

3. **Spectral Analysis** — The Laplacian of the topology encodes connectivity, bottlenecks, and phase transitions in its eigenvalues. The spectral gap is the system's "health" — a shrinking gap signals an approaching phase transition.

The synthesis: **a CSTSystem** places a conserved quantity on a topological space and uses spectral analysis to monitor system health, detect phase transitions, and verify conservation.

## Architecture

```
include/cst_core.h     — Public API header
src/conservation.c     — Conservation laws, Noether charge verification
src/topology.c         — Simplicial complexes, Betti numbers, persistent homology
src/spectral.c         — Sparse Laplacian, power iteration eigenvalues, Cheeger constant
src/cst_unified.c      — CSTSystem: unified conservation + topology + spectral analysis
tests/test_cst.c       — 37 tests covering all components
```

### Component Map to CST Thesis

| Component | Maps to | Key Function |
|-----------|---------|-------------|
| `ConservedQuantity` | Conservation | `conservation_check()` |
| `NoetherMap` | Symmetry ↔ Conservation | `discrete_noether_check()` |
| `SimplicialComplex` | Topological Space | `compute_betti_numbers()` |
| `PersistenceDiagram` | Topology over scale | `persistence_via_filtration()` |
| `SparseMatrix` (Laplacian) | Spectral substrate | `laplacian_from_complex()` |
| `power_iteration` | Spectral decomposition | `spectral_gap()` |
| `CSTSystem` | **Unified framework** | `cst_verify()`, `cst_phase_detect()`, `cst_predict()` |

## Building

```bash
make          # builds libcstcore.a
make test     # builds and runs all tests
make clean    # remove build artifacts
```

## API Overview

### Conservation

```c
ConservedQuantity q = {CST_QTYPE_ENERGY, .value=100, .flux=5, .source=10, .sink=5};
bool ok = conservation_check(&q, 1e-9);  // |flux - (source-sink)| ≤ tol?
```

### Discrete Noether

```c
// Verify charge conservation under diffusion on a graph
bool ok = discrete_noether_check(adjacency, charges, n, dt, tol);
```

### Topology

```c
SimplicialComplex sc = build_my_complex();
BettiNumbers bn = compute_betti_numbers(&sc);
// bn.b0 = connected components, bn.b1 = loops, bn.b2 = voids

PersistenceDiagram pd = persistence_via_filtration(&sc, vertex_values, n);
// pd.pairs[i] = {birth, death, dimension}
```

### Spectral

```c
SparseMatrix L = laplacian_from_complex(&sc);
double evals[k], evecs[k*n];
power_iteration(&L, k, evals, evecs, 1000, 1e-8);
double gap = spectral_gap(evals, k);
double h = cheeger_constant(fiedler_vector, n);
```

### Unified CST

```c
CSTSystem sys = { ... };  // wire up quantity + topology + eigenvalues
bool conserved  = cst_verify(&sys);
double health   = cst_spectral_health(&sys);    // spectral gap
CstPhase phase  = cst_phase_detect(&sys);        // Normal/Warning/Critical
double distance = cst_predict(&sys);             // 0=transition imminent, 1=safe
```

## Test Coverage

37 tests across all components:

- **Conservation (6):** trivial conserved, source-sink balance, leak detection, discrete Noether on triangle/path/single-vertex graphs
- **Topology (9):** filled/unfilled triangle, tetrahedron, octahedral sphere, cycle graph, disconnected components, persistent homology
- **Spectral (7):** Laplacian construction, eigenvalue verification for path/complete graphs, spectral gap, Cheeger constant
- **Unified (8):** verification, spectral health, phase detection (Normal/Warning/Critical), transition prediction
- **Edge cases (7):** empty system, single vertex, null inputs, degenerate filtration, zero tolerance

## Mathematics

### Conservation
dQ/dt = sources − sinks (within tolerance)

### Topology (mod-2 homology)
B\_k = dim(C\_k) − rank(d\_k) − rank(d\_{k+1})

where d\_k are boundary operators with d² = 0.

### Spectral
L = D − A (graph Laplacian)

λ₂ (Fiedler value) controls connectivity; the spectral gap λ₂ − λ₁ signals phase transitions.

### Phase Detection
- **Normal:** conservation holds, spectral gap large
- **Warning:** conservation violation OR small spectral gap
- **Critical:** conservation violation AND tiny spectral gap

## License

MIT
