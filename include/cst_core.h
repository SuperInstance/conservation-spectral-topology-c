#ifndef CST_CORE_H
#define CST_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================
 * 1. CONSERVATION LAWS
 * ======================================================================== */

typedef enum {
    CST_QTYPE_ENERGY   = 0,
    CST_QTYPE_MASS     = 1,
    CST_QTYPE_MOMENTUM = 2,
    CST_QTYPE_CHARGE   = 3,
    CST_QTYPE_CUSTOM   = 4
} CstQuantityType;

typedef struct {
    CstQuantityType type;
    double value;      /* current amount */
    double flux;       /* dQ/dt observed */
    double source;     /* inflow rate */
    double sink;       /* outflow rate */
} ConservedQuantity;

typedef struct {
    double tolerance;  /* allowed violation */
} ConservationLaw;

typedef struct {
    int symmetry_id;
    int quantity_id;
} NoetherPair;

typedef struct {
    NoetherPair *pairs;
    size_t       count;
} NoetherMap;

/* Check conservation: |flux - (source - sink)| <= tolerance */
bool conservation_check(const ConservedQuantity *q, double tolerance);

/* Discrete Noether charge conservation check on a graph.
   graph: adjacency matrix flattened row-major, n x n.
   charge_per_vertex: array of length n.
   Returns true if total charge is conserved (sum invariant). */
bool discrete_noether_check(const double *adjacency, const double *charge_per_vertex,
                            size_t n, double dt, double tolerance);

/* ========================================================================
 * 2. TOPOLOGICAL SUBSTRATE
 * ======================================================================== */

typedef struct {
    int v[3]; /* vertex indices, sorted */
} Triangle;

typedef struct {
    int v[4]; /* vertex indices, sorted */
} Tetrahedron;

typedef struct {
    int    num_vertices;
    int   *edges;        /* flat [2*num_edges] */
    size_t num_edges;
    Triangle   *triangles;
    size_t      num_triangles;
    Tetrahedron *tetrahedra;
    size_t       num_tetrahedra;
} SimplicialComplex;

/* Chain complex: boundary operators as sparse maps */
typedef struct {
    int    *boundary_map; /* for each simplex, list of boundary simplex indices; -1 terminated */
    size_t  num_simplices;
    int     dimension;    /* 0=vertices, 1=edges, 2=triangles, ... */
} ChainGroup;

typedef struct {
    ChainGroup *groups;
    size_t      num_groups;
} ChainComplex;

/* Betti numbers: B0=components, B1=loops, B2=voids */
typedef struct {
    int b0;
    int b1;
    int b2;
} BettiNumbers;

/* Compute Betti numbers for a simplicial complex */
BettiNumbers compute_betti_numbers(const SimplicialComplex *sc);

/* Persistence pair */
typedef struct {
    double birth;
    double death;     /* INFINITY if never dies */
    int    dimension;
} PersistencePair;

typedef struct {
    PersistencePair *pairs;
    size_t           count;
} PersistenceDiagram;

/* Simple filtration: filter vertices by value, compute persistence */
PersistenceDiagram persistence_via_filtration(const SimplicialComplex *sc,
                                              const double *vertex_values,
                                              size_t num_values);

/* Free functions */
void simplicial_complex_free(SimplicialComplex *sc);
void persistence_diagram_free(PersistenceDiagram *pd);

/* ========================================================================
 * 3. SPECTRAL ANALYSIS
 * ======================================================================== */

typedef struct {
    int    *row_ptr;
    int    *col_idx;
    double *values;
    int     nrows;
    int     ncols;
    int     nnz;
} SparseMatrix;

/* Build sparse matrix (allocates, caller must free with sparse_matrix_free) */
void sparse_matrix_free(SparseMatrix *m);

/* Build graph Laplacian from simplicial complex (1-skeleton) */
SparseMatrix laplacian_from_complex(const SimplicialComplex *sc);

/* Power iteration: compute top-k eigenvalues/eigenvectors.
   eigenvalues: output array of length k
   eigenvectors: output flat array, k rows of length n each (row-major)
   Returns 0 on success. */
int power_iteration(const SparseMatrix *A, int k,
                    double *eigenvalues, double *eigenvectors,
                    int max_iter, double tol);

/* Spectral gap: smallest nonzero eigenvalue minus zero (or gap between smallest two) */
double spectral_gap(const double *eigenvalues, int n);

/* Cheeger constant estimate from Fiedler vector (second eigenvector) */
double cheeger_constant(const double *fiedler, int n);

/* ========================================================================
 * 4. UNIFIED CST FRAMEWORK
 * ======================================================================== */

typedef enum {
    CST_PHASE_NORMAL   = 0,
    CST_PHASE_WARNING  = 1,
    CST_PHASE_CRITICAL = 2
} CstPhase;

typedef struct {
    ConservedQuantity  quantity;
    SimplicialComplex  topology;
    /* Spectral state */
    double  spectral_gap_val;
    double  cheeger_val;
    double *eigenvalues;
    int     num_eigenvalues;
    /* Thresholds */
    double  conservation_tolerance;
    double  warning_gap;
    double  critical_gap;
} CSTSystem;

/* Verify conservation law holds on the topology */
bool cst_verify(const CSTSystem *sys);

/* Spectral health: higher gap = healthier */
double cst_spectral_health(const CSTSystem *sys);

/* Phase detection based on conservation + spectral gap */
CstPhase cst_phase_detect(const CSTSystem *sys);

/* Predict phase transition: returns estimated "distance" to transition (0=imminent) */
double cst_predict(const CSTSystem *sys);

/* Free CST system */
void cst_system_free(CSTSystem *sys);

#ifdef __cplusplus
}
#endif

#endif /* CST_CORE_H */
