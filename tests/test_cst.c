#include "cst_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT_EQ_DBL(a, b, tol) do { \
    if (fabs((a) - (b)) > (tol)) { \
        printf("  FAIL: %s:%d: %.6f != %.6f (tol %.2e)\n", __FILE__, __LINE__, (double)(a), (double)(b), (double)(tol)); \
        tests_failed++; return; \
    } \
} while(0)

#define ASSERT_TRUE(cond) do { \
    if (!(cond)) { \
        printf("  FAIL: %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        tests_failed++; return; \
    } \
} while(0)

#define ASSERT_FALSE(cond) do { \
    if (cond) { \
        printf("  FAIL: %s:%d: expected false: %s\n", __FILE__, __LINE__, #cond); \
        tests_failed++; return; \
    } \
} while(0)

#define TEST(name) static void test_##name(void)

#define RUN(name) do { \
    printf("  %-50s", #name); \
    int f_before = tests_failed; \
    test_##name(); \
    if (tests_failed == f_before) { printf("PASS\n"); tests_passed++; } \
} while(0)

/* Helper: build a simple triangle complex */
static SimplicialComplex make_triangle(int a, int b, int c)
{
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 3;
    sc.edges = (int *)malloc(6 * sizeof(int));
    sc.edges[0]=a; sc.edges[1]=b;
    sc.edges[2]=a; sc.edges[3]=c;
    sc.edges[4]=b; sc.edges[5]=c;
    sc.num_edges = 3;
    sc.triangles = (Triangle *)malloc(sizeof(Triangle));
    sc.triangles[0].v[0]=a; sc.triangles[0].v[1]=b; sc.triangles[0].v[2]=c;
    sc.num_triangles = 1;
    return sc;
}

/* ================================================================== */
/* CONSERVATION TESTS */
/* ================================================================== */

TEST(conserved_trivial)
{
    ConservedQuantity q = {CST_QTYPE_ENERGY, 100.0, 0.0, 0.0, 0.0};
    ASSERT_TRUE(conservation_check(&q, 1e-9));
}

TEST(conserved_with_source_sink)
{
    ConservedQuantity q = {CST_QTYPE_ENERGY, 100.0, 5.0, 10.0, 5.0};
    ASSERT_TRUE(conservation_check(&q, 1e-9));
}

TEST(leaking_detected)
{
    ConservedQuantity q = {CST_QTYPE_ENERGY, 100.0, 5.0, 10.0, 3.0};
    /* flux=5, source-sink=7, violation of 2 */
    ASSERT_FALSE(conservation_check(&q, 1.0));
    ASSERT_TRUE(conservation_check(&q, 3.0)); /* within larger tolerance */
}

TEST(discrete_noether_simple_graph)
{
    /* Triangle graph: fully connected 3 vertices, uniform charge */
    double adj[] = {0, 1, 1,
                    1, 0, 1,
                    1, 1, 0};
    double charge[] = {1.0, 2.0, 3.0};
    /* Diffusion on symmetric adjacency should conserve total */
    ASSERT_TRUE(discrete_noether_check(adj, charge, 3, 0.01, 1e-9));
}

TEST(discrete_noether_path_graph)
{
    /* Path: 0-1-2 */
    double adj[] = {0, 1, 0,
                    1, 0, 1,
                    0, 1, 0};
    double charge[] = {1.0, 2.0, 3.0};
    ASSERT_TRUE(discrete_noether_check(adj, charge, 3, 0.1, 1e-9));
}

TEST(discrete_noether_single_vertex)
{
    double adj[] = {0};
    double charge[] = {5.0};
    ASSERT_TRUE(discrete_noether_check(adj, charge, 1, 0.1, 1e-9));
}

/* ================================================================== */
/* TOPOLOGY TESTS */
/* ================================================================== */

TEST(betti_triangle)
{
    /* Filled triangle = disk: B0=1, B1=0, B2=0 */
    SimplicialComplex sc = make_triangle(0, 1, 2);
    BettiNumbers bn = compute_betti_numbers(&sc);
    ASSERT_EQ_DBL(bn.b0, 1, 0.01);
    ASSERT_EQ_DBL(bn.b1, 0, 0.01); /* filled triangle, no loop */
    ASSERT_EQ_DBL(bn.b2, 0, 0.01);
    simplicial_complex_free(&sc);
}

TEST(betti_two_triangles_shared_edge)
{
    /* Two filled triangles sharing edge (1,2) = disk: B0=1, B1=0 */
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 4;
    sc.num_edges = 5;
    sc.edges = (int *)malloc(10 * sizeof(int));
    sc.edges[0]=0; sc.edges[1]=1;
    sc.edges[2]=0; sc.edges[3]=2;
    sc.edges[4]=1; sc.edges[5]=2;
    sc.edges[6]=1; sc.edges[7]=3;
    sc.edges[8]=2; sc.edges[9]=3;
    sc.num_triangles = 2;
    sc.triangles = (Triangle *)malloc(2 * sizeof(Triangle));
    sc.triangles[0].v[0]=0; sc.triangles[0].v[1]=1; sc.triangles[0].v[2]=2;
    sc.triangles[1].v[0]=1; sc.triangles[1].v[1]=2; sc.triangles[1].v[2]=3;
    BettiNumbers bn = compute_betti_numbers(&sc);
    ASSERT_EQ_DBL(bn.b0, 1, 0.01);
    ASSERT_EQ_DBL(bn.b1, 0, 0.01); /* filled disk, no loop */
    simplicial_complex_free(&sc);
}

TEST(betti_single_edge)
{
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 2;
    sc.num_edges = 1;
    sc.edges = (int *)malloc(2 * sizeof(int));
    sc.edges[0] = 0; sc.edges[1] = 1;
    BettiNumbers bn = compute_betti_numbers(&sc);
    ASSERT_EQ_DBL(bn.b0, 1, 0.01);
    ASSERT_EQ_DBL(bn.b1, 0, 0.01);
    simplicial_complex_free(&sc);
}

TEST(betti_disconnected)
{
    /* Two disconnected edges */
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 4;
    sc.num_edges = 2;
    sc.edges = (int *)malloc(4 * sizeof(int));
    sc.edges[0]=0; sc.edges[1]=1;
    sc.edges[2]=2; sc.edges[3]=3;
    BettiNumbers bn = compute_betti_numbers(&sc);
    ASSERT_EQ_DBL(bn.b0, 2, 0.01);
    ASSERT_EQ_DBL(bn.b1, 0, 0.01);
    simplicial_complex_free(&sc);
}

TEST(betti_tetrahedron)
{
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 4;
    sc.num_edges = 6;
    sc.edges = (int *)malloc(12 * sizeof(int));
    int e[][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    for (int i = 0; i < 6; i++) { sc.edges[2*i]=e[i][0]; sc.edges[2*i+1]=e[i][1]; }
    sc.num_triangles = 4;
    sc.triangles = (Triangle *)malloc(4 * sizeof(Triangle));
    int t[][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
    for (int i = 0; i < 4; i++) { sc.triangles[i].v[0]=t[i][0]; sc.triangles[i].v[1]=t[i][1]; sc.triangles[i].v[2]=t[i][2]; }
    sc.num_tetrahedra = 1;
    sc.tetrahedra = (Tetrahedron *)malloc(sizeof(Tetrahedron));
    sc.tetrahedra[0].v[0]=0; sc.tetrahedra[0].v[1]=1; sc.tetrahedra[0].v[2]=2; sc.tetrahedra[0].v[3]=3;

    BettiNumbers bn = compute_betti_numbers(&sc);
    ASSERT_EQ_DBL(bn.b0, 1, 0.01);   /* connected */
    ASSERT_EQ_DBL(bn.b1, 0, 0.01);   /* no loops */
    ASSERT_EQ_DBL(bn.b2, 0, 0.01);   /* no voids (filled by tetrahedron) */
    simplicial_complex_free(&sc);
}

TEST(betti_sphere_8_triangle)
{
    /* Octahedron surface = S² sphere, 6 vertices, 8 triangles, 12 edges */
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 6;
    sc.num_edges = 12;
    sc.edges = (int *)malloc(24 * sizeof(int));
    /* Octahedron edges: top(0) connected to equator(1,2,3,4), bottom(5) connected to equator, equator cycle */
    int edges[][2] = {
        {0,1},{0,2},{0,3},{0,4},  /* top to equator */
        {5,1},{5,2},{5,3},{5,4},  /* bottom to equator */
        {1,2},{2,3},{3,4},{4,1}   /* equator */
    };
    for (int i = 0; i < 12; i++) { sc.edges[2*i]=edges[i][0]; sc.edges[2*i+1]=edges[i][1]; }
    sc.num_triangles = 8;
    sc.triangles = (Triangle *)malloc(8 * sizeof(Triangle));
    int tris[][3] = {
        {0,1,2},{0,2,3},{0,3,4},{0,4,1},  /* top hemi */
        {5,2,1},{5,3,2},{5,4,3},{5,1,4}   /* bottom hemi */
    };
    for (int i = 0; i < 8; i++) { sc.triangles[i].v[0]=tris[i][0]; sc.triangles[i].v[1]=tris[i][1]; sc.triangles[i].v[2]=tris[i][2]; }

    BettiNumbers bn = compute_betti_numbers(&sc);
    ASSERT_EQ_DBL(bn.b0, 1, 0.01);
    ASSERT_EQ_DBL(bn.b1, 0, 0.01);  /* sphere has no tunnels */
    ASSERT_EQ_DBL(bn.b2, 1, 0.01);  /* sphere has 1 void */
    simplicial_complex_free(&sc);
}

TEST(betti_cycle_graph)
{
    /* Square: 4 vertices, 4 edges forming a cycle */
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 4;
    sc.num_edges = 4;
    sc.edges = (int *)malloc(8 * sizeof(int));
    sc.edges[0]=0; sc.edges[1]=1;
    sc.edges[2]=1; sc.edges[3]=2;
    sc.edges[4]=2; sc.edges[5]=3;
    sc.edges[6]=3; sc.edges[7]=0;
    BettiNumbers bn = compute_betti_numbers(&sc);
    ASSERT_EQ_DBL(bn.b0, 1, 0.01);
    ASSERT_EQ_DBL(bn.b1, 1, 0.01); /* one loop */
    simplicial_complex_free(&sc);
}

TEST(betti_triangle_unfilled)
{
    /* Triangle edges only (no face) = loop */
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 3;
    sc.num_edges = 3;
    sc.edges = (int *)malloc(6 * sizeof(int));
    sc.edges[0]=0; sc.edges[1]=1;
    sc.edges[2]=1; sc.edges[3]=2;
    sc.edges[4]=2; sc.edges[5]=0;
    BettiNumbers bn = compute_betti_numbers(&sc);
    ASSERT_EQ_DBL(bn.b0, 1, 0.01);
    ASSERT_EQ_DBL(bn.b1, 1, 0.01); /* unfilled triangle = 1 loop */
    simplicial_complex_free(&sc);
}

TEST(persistence_basic)
{
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 3;
    sc.num_edges = 3;
    sc.edges = (int *)malloc(6 * sizeof(int));
    sc.edges[0]=0; sc.edges[1]=1;
    sc.edges[2]=0; sc.edges[3]=2;
    sc.edges[4]=1; sc.edges[5]=2;

    double vals[] = {0.0, 1.0, 2.0};
    PersistenceDiagram pd = persistence_via_filtration(&sc, vals, 3);
    ASSERT_TRUE(pd.count >= 1); /* At least one pair */
    /* First component born at 0.0, dies when edge connects */
    persistence_diagram_free(&pd);
    simplicial_complex_free(&sc);
}

/* ================================================================== */
/* SPECTRAL TESTS */
/* ================================================================== */

TEST(laplacian_single_vertex)
{
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 1;
    sc.num_edges = 0;
    sc.edges = NULL;
    SparseMatrix L = laplacian_from_complex(&sc);
    ASSERT_EQ_DBL(L.nrows, 1, 0.1);
    ASSERT_EQ_DBL(L.nnz, 1, 0.1);
    /* L[0,0] = 0 (isolated vertex, degree 0) */
    ASSERT_EQ_DBL(L.values[0], 0.0, 1e-12);
    sparse_matrix_free(&L);
    simplicial_complex_free(&sc);
}

TEST(laplacian_edge)
{
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 2;
    sc.num_edges = 1;
    sc.edges = (int *)malloc(2 * sizeof(int));
    sc.edges[0] = 0; sc.edges[1] = 1;
    SparseMatrix L = laplacian_from_complex(&sc);
    /* L = [[1,-1],[-1,1]] */
    ASSERT_EQ_DBL(L.nnz, 4, 0.1);
    sparse_matrix_free(&L);
    simplicial_complex_free(&sc);
}

TEST(laplacian_eigenvalues_path)
{
    /* Path graph 0-1-2: L = [[1,-1,0],[-1,2,-1],[0,-1,1]] */
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 3;
    sc.num_edges = 2;
    sc.edges = (int *)malloc(4 * sizeof(int));
    sc.edges[0]=0; sc.edges[1]=1;
    sc.edges[2]=1; sc.edges[3]=2;
    SparseMatrix L = laplacian_from_complex(&sc);

    double evals[3];
    double evecs[9];
    int rc = power_iteration(&L, 3, evals, evecs, 1000, 1e-8);
    ASSERT_EQ_DBL(rc, 0, 0.1);

    /* Sort eigenvalues */
    for (int i = 0; i < 2; i++)
        for (int j = i+1; j < 3; j++)
            if (evals[j] < evals[i]) { double t = evals[i]; evals[i] = evals[j]; evals[j] = t; }

    /* Known: eigenvalues of path-3 are 0, 1, 3 */
    ASSERT_EQ_DBL(evals[0], 0.0, 0.01);
    ASSERT_EQ_DBL(evals[1], 1.0, 0.01);
    ASSERT_EQ_DBL(evals[2], 3.0, 0.01);

    sparse_matrix_free(&L);
    simplicial_complex_free(&sc);
}

TEST(spectral_gap_connected)
{
    double evals[] = {0.0, 1.0, 3.0};
    double gap = spectral_gap(evals, 3);
    ASSERT_EQ_DBL(gap, 1.0, 0.1);
}

TEST(spectral_gap_single)
{
    double evals[] = {0.0};
    double gap = spectral_gap(evals, 1);
    ASSERT_EQ_DBL(gap, 0.0, 1e-12);
}

TEST(laplacian_complete_graph)
{
    /* K3: complete graph on 3 vertices */
    SimplicialComplex sc = make_triangle(0, 1, 2);
    SparseMatrix L = laplacian_from_complex(&sc);

    double evals[3];
    double evecs[9];
    power_iteration(&L, 3, evals, evecs, 1000, 1e-8);

    /* Sort */
    for (int i = 0; i < 2; i++)
        for (int j = i+1; j < 3; j++)
            if (evals[j] < evals[i]) { double t = evals[i]; evals[i] = evals[j]; evals[j] = t; }

    /* K3 eigenvalues: 0, 3, 3 */
    ASSERT_EQ_DBL(evals[0], 0.0, 0.01);
    ASSERT_EQ_DBL(evals[1], 3.0, 0.01);

    sparse_matrix_free(&L);
    simplicial_complex_free(&sc);
}

TEST(cheeger_nonzero)
{
    /* Fiedler vector of path-3: [1, 0, -1] roughly */
    double fiedler[] = {1.0, 0.0, -1.0};
    double h = cheeger_constant(fiedler, 3);
    ASSERT_TRUE(h > 0.0);
}

TEST(degenerate_eigenvalues)
{
    /* K4: complete graph on 4 vertices — eigenvalues 0, 4, 4, 4 (degenerate) */
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 4;
    sc.num_edges = 6;
    sc.edges = (int *)malloc(12 * sizeof(int));
    int e[][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    for (int i = 0; i < 6; i++) { sc.edges[2*i]=e[i][0]; sc.edges[2*i+1]=e[i][1]; }
    SparseMatrix L = laplacian_from_complex(&sc);

    double evals[4];
    double evecs[16];
    int rc = power_iteration(&L, 4, evals, evecs, 2000, 1e-10);
    ASSERT_EQ_DBL(rc, 0, 0.1);

    /* Sort eigenvalues */
    for (int i = 0; i < 3; i++)
        for (int j = i+1; j < 4; j++)
            if (evals[j] < evals[i]) { double t = evals[i]; evals[i] = evals[j]; evals[j] = t; }

    /* K4 eigenvalues: 0, 4, 4, 4 */
    ASSERT_EQ_DBL(evals[0], 0.0, 0.01);
    ASSERT_EQ_DBL(evals[1], 4.0, 0.5); /* degenerate, may be less precise */
    ASSERT_EQ_DBL(evals[3], 4.0, 0.5);

    /* Spectral gap should be 4.0 */
    double gap = spectral_gap(evals, 4);
    ASSERT_EQ_DBL(gap, 4.0, 0.5);

    sparse_matrix_free(&L);
    simplicial_complex_free(&sc);
}

/* ================================================================== */
/* UNIFIED CST TESTS */
/* ================================================================== */

TEST(cst_verify_good)
{
    CSTSystem sys;
    memset(&sys, 0, sizeof(sys));
    sys.quantity = (ConservedQuantity){CST_QTYPE_ENERGY, 100.0, 5.0, 10.0, 5.0};
    sys.conservation_tolerance = 1e-9;
    ASSERT_TRUE(cst_verify(&sys));
}

TEST(cst_verify_bad)
{
    CSTSystem sys;
    memset(&sys, 0, sizeof(sys));
    sys.quantity = (ConservedQuantity){CST_QTYPE_ENERGY, 100.0, 5.0, 10.0, 3.0};
    sys.conservation_tolerance = 0.1;
    ASSERT_FALSE(cst_verify(&sys));
}

TEST(cst_spectral_health_triangle)
{
    SimplicialComplex sc = make_triangle(0, 1, 2);
    SparseMatrix L = laplacian_from_complex(&sc);

    CSTSystem sys;
    memset(&sys, 0, sizeof(sys));
    sys.topology = sc;
    sys.num_eigenvalues = 3;
    sys.eigenvalues = (double *)malloc(3 * sizeof(double));
    double evecs[9];
    power_iteration(&L, 3, sys.eigenvalues, evecs, 1000, 1e-8);

    double health = cst_spectral_health(&sys);
    ASSERT_TRUE(health > 0.0);

    sparse_matrix_free(&L);
    cst_system_free(&sys);
}

TEST(cst_phase_normal)
{
    CSTSystem sys;
    memset(&sys, 0, sizeof(sys));
    sys.quantity = (ConservedQuantity){CST_QTYPE_ENERGY, 100.0, 0.0, 0.0, 0.0};
    sys.conservation_tolerance = 1e-9;
    double evals[] = {0.0, 3.0, 3.0};
    sys.eigenvalues = evals;
    sys.num_eigenvalues = 3;
    sys.warning_gap = 0.5;
    sys.critical_gap = 0.1;
    ASSERT_TRUE(cst_phase_detect(&sys) == CST_PHASE_NORMAL);
}

TEST(cst_phase_warning)
{
    CSTSystem sys;
    memset(&sys, 0, sizeof(sys));
    sys.quantity = (ConservedQuantity){CST_QTYPE_ENERGY, 100.0, 0.0, 0.0, 0.0};
    sys.conservation_tolerance = 1e-9;
    double evals[] = {0.0, 0.01, 3.0}; /* very small gap */
    sys.eigenvalues = evals;
    sys.num_eigenvalues = 3;
    sys.warning_gap = 0.5;
    sys.critical_gap = 0.1;
    ASSERT_TRUE(cst_phase_detect(&sys) == CST_PHASE_WARNING);
}

TEST(cst_phase_critical)
{
    CSTSystem sys;
    memset(&sys, 0, sizeof(sys));
    sys.quantity = (ConservedQuantity){CST_QTYPE_ENERGY, 100.0, 5.0, 0.0, 0.0}; /* leaking */
    sys.conservation_tolerance = 0.01;
    double evals[] = {0.0, 0.001, 3.0}; /* tiny gap */
    sys.eigenvalues = evals;
    sys.num_eigenvalues = 3;
    sys.warning_gap = 0.5;
    sys.critical_gap = 0.1;
    ASSERT_TRUE(cst_phase_detect(&sys) == CST_PHASE_CRITICAL);
}

TEST(cst_predict_far)
{
    CSTSystem sys;
    memset(&sys, 0, sizeof(sys));
    double evals[] = {0.0, 5.0, 10.0};
    sys.eigenvalues = evals;
    sys.num_eigenvalues = 3;
    sys.warning_gap = 0.5;
    double dist = cst_predict(&sys);
    ASSERT_TRUE(dist > 0.9); /* far from transition */
}

TEST(cst_predict_near)
{
    CSTSystem sys;
    memset(&sys, 0, sizeof(sys));
    double evals[] = {0.0, 0.001, 10.0};
    sys.eigenvalues = evals;
    sys.num_eigenvalues = 3;
    sys.warning_gap = 0.5;
    double dist = cst_predict(&sys);
    ASSERT_TRUE(dist < 0.1); /* near transition */
}

/* ================================================================== */
/* EDGE CASES */
/* ================================================================== */

TEST(empty_system)
{
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    BettiNumbers bn = compute_betti_numbers(&sc);
    ASSERT_EQ_DBL(bn.b0, 0, 0.01);
    ASSERT_EQ_DBL(bn.b1, 0, 0.01);
}

TEST(single_vertex_laplacian)
{
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 1;
    SparseMatrix L = laplacian_from_complex(&sc);
    ASSERT_EQ_DBL(L.nrows, 1, 0.1);
    ASSERT_EQ_DBL(L.ncols, 1, 0.1);
    sparse_matrix_free(&L);
}

TEST(null_persistence)
{
    PersistenceDiagram pd = persistence_via_filtration(NULL, NULL, 0);
    ASSERT_TRUE(pd.count == 0);
}

TEST(degenerate_filtration)
{
    /* All same values */
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 2;
    sc.num_edges = 1;
    sc.edges = (int *)malloc(2 * sizeof(int));
    sc.edges[0]=0; sc.edges[1]=1;
    double vals[] = {1.0, 1.0};
    PersistenceDiagram pd = persistence_via_filtration(&sc, vals, 2);
    ASSERT_TRUE(pd.count >= 1);
    persistence_diagram_free(&pd);
    simplicial_complex_free(&sc);
}

TEST(power_iteration_single)
{
    SimplicialComplex sc;
    memset(&sc, 0, sizeof(sc));
    sc.num_vertices = 1;
    SparseMatrix L = laplacian_from_complex(&sc);
    double evals[1], evecs[1];
    int rc = power_iteration(&L, 1, evals, evecs, 100, 1e-8);
    ASSERT_EQ_DBL(rc, 0, 0.1);
    ASSERT_EQ_DBL(evals[0], 0.0, 0.1);
    sparse_matrix_free(&L);
}

TEST(conservation_zero_tolerance)
{
    ConservedQuantity q = {CST_QTYPE_ENERGY, 100.0, 0.0, 0.0, 0.0};
    ASSERT_TRUE(conservation_check(&q, 0.0));
}

TEST(cst_null_returns)
{
    ASSERT_FALSE(cst_verify(NULL));
    ASSERT_EQ_DBL(cst_spectral_health(NULL), 0.0, 1e-12);
    ASSERT_TRUE(cst_phase_detect(NULL) == CST_PHASE_CRITICAL);
    ASSERT_EQ_DBL(cst_predict(NULL), 0.0, 1e-12);
}

/* ================================================================== */

int main(void)
{
    printf("CST Core Test Suite\n");
    printf("====================\n\n");

    printf("[Conservation]\n");
    RUN(conserved_trivial);
    RUN(conserved_with_source_sink);
    RUN(leaking_detected);
    RUN(discrete_noether_simple_graph);
    RUN(discrete_noether_path_graph);
    RUN(discrete_noether_single_vertex);

    printf("\n[Topology]\n");
    RUN(betti_triangle);
    RUN(betti_two_triangles_shared_edge);
    RUN(betti_single_edge);
    RUN(betti_disconnected);
    RUN(betti_tetrahedron);
    RUN(betti_sphere_8_triangle);
    RUN(betti_cycle_graph);
    RUN(betti_triangle_unfilled);
    RUN(persistence_basic);

    printf("\n[Spectral]\n");
    RUN(laplacian_single_vertex);
    RUN(laplacian_edge);
    RUN(laplacian_eigenvalues_path);
    RUN(spectral_gap_connected);
    RUN(spectral_gap_single);
    RUN(laplacian_complete_graph);
    RUN(cheeger_nonzero);
    RUN(degenerate_eigenvalues);

    printf("\n[Unified CST]\n");
    RUN(cst_verify_good);
    RUN(cst_verify_bad);
    RUN(cst_spectral_health_triangle);
    RUN(cst_phase_normal);
    RUN(cst_phase_warning);
    RUN(cst_phase_critical);
    RUN(cst_predict_far);
    RUN(cst_predict_near);

    printf("\n[Edge Cases]\n");
    RUN(empty_system);
    RUN(single_vertex_laplacian);
    RUN(null_persistence);
    RUN(degenerate_filtration);
    RUN(power_iteration_single);
    RUN(conservation_zero_tolerance);
    RUN(cst_null_returns);

    printf("\n====================\n");
    printf("Results: %d passed, %d failed, %d total\n",
           tests_passed, tests_failed, tests_passed + tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
