#include "cst_core.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ---- helpers ---- */

static int cmp_int(const void *a, const void *b)
{
    return *(const int *)a - *(const int *)b;
}

/* Comparators for qsort in persistence */
typedef struct { int idx; double val; } idx_val;

static int cmp_idx_val(const void *a, const void *b)
{
    double da = ((const idx_val *)a)->val, db = ((const idx_val *)b)->val;
    return (da > db) - (da < db);
}

typedef struct { size_t idx; double filt; } edge_filt;

static int cmp_edge_filt(const void *a, const void *b)
{
    double da = ((const edge_filt *)a)->filt, db = ((const edge_filt *)b)->filt;
    return (da > db) - (da < db);
}

/* Union-Find for connected components */
static void uf_init(int *parent, int n)
{
    for (int i = 0; i < n; i++) parent[i] = i;
}

static int uf_find(int *parent, int x)
{
    while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
}

static void uf_union(int *parent, int a, int b)
{
    a = uf_find(parent, a);
    b = uf_find(parent, b);
    if (a != b) parent[a] = b;
}

/* Mod-2 rank via Gaussian elimination on a binary matrix.
   matrix: rows × cols, stored row-major, each entry 0 or 1.
   Returns rank over Z/2. */
static int mod2_rank(int rows, int cols, const int *matrix)
{
    if (rows == 0 || cols == 0) return 0;
    /* Overflow check for allocation size */
    size_t alloc_size;
    if (rows > 0 && cols > 0 && (size_t)rows <= SIZE_MAX / (size_t)cols) {
        alloc_size = (size_t)rows * (size_t)cols;
        if (alloc_size > SIZE_MAX / sizeof(int)) return -1;
    } else {
        return -1;
    }
    int *m = (int *)malloc(alloc_size * sizeof(int));
    if (!m) return -1;
    memcpy(m, matrix, rows * cols * sizeof(int));

    int rank = 0;
    for (int col = 0; col < cols && rank < rows; col++) {
        /* Find pivot */
        int pivot = -1;
        for (int row = rank; row < rows; row++) {
            if (m[row * cols + col]) { pivot = row; break; }
        }
        if (pivot < 0) continue;
        /* Swap */
        if (pivot != rank) {
            for (int c = 0; c < cols; c++) {
                int tmp = m[rank * cols + c];
                m[rank * cols + c] = m[pivot * cols + c];
                m[pivot * cols + c] = tmp;
            }
        }
        /* Eliminate */
        for (int row = 0; row < rows; row++) {
            if (row != rank && m[row * cols + col]) {
                for (int c = 0; c < cols; c++)
                    m[row * cols + c] ^= m[rank * cols + c];
            }
        }
        rank++;
    }
    free(m);
    return rank;
}

/* Free */
void simplicial_complex_free(SimplicialComplex *sc)
{
    if (!sc) return;
    free(sc->edges);       sc->edges = NULL;
    free(sc->triangles);   sc->triangles = NULL;
    free(sc->tetrahedra);  sc->tetrahedra = NULL;
}

void persistence_diagram_free(PersistenceDiagram *pd)
{
    if (!pd) return;
    free(pd->pairs); pd->pairs = NULL;
}

/* ---- Betti numbers ---- */
/* Using mod-2 homology with proper rank computation:
   B_k = dim(C_k) - rank(d_k) - rank(d_{k+1}) */

static int count_connected_components(int num_vertices, const int *edges, size_t num_edges)
{
    if (num_vertices == 0) return 0;
    int *parent = (int *)malloc(num_vertices * sizeof(int));
    uf_init(parent, num_vertices);
    for (size_t e = 0; e < num_edges; e++) {
        uf_union(parent, edges[2*e], edges[2*e+1]);
    }
    int count = 0;
    for (int i = 0; i < num_vertices; i++) {
        if (uf_find(parent, i) == i) count++;
    }
    free(parent);
    return count;
}

/* Build boundary matrix d1: edges → vertices. Matrix is V × E.
   Entry (v, e) = 1 iff vertex v is endpoint of edge e. */
static int *build_d1(const SimplicialComplex *sc)
{
    int V = sc->num_vertices;
    int E = (int)sc->num_edges;
    int *d1 = (int *)calloc(V * E, sizeof(int));
    for (int e = 0; e < E; e++) {
        d1[sc->edges[2*e] * E + e] = 1;
        d1[sc->edges[2*e+1] * E + e] = 1;
    }
    return d1;
}

/* Find edge index by endpoints */
static int find_edge(const SimplicialComplex *sc, int a, int b)
{
    for (size_t e = 0; e < sc->num_edges; e++) {
        int ea = sc->edges[2*e], eb = sc->edges[2*e+1];
        if ((ea == a && eb == b) || (ea == b && eb == a)) return (int)e;
    }
    return -1;
}

/* Build boundary matrix d2: triangles → edges. Matrix is E × F.
   Entry (e, t) = 1 iff edge e is in boundary of triangle t. */
static int *build_d2(const SimplicialComplex *sc)
{
    int E = (int)sc->num_edges;
    int F = (int)sc->num_triangles;
    int *d2 = (int *)calloc(E * F, sizeof(int));
    for (int t = 0; t < F; t++) {
        int a = sc->triangles[t].v[0];
        int b = sc->triangles[t].v[1];
        int c = sc->triangles[t].v[2];
        int e1 = find_edge(sc, a, b);
        int e2 = find_edge(sc, a, c);
        int e3 = find_edge(sc, b, c);
        if (e1 >= 0) d2[e1 * F + t] = 1;
        if (e2 >= 0) d2[e2 * F + t] = 1;
        if (e3 >= 0) d2[e3 * F + t] = 1;
    }
    return d2;
}

/* Build boundary matrix d3: tetrahedra → triangles. Matrix is F × T.
   Entry (t, tet) = 1 iff triangle t is in boundary of tetrahedron tet. */
static int *build_d3(const SimplicialComplex *sc)
{
    int F = (int)sc->num_triangles;
    int T = (int)sc->num_tetrahedra;
    int *d3 = (int *)calloc(F * T, sizeof(int));
    for (int tet = 0; tet < T; tet++) {
        int v0 = sc->tetrahedra[tet].v[0];
        int v1 = sc->tetrahedra[tet].v[1];
        int v2 = sc->tetrahedra[tet].v[2];
        int v3 = sc->tetrahedra[tet].v[3];
        /* Boundary: 4 faces */
        int faces[4][3] = {{v1,v2,v3},{v0,v2,v3},{v0,v1,v3},{v0,v1,v2}};
        for (int f = 0; f < 4; f++) {
            /* Find matching triangle */
            for (int t = 0; t < F; t++) {
                int *tv = sc->triangles[t].v;
                int match = 1;
                int used[3] = {0,0,0};
                for (int i = 0; i < 3; i++) {
                    int found = 0;
                    for (int j = 0; j < 3; j++) {
                        if (!used[j] && tv[j] == faces[f][i]) { used[j]=1; found=1; break; }
                    }
                    if (!found) { match=0; break; }
                }
                if (match) { d3[t * T + tet] = 1; break; }
            }
        }
    }
    return d3;
}

BettiNumbers compute_betti_numbers(const SimplicialComplex *sc)
{
    BettiNumbers bn = {0, 0, 0};
    if (!sc || sc->num_vertices == 0) return bn;

    int V = sc->num_vertices;
    int E = (int)sc->num_edges;
    int F = (int)sc->num_triangles;
    int T = (int)sc->num_tetrahedra;

    /* B0 via union-find */
    bn.b0 = count_connected_components(V, sc->edges, sc->num_edges);

    /* Compute ranks of boundary operators over Z/2 */
    int rank_d1 = 0, rank_d2 = 0, rank_d3 = 0;

    if (E > 0) {
        int *d1 = build_d1(sc);
        rank_d1 = mod2_rank(V, E, d1);
        free(d1);
    }

    if (F > 0 && E > 0) {
        int *d2 = build_d2(sc);
        rank_d2 = mod2_rank(E, F, d2);
        free(d2);
    }

    if (T > 0 && F > 0) {
        int *d3 = build_d3(sc);
        rank_d3 = mod2_rank(F, T, d3);
        free(d3);
    }

    /* B_k = dim(C_k) - rank(d_k) - rank(d_{k+1}) */
    /* B_0 = V - rank(d1) - 0 = V - rank(d1) -- but this should equal components */
    /* Actually with d_0=0: B_0 = dim(ker d_0) - dim(im d_1) = V - rank(d1) */
    /* This equals components count, so it's a consistency check */

    bn.b0 = V - rank_d1; /* should match union-find result */
    bn.b1 = E - rank_d1 - rank_d2;
    bn.b2 = F - rank_d2 - rank_d3;

    if (bn.b1 < 0) bn.b1 = 0;
    if (bn.b2 < 0) bn.b2 = 0;

    /* Use union-find B0 as authoritative (more reliable) */
    bn.b0 = count_connected_components(V, sc->edges, sc->num_edges);

    return bn;
}

/* ---- Persistence via filtration ---- */

PersistenceDiagram persistence_via_filtration(const SimplicialComplex *sc,
                                              const double *vertex_values,
                                              size_t num_values)
{
    PersistenceDiagram pd;
    pd.pairs = NULL;
    pd.count = 0;

    if (!sc || !vertex_values || num_values == 0) return pd;

    /* Create filtration: sort vertices by value */
    idx_val *iv = (idx_val *)malloc(num_values * sizeof(idx_val));
    for (size_t i = 0; i < num_values; i++) { iv[i].idx = (int)i; iv[i].val = vertex_values[i]; }
    qsort(iv, num_values, sizeof(idx_val), cmp_idx_val);
    int *order = (int *)malloc(num_values * sizeof(int));
    for (size_t i = 0; i < num_values; i++) order[i] = iv[i].idx;
    free(iv);

    /* Allocate max possible pairs */
    size_t max_pairs = num_values + sc->num_edges;
    pd.pairs = (PersistencePair *)malloc(max_pairs * sizeof(PersistencePair));

    /* Union-find for component tracking */
    int *parent = (int *)malloc(sc->num_vertices * sizeof(int));
    uf_init(parent, sc->num_vertices);

    /* birth value of each component root */
    double *birth_val = (double *)malloc(sc->num_vertices * sizeof(double));
    for (int i = 0; i < sc->num_vertices; i++) {
        birth_val[i] = (i < (int)num_values) ? vertex_values[i] : 0.0;
    }

    /* Process vertices in filtration order, then edges */
    /* Sort edges by max birth time of endpoints */
    edge_filt *ef = (edge_filt *)malloc(sc->num_edges * sizeof(edge_filt));
    for (size_t e = 0; e < sc->num_edges; e++) {
        int a = sc->edges[2*e], b = sc->edges[2*e+1];
        double fa = (a < (int)num_values) ? vertex_values[a] : 0.0;
        double fb = (b < (int)num_values) ? vertex_values[b] : 0.0;
        ef[e].idx = e;
        ef[e].filt = fa > fb ? fa : fb;
    }
    qsort(ef, sc->num_edges, sizeof(edge_filt), cmp_edge_filt);

    /* Track which vertices have appeared */
    int *appeared = (int *)calloc(sc->num_vertices, sizeof(int));

    int vidx = 0;
    int eidx = 0;

    while (vidx < (int)num_values || eidx < (int)sc->num_edges) {
        double v_time = (vidx < (int)num_values) ? vertex_values[order[vidx]] : DBL_MAX;
        double e_time = (eidx < (int)sc->num_edges) ? ef[eidx].filt : DBL_MAX;

        if (v_time <= e_time && vidx < (int)num_values) {
            int v = order[vidx];
            appeared[v] = 1;
            vidx++;
        } else if (eidx < (int)sc->num_edges) {
            size_t ei = ef[eidx].idx;
            int a = sc->edges[2*ei], b = sc->edges[2*ei+1];
            if (appeared[a] && appeared[b]) {
                int ra = uf_find(parent, a);
                int rb = uf_find(parent, b);
                if (ra != rb) {
                    /* Merge: kill the component born later */
                    double ba = birth_val[ra];
                    double bb = birth_val[rb];
                    int killed, survivor;
                    if (ba > bb) { killed = ra; survivor = rb; }
                    else { killed = rb; survivor = ra; }
                    pd.pairs[pd.count].birth = birth_val[killed];
                    pd.pairs[pd.count].death = e_time;
                    pd.pairs[pd.count].dimension = 0;
                    pd.count++;
                    uf_union(parent, a, b);
                    /* Update survivor's birth */
                    int root = uf_find(parent, a);
                    birth_val[root] = birth_val[survivor];
                }
            }
            eidx++;
        }
    }

    /* Remaining alive components: death = infinity */
    int *seen = (int *)calloc(sc->num_vertices, sizeof(int));
    for (int i = 0; i < sc->num_vertices; i++) {
        if (!appeared[i]) continue;
        int root = uf_find(parent, i);
        if (!seen[root]) {
            seen[root] = 1;
            pd.pairs[pd.count].birth = birth_val[root];
            pd.pairs[pd.count].death = INFINITY;
            pd.pairs[pd.count].dimension = 0;
            pd.count++;
        }
    }

    free(order);
    free(parent);
    free(birth_val);
    free(ef);
    free(appeared);
    free(seen);
    return pd;
}
