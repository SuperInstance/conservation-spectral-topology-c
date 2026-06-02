#include "cst_core.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ---- Sparse matrix utilities ---- */

void sparse_matrix_free(SparseMatrix *m)
{
    if (!m) return;
    free(m->row_ptr);  m->row_ptr  = NULL;
    free(m->col_idx);  m->col_idx  = NULL;
    free(m->values);   m->values   = NULL;
}

/* Build graph Laplacian: L = D - A */
SparseMatrix laplacian_from_complex(const SimplicialComplex *sc)
{
    SparseMatrix L;
    memset(&L, 0, sizeof(L));
    L.nrows = sc->num_vertices;
    L.ncols = sc->num_vertices;

    if (sc->num_vertices == 0) return L;

    int *deg = (int *)calloc(sc->num_vertices, sizeof(int));
    int *nnz_per_row = (int *)calloc(sc->num_vertices, sizeof(int));

    for (size_t e = 0; e < sc->num_edges; e++) {
        int a = sc->edges[2*e], b = sc->edges[2*e+1];
        if (a == b) continue;
        deg[a]++;
        deg[b]++;
        nnz_per_row[a]++;
        nnz_per_row[b]++;
    }

    L.nnz = 0;
    for (int i = 0; i < sc->num_vertices; i++)
        L.nnz += 1 + nnz_per_row[i];

    L.row_ptr = (int *)malloc((sc->num_vertices + 1) * sizeof(int));
    L.col_idx = (int *)malloc(L.nnz * sizeof(int));
    L.values  = (double *)malloc(L.nnz * sizeof(double));

    L.row_ptr[0] = 0;
    for (int i = 0; i < sc->num_vertices; i++)
        L.row_ptr[i+1] = L.row_ptr[i] + 1 + nnz_per_row[i];

    int *pos = (int *)malloc(sc->num_vertices * sizeof(int));
    for (int i = 0; i < sc->num_vertices; i++) {
        pos[i] = L.row_ptr[i];
        L.col_idx[pos[i]] = i;
        L.values[pos[i]] = (double)deg[i];
        pos[i]++;
    }

    for (size_t e = 0; e < sc->num_edges; e++) {
        int a = sc->edges[2*e], b = sc->edges[2*e+1];
        if (a == b) continue;
        L.col_idx[pos[a]] = b;
        L.values[pos[a]] = -1.0;
        pos[a]++;
        L.col_idx[pos[b]] = a;
        L.values[pos[b]] = -1.0;
        pos[b]++;
    }

    free(deg);
    free(nnz_per_row);
    free(pos);
    return L;
}

/* Sparse matrix-vector multiply: y = A * x */
static void spmv(const SparseMatrix *A, const double *x, double *y)
{
    for (int i = 0; i < A->nrows; i++) {
        double sum = 0.0;
        for (int j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
            sum += A->values[j] * x[A->col_idx[j]];
        y[i] = sum;
    }
}

/* ---- Power iteration with proper deflation ---- */
int power_iteration(const SparseMatrix *A, int k,
                    double *eigenvalues, double *eigenvectors,
                    int max_iter, double tol)
{
    if (!A || A->nrows == 0 || k <= 0) return -1;
    int n = A->nrows;

    double *v  = (double *)malloc(n * sizeof(double));
    double *w  = (double *)malloc(n * sizeof(double));

    for (int eig = 0; eig < k; eig++) {
        /* Random initial vector */
        unsigned int seed = (unsigned int)(42 + eig * 1337);
        for (int i = 0; i < n; i++) {
            seed = seed * 1103515245 + 12345;
            v[i] = ((double)((seed >> 16) & 0x7fff) / 32768.0) - 0.5;
        }

        /* Orthogonalize against previous eigenvectors */
        for (int prev = 0; prev < eig; prev++) {
            double dot = 0.0;
            for (int i = 0; i < n; i++)
                dot += v[i] * eigenvectors[prev * n + i];
            for (int i = 0; i < n; i++)
                v[i] -= dot * eigenvectors[prev * n + i];
        }

        /* Normalize */
        double norm = 0.0;
        for (int i = 0; i < n; i++) norm += v[i] * v[i];
        norm = sqrt(norm);
        if (norm < 1e-15) { v[0] = 1.0; norm = 1.0; }
        for (int i = 0; i < n; i++) v[i] /= norm;

        double lambda = 0.0;
        for (int iter = 0; iter < max_iter; iter++) {
            spmv(A, v, w);

            /* Deflate: remove projections onto previous eigenvectors */
            for (int prev = 0; prev < eig; prev++) {
                double dot = 0.0;
                for (int i = 0; i < n; i++)
                    dot += w[i] * eigenvectors[prev * n + i];
                for (int i = 0; i < n; i++)
                    w[i] -= dot * eigenvectors[prev * n + i];
            }

            /* Rayleigh quotient */
            lambda = 0.0;
            for (int i = 0; i < n; i++) lambda += w[i] * v[i];

            /* Normalize w */
            norm = 0.0;
            for (int i = 0; i < n; i++) norm += w[i] * w[i];
            norm = sqrt(norm);
            if (norm < 1e-15) {
                /* Eigenvalue is 0, w collapsed — set v to arbitrary orthonormal */
                for (int i = 0; i < n; i++) w[i] = 0.0;
                /* Find a component orthogonal to previous eigenvectors */
                for (int i = 0; i < n; i++) {
                    w[i] = 1.0;
                    for (int prev = 0; prev < eig; prev++) {
                        double dot = 0.0;
                        for (int j = 0; j < n; j++)
                            dot += w[j] * eigenvectors[prev * n + j];
                        for (int j = 0; j < n; j++)
                            w[j] -= dot * eigenvectors[prev * n + j];
                    }
                    norm = 0.0;
                    for (int j = 0; j < n; j++) norm += w[j] * w[j];
                    if (norm > 1e-10) break;
                    for (int j = 0; j < n; j++) w[j] = 0.0;
                }
                norm = sqrt(norm);
                if (norm < 1e-15) { lambda = 0.0; break; }
            }
            for (int i = 0; i < n; i++) w[i] /= norm;

            /* Convergence check */
            double diff = 0.0, diff2 = 0.0;
            for (int i = 0; i < n; i++) {
                double d = w[i] - v[i]; diff += d * d;
                d = w[i] + v[i]; diff2 += d * d;
            }
            memcpy(v, w, n * sizeof(double));
            if (diff < tol * tol || diff2 < tol * tol) break;
        }

        eigenvalues[eig] = lambda;
        memcpy(&eigenvectors[eig * n], v, n * sizeof(double));
    }

    free(v);
    free(w);
    return 0;
}

/* ---- Spectral gap ---- */
double spectral_gap(const double *eigenvalues, int n)
{
    if (n < 2) return 0.0;

    double *sorted = (double *)malloc(n * sizeof(double));
    memcpy(sorted, eigenvalues, n * sizeof(double));
    for (int i = 0; i < n - 1; i++)
        for (int j = i + 1; j < n; j++)
            if (sorted[j] < sorted[i]) { double t = sorted[i]; sorted[i] = sorted[j]; sorted[j] = t; }

    double gap = sorted[1] - sorted[0];
    free(sorted);
    return gap;
}

/* ---- Cheeger constant from Fiedler vector ---- */
double cheeger_constant(const double *fiedler, int n)
{
    if (n < 2) return 0.0;

    double norm_sq = 0.0;
    for (int i = 0; i < n; i++) norm_sq += fiedler[i] * fiedler[i];
    if (norm_sq < 1e-15) return 0.0;
    double norm = sqrt(norm_sq);

    int pos_count = 0;
    for (int i = 0; i < n; i++)
        if (fiedler[i] >= 0) pos_count++;
    int neg_count = n - pos_count;
    double min_side = (pos_count < neg_count) ? (double)pos_count : (double)neg_count;
    if (min_side == 0) return 0.0;

    return norm / (2.0 * min_side);
}
