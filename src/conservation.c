#include "cst_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* ------------------------------------------------------------------ */
bool conservation_check(const ConservedQuantity *q, double tolerance)
{
    double expected_flux = q->source - q->sink;
    double diff = fabs(q->flux - expected_flux);
    return diff <= tolerance;
}

/* ------------------------------------------------------------------ */
/* Discrete Noether: sum charges over all vertices, simulate one Euler
   step with adjacency-weighted diffusion, check total charge conserved. */
bool discrete_noether_check(const double *adjacency, const double *charge_per_vertex,
                            size_t n, double dt, double tolerance)
{
    if (n == 0) return true;

    double total_before = 0.0;
    for (size_t i = 0; i < n; i++)
        total_before += charge_per_vertex[i];

    /* Diffusion step: dQ_i/dt = sum_j A_ij * (Q_j - Q_i) */
    double *next = (double *)malloc(n * sizeof(double));
    if (!next) return false;
    memcpy(next, charge_per_vertex, n * sizeof(double));

    for (size_t i = 0; i < n; i++) {
        double delta = 0.0;
        for (size_t j = 0; j < n; j++) {
            delta += adjacency[i * n + j] * (charge_per_vertex[j] - charge_per_vertex[i]);
        }
        next[i] += dt * delta;
    }

    double total_after = 0.0;
    for (size_t i = 0; i < n; i++)
        total_after += next[i];

    free(next);
    return fabs(total_after - total_before) <= tolerance;
}
