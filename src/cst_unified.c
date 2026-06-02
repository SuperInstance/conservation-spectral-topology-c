#include "cst_core.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

bool cst_verify(const CSTSystem *sys)
{
    if (!sys) return false;
    return conservation_check(&sys->quantity, sys->conservation_tolerance);
}

double cst_spectral_health(const CSTSystem *sys)
{
    if (!sys || sys->num_eigenvalues < 2) return 0.0;
    return spectral_gap(sys->eigenvalues, sys->num_eigenvalues);
}

CstPhase cst_phase_detect(const CSTSystem *sys)
{
    if (!sys) return CST_PHASE_CRITICAL;

    bool conserved = cst_verify(sys);
    double gap = cst_spectral_health(sys);

    if (!conserved && gap < sys->critical_gap)
        return CST_PHASE_CRITICAL;
    if (!conserved || gap < sys->warning_gap)
        return CST_PHASE_WARNING;
    return CST_PHASE_NORMAL;
}

double cst_predict(const CSTSystem *sys)
{
    if (!sys) return 0.0;
    double gap = cst_spectral_health(sys);
    /* Distance to phase transition: normalized gap.
       0 = gap is zero, transition imminent.
       1 = gap is very large, far from transition. */
    if (gap <= 0.0) return 0.0;
    return gap / (gap + sys->warning_gap);
}

void cst_system_free(CSTSystem *sys)
{
    if (!sys) return;
    free(sys->eigenvalues);
    sys->eigenvalues = NULL;
    simplicial_complex_free(&sys->topology);
}
