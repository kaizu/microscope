#ifndef __MICROSCOPE__INTEGRATION
#define __MICROSCOPE__INTEGRATION

#include <stdexcept>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>

#ifdef CUBATURE
#include "cubature.h"
#endif /* CUBATURE */

namespace microscope
{

double integrate1d_simpson(
    double (* function)(double const, void*), void * params,
    double const xmin, double const xmax, const unsigned int N = 10)
{
    const double delta((xmax - xmin) / N);

    double result;

    result = (* function)(xmin, params);
    result += (* function)(xmin + delta * N, params);

    for (unsigned int i(1); i < N; i += 2)
    {
        result += 4.0 * (* function)(xmin + delta * i, params);
    }

    for(unsigned int i(2); i < N; i += 2)
    {
        result += 2.0 * (* function)(xmin + delta * i, params);
    }

    return result * delta / 3.0;
}

double integrate1d_gsl_qng(
    double (* function)(double const, void*), void * params,
    double const xmin, double const xmax)
{
    gsl_function F;
    F.function = function;
    F.params = params;

    double result, abserr;
    size_t neval;
    const double epsabs(1e-8);
    const double epsrel(1e-8);

    gsl_integration_qng(&F, xmin, xmax, epsabs, epsrel, &result, &abserr, &neval);
    return result;
}

double integrate1d_gsl_qags(
    double (* function)(double const, void*), void * params,
    double const xmin, double const xmax)
{
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    gsl_function F;
    F.function = function;
    F.params = params;

    double result, abserr;
    const double epsabs(1e-8);
    const double epsrel(1e-8);

    gsl_integration_qags(&F, xmin, xmax, epsabs, epsrel, 1000, w, &result, &abserr);
    gsl_integration_workspace_free(w);
    return result;
}

inline double integrate1d(
    double (* function)(double const, void*), void * params,
    double const xmin, double const xmax)
{
    return integrate1d_gsl_qags(function, params, xmin, xmax);
    // return integrate1d_gsl_qng(function, params, xmin, xmax);
}

struct Fxy_params
{
    double x;
    double (* function) (double const x, double const y, void * params);
    void * params;
};

double Fxy(double const y, void *params)
{
    struct Fxy_params *p = (struct Fxy_params *) params;
    return p->function(p->x, y, p->params);
}

struct Fx_params
{
    double ymin, ymax;
    double (* function) (double const x, double const y, void * params);
    void* params;
};

double Fx_gsl_qags(double const x, void *params)
{
    struct Fx_params *p = (struct Fx_params *) params;
    struct Fxy_params newp = {x, p->function, p->params};
    return integrate1d_gsl_qags(&Fxy, &newp, p->ymin, p->ymax);
}

double integrate2d_gsl_qags(
    double (* function)(double const, double const, void*), void * params,
    double const xmin, double const xmax, double const ymin, double const ymax)
{
    struct Fx_params p = {ymin, ymax, function, params};
    return integrate1d_gsl_qags(&Fx_gsl_qags, &p, xmin, xmax);
}

double Fx_simpson(double const x, void *params)
{
    struct Fx_params *p = (struct Fx_params *) params;
    struct Fxy_params newp = {x, p->function, p->params};
    return integrate1d_simpson(&Fxy, &newp, p->ymin, p->ymax);
}

double integrate2d_simpson(
    double (* function)(double const, double const, void*), void * params,
    double const xmin, double const xmax, double const ymin, double const ymax)
{
    struct Fx_params p = {ymin, ymax, function, params};
    return integrate1d_simpson(&Fx_simpson, &p, xmin, xmax);
}

#ifdef CUBATURE

struct hcubature_params
{
    double (* function)(double const, double const, void*);
    void *params;
};

int F_hcubature(
    unsigned const ndim, const double *x, void *fdata, unsigned const fdim, double *fval)
{
    struct hcubature_params *p = (struct hcubature_params *) fdata;
    fval[0] = p->function(x[0], x[1], p->params);
    return 0; // success
}

double integrate2d_hcubature(
    double (* function)(double const, double const, void*), void * params,
    double const x0min, double const x0max, double const x1min, double const x1max)
{
    struct hcubature_params p = {function, params};
    double xmin[2] = {x0min, x1min}, xmax[2] = {x0max, x1max}, val, err;
    const double epsabs(1e-8);
    const double epsrel(1e-8);
    hcubature(1, F_hcubature, &p, 2, xmin, xmax,
        0, epsabs, epsrel, ERROR_INDIVIDUAL, &val, &err);
    return val;
}

#endif /* CUBATURE */

inline double integrate2d(
    double (* function)(double const, double const, void*), void * params,
    double const xmin, double const xmax, double const ymin, double const ymax)
{
    return integrate2d_gsl_qags(function, params, xmin, xmax, ymin, ymax);
    // return integrate2d_hcubature(function, params, xmin, xmax, ymin, ymax);
}

double find_root(double (* function)(double const, void*), void * params,
                 gsl_root_fsolver* solver, double const low, double const high,
                 double const tol_abs, double const tol_rel)
{
    gsl_function F;
    F.function = function;
    F.params = params;

    double l(low), h(high);
    gsl_root_fsolver_set(solver, const_cast<gsl_function*>(&F), l, h);

    const unsigned int itermax(100);
    unsigned int i(0);
    for (; ; )
    {
        gsl_root_fsolver_iterate(solver);
        l = gsl_root_fsolver_x_lower(solver);
        h = gsl_root_fsolver_x_upper(solver);

        const int status(gsl_root_test_interval(l, h, tol_abs, tol_rel));
        if (status == GSL_CONTINUE)
        {
            if (i >= itermax)
            {
                gsl_root_fsolver_free(solver);
                throw std::runtime_error("failed to converge");
            }
        }
        else
        {
            break;
        }

        ++i;
    }

    const double root(gsl_root_fsolver_root(solver));
    return root;
}

} // microscope

#endif /* __MICROSCOPE__INTEGRATION */
