#ifndef __MICROSCOPE__INTEGRATION
#define __MICROSCOPE__INTEGRATION

#include <gsl/gsl_integration.h>
#include "cubature.h"


namespace microscope
{

double integrate1d_gsl_qng(
    double (* function)(double, void*), void * params,
    double xmin, double xmax)
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
    double (* function)(double, void*), void * params,
    double xmin, double xmax)
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
    double (* function)(double, void*), void * params,
    double xmin, double xmax)
{
    return integrate1d_gsl_qags(function, params, xmin, xmax);
}

struct Fxy_params
{
    double x;
    double (* function) (double x, double y, void * params);
    void * params;
};

double Fxy(double y, void *params)
{
    struct Fxy_params *p = (struct Fxy_params *) params;
    return p->function(p->x, y, p->params);
}

struct Fx_params
{
    double ymin, ymax;
    double (* function) (double x, double y, void * params);
    void* params;
};

double Fx(double x, void *params)
{
    struct Fx_params *p = (struct Fx_params *) params;
    struct Fxy_params newp = {x, p->function, p->params};
    return integrate1d(&Fxy, &newp, p->ymin, p->ymax);
}

double integrate2d_gsl_qags(
    double (* function)(double, double, void*), void * params,
    double xmin, double xmax, double ymin, double ymax)
{
    struct Fx_params p = {ymin, ymax, function, params};
    return integrate1d_gsl_qags(&Fx, &p, xmin, xmax);
}

struct hcubature_params
{
    double (* function)(double, double, void*);
    void *params;
};

int F_hcubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
    struct hcubature_params *p = (struct hcubature_params *) fdata;
    fval[0] = p->function(x[0], x[1], p->params);
    return 0; // success
}

double integrate2d_hcubature(
    double (* function)(double, double, void*), void * params,
    double x0min, double x0max, double x1min, double x1max)
{
    struct hcubature_params p = {function, params};
    double xmin[2] = {x0min, x1min}, xmax[2] = {x0max, x1max}, val, err;
    hcubature(1, F_hcubature, &p, 2, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
    return val;
}

inline double integrate2d(
    double (* function)(double, double, void*), void * params,
    double xmin, double xmax, double ymin, double ymax)
{
    return integrate2d_gsl_qags(function, params, xmin, xmax, ymin, ymax);
}

} // microscope

#endif /* __MICROSCOPE__INTEGRATION */
