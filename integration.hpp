#ifndef __MICROSCOPE__INTEGRATION
#define __MICROSCOPE__INTEGRATION

#include <gsl/gsl_integration.h>


namespace microscope
{

double integrate1d(
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

double integrate2d(
    double (* function)(double, double, void*), void * params,
    double xmin, double xmax, double ymin, double ymax)
{
    struct Fx_params p = {ymin, ymax, function, params};
    return integrate1d(&Fx, &p, xmin, xmax);
}

} // microscope

#endif /* __MICROSCOPE__INTEGRATION */
