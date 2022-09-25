#include "m_pd.h"
#include <stdlib.h>

// A rossler chaotic attractor
//
// from MCLDChaosUgen
// https://github.com/supercollider/sc3-plugins/blob/main/source/MCLDUGens/MCLDChaosUGens.cpp

#define sc_max(a, b) (((a) > (b)) ? (a) : (b))
#define ONESIXTH 0.1666666666666667

static t_class* ross_class;

typedef struct _ross_tilde {
    t_object x_obj;

    double x0, y0, xn, yn, xnm1, ynm1;
    float counter;
    double z0, zn, znm1, frac;
    double sr;

    // interface
    float freq;
    float c;
    double a, b, h;

    t_inlet* c_in;

    t_outlet* x_out;
    t_outlet* y_out;
    t_outlet* z_out;
} t_ross;

// --- interface

void ross_float(t_ross* x, float f)
{
    x->freq = f;
    post("freq %f", x->freq);
}

// --- DSP

static t_int* ross_perform(t_int* w)
{
    t_ross* x = (t_ross*)(w[1]);
    int frames = w[2];

    t_sample* xout = (t_sample*)w[3];
    t_sample* yout = (t_sample*)w[4];
    t_sample* zout = (t_sample*)w[5];

    float freq = x->freq;
    double a = x->a;
    double b = x->b;
    double c = x->c;
    double h = x->h;
    double x0 = x->x0;
    double y0 = x->y0;
    double z0 = x->z0;

    double xn = x->xn;
    double yn = x->yn;
    double zn = x->zn;
    float counter = x->counter;
    double xnm1 = x->xnm1;
    double ynm1 = x->ynm1;
    double znm1 = x->znm1;
    double frac = x->frac;

    float samplesPerCycle;
    double slope;
    if (freq < x->sr) {
        samplesPerCycle = x->sr / sc_max(freq, 0.001f);
        slope = 1.f / samplesPerCycle;
    } else {
        samplesPerCycle = 1.f;
        slope = 1.f;
    }

    if ((x->x0 != x0) || (x->y0 != y0) || (x->z0 != z0)) {
        xnm1 = xn;
        ynm1 = yn;
        znm1 = zn;
        x->x0 = xn = x0;
        x->y0 = yn = y0;
        x->z0 = zn = z0;
    }

    double dx = xn - xnm1;
    double dy = yn - ynm1;
    double dz = zn - znm1;

    for (int i = 0; i < frames; ++i) {
        if (counter >= samplesPerCycle) {
            counter -= samplesPerCycle;
            frac = 0.f;

            xnm1 = xn;
            ynm1 = yn;
            znm1 = zn;

            double k1x, k2x, k3x, k4x, k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z, kxHalf, kyHalf,
                kzHalf;

            // 4th order Runge-Kutta approximate solution for differential equations
            k1x = -(h * (ynm1 + znm1));
            k1y = h * (xnm1 + a * ynm1);
            k1z = h * (b + znm1 * (xnm1 - c));
            kxHalf = k1x * 0.5;
            kyHalf = k1y * 0.5;
            kzHalf = k1z * 0.5;

            k2x = -(h * (ynm1 + kyHalf + znm1 + kzHalf));
            k2y = h * (xnm1 + kxHalf + a * (ynm1 + kyHalf));
            k2z = h * (b + (znm1 + kzHalf) * (xnm1 + kxHalf - c));
            kxHalf = k2x * 0.5;
            kyHalf = k2y * 0.5;
            kzHalf = k2z * 0.5;

            k3x = -(h * (ynm1 + kyHalf + znm1 + kzHalf));
            k3y = h * (xnm1 + kxHalf + a * (ynm1 + kyHalf));
            k3z = h * (b + (znm1 + kzHalf) * (xnm1 + kxHalf - c));

            k4x = -(h * (ynm1 + k3y + znm1 + k3z));
            k4y = h * (xnm1 + k3x + a * (ynm1 + k3y));
            k4z = h * (b + (znm1 + k3z) * (xnm1 + k3x - c));

            xn = xn + (k1x + 2.0 * (k2x + k3x) + k4x) * ONESIXTH;
            yn = yn + (k1y + 2.0 * (k2y + k3y) + k4y) * ONESIXTH;
            zn = zn + (k1z + 2.0 * (k2z + k3z) + k4z) * ONESIXTH;

            dx = xn - xnm1;
            dy = yn - ynm1;
            dz = zn - znm1;
        }
        counter++;
        xout[i] = (xnm1 + dx * frac) * 0.5f;
        yout[i] = (ynm1 + dy * frac) * 0.5f;
        zout[i] = (znm1 + dz * frac) * 1.0f;

        frac += slope;
    }

    x->xn = xn;
    x->yn = yn;
    x->zn = zn;
    x->counter = counter;
    x->xnm1 = xnm1;
    x->ynm1 = ynm1;
    x->znm1 = znm1;
    x->frac = frac;

    return (w + 6);
}

static void ross_dsp(t_ross* x, t_signal** sp)
{
    x->sr = sys_getsr();
    dsp_add(ross_perform, 5, x, sp[0]->s_n, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec);
}

// --- init

static void* ross_new()
{
    t_ross* x = (t_ross*)pd_new(ross_class);
    x->x0 = x->xn = x->xnm1 = 0.1;
    x->y0 = x->yn = x->ynm1 = 0;
    x->z0 = x->zn = x->znm1 = 0;
    x->counter = 0.f;
    x->frac = 0.f;

    x->freq = 5000;
    x->a = 0.2;
    x->b = 0.2;
    x->c = 5.7;
    x->h = 0.05;

    x->c_in = floatinlet_new(&x->x_obj, &x->c);
    x->x_out = outlet_new(&x->x_obj, &s_signal);
    x->y_out = outlet_new(&x->x_obj, &s_signal);
    x->z_out = outlet_new(&x->x_obj, &s_signal);

    return (void*)x;
}

static void* ross_free(t_ross* x)
{
    inlet_free(x->c_in);

    outlet_free(x->x_out);
    outlet_free(x->y_out);
    outlet_free(x->z_out);
    return (void*)x;
}

void ross_tilde_setup()
{
    ross_class = class_new(gensym("ross~"), (t_newmethod)ross_new, (t_method)ross_free,
        sizeof(t_ross), CLASS_DEFAULT, 0);

    class_addfloat(ross_class, ross_float);
    class_addmethod(ross_class, (t_method)ross_dsp, gensym("dsp"), A_CANT, 0);
}
