#include "m_pd.h"
#include <math.h>

// A rossler chaotic attractor
//
// from ZetaCarinaeModules
// https://github.com/mhampton/ZetaCarinaeModules/blob/master/src/RosslerRustler.cpp

#define br_minimum(x, y) (y < x ? y : x)
#define br_maximum(x, y) (x < y ? y : x)
#define br_clamp(x, minVal, maxVal) (br_minimum(br_maximum(x, minVal), maxVal))

static const float FREQ_C4 = 261.6256f;

static t_class* ross_class;

typedef struct _ross_tilde {
    t_object x_obj;

    float pitch;
    float gain;
    float mix;

    float a_param;
    float b_param;
    float c_param;

    float xout;
    float yout;
    float zout;

    float sampletime;

    t_outlet* x_outlet;
} t_ross;

// --- interface

static void ross_a(t_ross* x, float a)
{
    x->a_param = br_clamp(a, 0, 1);
}

static void ross_b(t_ross* x, float b)
{
    x->b_param = br_clamp(b, 0, 1);
}

static void ross_c(t_ross* x, float c)
{
    x->c_param = br_clamp(c, 0, 30);
}

static void ross_pitch(t_ross* x, float pitch)
{
    x->pitch = br_clamp(pitch, -10, 10);
}

static void ross_mix(t_ross* x, float mix)
{
    x->mix = br_clamp(mix, 0.f, 1.f);
}

static void ross_gain(t_ross* x, float gain)
{
    x->gain = br_clamp(gain, 0.f, 10.f);
}
// --- DSP

static void rossler_slope(float x, float y, float z, float a, float b, float c, float pert, float* output)
{
    output[0] = -y - z;
    output[1] = x + a * y + pert;
    output[2] = b + z * (x - c);
}

static t_int* ross_perform(t_int* w)
{
    t_ross* x = (t_ross*)(w[1]);
    int frames = w[2];
    t_sample* extin = (t_sample*)w[3];
    t_sample* out = (t_sample*)w[4];

    float gain = x->gain;
    float mix = x->mix;

    float A = x->a_param;
    float B = x->b_param;
    float C = x->c_param;

    float pitch = FREQ_C4 * powf(2.f, x->pitch) * 6.2831853f;
    float dt = x->sampletime * pitch / 2.0f;

    for (int i = 0; i < frames; i++) {
        float ext = extin[i];
        float k[3];
        float k2[3];

        rossler_slope(x->xout, x->yout, x->zout, A, B, C, ext * gain, k);
        rossler_slope(x->xout + k[0] * dt, x->yout + k[1] * dt, x->zout + k[2] * dt, A, B, C, ext * gain, k2);

        x->xout += (k[0] + k2[0]) * dt;
        x->yout += (k[1] + k2[1]) * dt;
        x->zout += (k[2] + k2[2]) * dt;

        x->xout = br_clamp(x->xout, -20.f, 20.f);
        x->yout = br_clamp(x->yout, -20.f, 20.f);
        x->zout = br_clamp(x->zout, -20.f, 20.f);

        out[i] = x->xout / 3.0f * (1 - mix) + mix * ext;
    }

    return (w + 5);
}

static void ross_dsp(t_ross* x, t_signal** sp)
{
    dsp_add(ross_perform, 4, x, sp[0]->s_n, sp[0]->s_vec, sp[1]->s_vec);
}

// --- init

static void ross_reset(t_ross* x)
{
    x->xout = 0.f;
    x->yout = 5.f;
    x->zout = 0.f;

    x->pitch = 0;
    x->gain = 0;
    x->mix = 0.5;

    x->a_param = 0.2;
    x->b_param = 0.2;
    x->c_param = 5.7;

    x->sampletime = 1.f / sys_getsr();
}

static void* ross_new()
{
    t_ross* x = (t_ross*)pd_new(ross_class);
    ross_reset(x);

    x->x_outlet = outlet_new(&x->x_obj, &s_signal);

    return (void*)x;
}

static void* ross_free(t_ross* x)
{
    outlet_free(x->x_outlet);
    return (void*)x;
}

void ross_tilde_setup()
{
    ross_class = class_new(gensym("ross~"), (t_newmethod)ross_new, (t_method)ross_free,
        sizeof(t_ross), CLASS_DEFAULT, 0);

    class_addmethod(ross_class, nullfn, gensym("signal"), 0);
    class_addmethod(ross_class, (t_method)ross_dsp, gensym("dsp"), A_CANT, 0);

    class_addmethod(ross_class, (t_method)ross_reset, gensym("reset"), 0);
    class_addmethod(ross_class, (t_method)ross_a, gensym("a"), A_FLOAT, 0);
    class_addmethod(ross_class, (t_method)ross_b, gensym("b"), A_FLOAT, 0);
    class_addmethod(ross_class, (t_method)ross_c, gensym("c"), A_FLOAT, 0);
    class_addmethod(ross_class, (t_method)ross_pitch, gensym("pitch"), A_FLOAT, 0);
    class_addmethod(ross_class, (t_method)ross_mix, gensym("mix"), A_FLOAT, 0);
    class_addmethod(ross_class, (t_method)ross_gain, gensym("gain"), A_FLOAT, 0);
}
