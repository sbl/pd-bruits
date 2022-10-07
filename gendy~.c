#include "m_pd.h"
#include <math.h>
#include <stdbool.h>
#include <time.h>

#include "bruits.h"

#include "mt19937ar/mt19937ar.h"

#define MAX_CONTROL_POINTS 64
#define RAMP_TIME 0.02

static t_class* gendy_class;

// gendy definition

typedef enum gendy_distro {
    uniform = 0,
    cauchy,
    lognormal,
    chisquared,
    exponential,
    extreme,
} gendy_distro;

typedef struct _gendy {
    t_object x_obj;

    uint8_t knum;
    bool active;
    double minfreq;
    double maxfreq;
    gendy_distro ampdist;
    double ampparam;
    double ampscale;
    gendy_distro durdist;
    double durparam;
    double durscale;

    // internal

    double phase;
    unsigned short index;
    double amp;
    double nextamp;
    double dur;
    double speed;

    double ampstep1[MAX_CONTROL_POINTS];
    double ampstep2[MAX_CONTROL_POINTS];

    double durstep1[MAX_CONTROL_POINTS];
    double durstep2[MAX_CONTROL_POINTS];

    double isamplerate;
} t_gendy;

static void gendy_init(t_gendy* x)
{
    x->knum = 12;
    x->active = true;
    x->minfreq = 220;
    x->maxfreq = 440;
    x->ampdist = uniform;
    x->ampparam = 0.5;
    x->ampscale = 0.5;
    x->durdist = uniform;
    x->durparam = 0.5;
    x->durscale = 0.5;

    // internal
    x->phase = 1;
    x->index = 0;
    x->amp = 0;
    x->nextamp = 0;
    x->dur = 1.0;
    x->speed = 1.0;

    x->isamplerate = 1 / sys_getsr();
}

static void gendy_knum(t_gendy* x, uint8_t knum)
{
    x->knum = br_clamp(knum, 1L, MAX_CONTROL_POINTS);
}

static void gendy_active(t_gendy* x, bool active)
{
    x->active = active;
}

static void gendy_minfreq(t_gendy* x, double minfreq)
{
    x->minfreq = minfreq;
}

static void gendy_maxfreq(t_gendy* x, double maxfreq)
{
    x->maxfreq = maxfreq;
}

static void gendy_ampdist(t_gendy* x, gendy_distro ampdist)
{
    x->ampdist = ampdist;
}

static void gendy_ampparam(t_gendy* x, double ampparam)
{
    x->ampparam = br_clamp(ampparam, 0., 1.);
}

static void gendy_ampscale(t_gendy* x, double ampscale)
{
    x->ampscale = br_clamp(ampscale, 0., 1.);
}

static void gendy_durdist(t_gendy* x, gendy_distro durdist)
{
    x->durdist = durdist;
}

static void gendy_durparam(t_gendy* x, double durparam)
{
    x->durparam = br_clamp(durparam, 0., 1.);
}

static void gendy_durscale(t_gendy* x, double durscale)
{
    x->durscale = br_clamp(durscale, 0., 1.);
}

// --- helpers

inline static double mirror(double input, double lower, double upper)
{

    if ((input >= lower) && (input <= upper))
        return input;

    double fold_range = 2.0 * fabs((double)(lower - upper));
    return fabs(remainder(input - lower, fold_range)) + lower;
}

// --- DSP

static t_int* gendy_perform(t_int* w)
{
    t_float* out = (t_float*)(w[1]);
    int n = (int)(w[2]);
    while (n--) {
        *out++ = genrand_real1();
    }
    return (w + 3);
}

static void gendy_dsp(t_gendy* x, t_signal** sp)
{
    x->isamplerate = 1 / sys_getsr();
    dsp_add(gendy_perform, 2, sp[0]->s_vec, (t_int)sp[0]->s_n);
}

static void* gendy_new(void)
{
    t_gendy* x = (t_gendy*)pd_new(gendy_class);

    gendy_init(x);
    init_genrand(time(NULL));

    outlet_new(&x->x_obj, gensym("signal"));
    return (x);
}

void gendy_tilde_setup(void)
{
    gendy_class = class_new(gensym("gendy~"), (t_newmethod)gendy_new, 0,
        sizeof(t_gendy), 0, A_DEFFLOAT, 0);

    class_addmethod(gendy_class, (t_method)gendy_active, gensym("active"), A_FLOAT, 0);
    class_addmethod(gendy_class, (t_method)gendy_knum, gensym("knum"), A_FLOAT, 0);
    class_addmethod(gendy_class, (t_method)gendy_minfreq, gensym("minfreq"), A_FLOAT, 0);
    class_addmethod(gendy_class, (t_method)gendy_maxfreq, gensym("maxfreq"), A_FLOAT, 0);
    class_addmethod(gendy_class, (t_method)gendy_ampparam, gensym("ampparam"), A_FLOAT, 0);
    class_addmethod(gendy_class, (t_method)gendy_durparam, gensym("durparam"), A_FLOAT, 0);
    class_addmethod(gendy_class, (t_method)gendy_ampscale, gensym("ampscale"), A_FLOAT, 0);
    class_addmethod(gendy_class, (t_method)gendy_durscale, gensym("durscale"), A_FLOAT, 0);
    class_addmethod(gendy_class, (t_method)gendy_ampdist, gensym("ampdist"), A_FLOAT, 0);
    class_addmethod(gendy_class, (t_method)gendy_durdist, gensym("durdist"), A_FLOAT, 0);

    class_addmethod(gendy_class, (t_method)gendy_dsp, gensym("dsp"), 0);
}
