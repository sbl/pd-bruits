#include "m_pd.h"
#include <math.h>
#include <stdbool.h>
#include <time.h>

#include "bruits.h"

#include "mt19937ar/mt19937ar.h"

/**
 * A gendy algorithm after Xenakis
 *
 * See Xenakis, Hoffmann, Lincoln and Serra for literature.
 */

#define MAX_CONTROL_POINTS 128
#define RAMP_TIME 0.02

static t_class* gendy_class;

// gendy definition

typedef enum gendy_distro {
    gendy_uniform = 0,
    gendy_cauchy,
    gendy_lognormal,
    gendy_chisquared,
    gendy_exponential,
    gendy_extreme,
} gendy_distro;

typedef struct _gendy {
    t_object x_obj;

    uint8_t knum;
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
    uint8_t index;
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
    x->minfreq = 220;
    x->maxfreq = 440;
    x->ampdist = gendy_uniform;
    x->ampparam = 0.5;
    x->ampscale = 0.5;
    x->durdist = gendy_uniform;
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

    for (int i = 0; i < MAX_CONTROL_POINTS; i++) {
        x->ampstep1[i] = 2 * genrand_real1() - 1;
        x->ampstep2[i] = 2 * genrand_real1() - 1;

        x->durstep1[i] = 2 * genrand_real1() - 1;
        x->durstep2[i] = genrand_real1();
    }
}

static inline double mirror(double input, double lower, double upper)
{

    if ((input >= lower) && (input <= upper))
        return input;

    double fold_range = 2.0 * fabs((double)(lower - upper));
    return fabs(remainder(input - lower, fold_range)) + lower;
}

static double gendy_distribution(gendy_distro d, double param)
{
    (void)param;
    double z = 0;

    z = genrand_real1();
    // TODO: implement
    switch (d) {
    case gendy_uniform:
        z = genrand_real1();
        break;
    case gendy_cauchy:
        /*z = rng_.cauchy(1, param);*/
        break;
    case gendy_lognormal:
        /*z = rng_.lognormal(0.5, param);*/
        break;
    case gendy_chisquared:
        /*z = rng_.chisquared(param);*/
        break;
    case gendy_exponential:
        /*z = rng_.exponential(param);*/
        break;
    case gendy_extreme:
        /*z = rng_.extreme(2 * param - 1, 1);*/
        break;
    default:
        break;
    }
    return 2 * z - 1.0;
}

static void gendy_debug(t_gendy* x)
{
    post("knum %d", x->knum);
    post("ampdist %d", x->ampdist);
    post("durdist %d", x->durdist);
    post("minfreq %f", x->minfreq);
    post("maxfreq %f", x->maxfreq);
    post("ampscale %f", x->ampscale);
    post("ampparam %f", x->ampparam);
    post("durscale %f", x->durscale);
    post("durparam %f", x->durparam);
}

static void gendy_knum(t_gendy* x, float knum)
{
    uint8_t k = (uint8_t)floorf(knum);
    x->knum = br_clamp(k, 1L, MAX_CONTROL_POINTS);
}

static void gendy_minfreq(t_gendy* x, float minfreq)
{
    x->minfreq = br_clamp(minfreq, 0.000001, 22000);
}

static void gendy_maxfreq(t_gendy* x, float maxfreq)
{
    x->maxfreq = br_clamp(maxfreq, 0.000001, 22000);
}

static void gendy_ampdist(t_gendy* x, float ampdist)
{
    gendy_distro dist = (gendy_distro)floorf(ampdist);
    x->ampdist = dist;
}

static void gendy_ampparam(t_gendy* x, float ampparam)
{
    x->ampparam = br_clamp(ampparam, 0., 1.);
}

static void gendy_ampscale(t_gendy* x, float ampscale)
{
    x->ampscale = br_clamp(ampscale, 0., 1.);
}

static void gendy_durdist(t_gendy* x, float durdist)
{
    gendy_distro dist = (gendy_distro)floorf(durdist);
    x->durdist = dist;
}

static void gendy_durparam(t_gendy* x, float durparam)
{
    x->durparam = br_clamp(durparam, 0., 1.);
}

static void gendy_durscale(t_gendy* x, float durscale)
{
    x->durscale = br_clamp(durscale, 0., 1.);
}

// --- DSP

static t_int* gendy_perform(t_int* w)
{
    t_gendy* x = (t_gendy*)(w[1]);
    t_float* out = (t_float*)(w[2]);
    int frames = (int)(w[3]);

    double rate = x->dur;
    double phase = x->phase;
    double amp = x->amp;
    double nextamp = x->nextamp;
    double speed = x->speed;

    for (int i = 0; i < frames; ++i) {
        double ampscale = x->ampscale;
        double durscale = x->durscale;
        int knum = x->knum;

        double z = 0.0;

        if (phase >= 1) {
            phase -= 1;

            uint8_t index = x->index;
            index = (index + 1) % knum;
            x->index = index;

            // amp
            amp = nextamp;

            gendy_distro ampdist = x->ampdist;
            gendy_distro ampstep = x->ampstep1[index] + gendy_distribution(ampdist, x->ampparam);
            ampstep = mirror(ampstep, -1.0, 1.0);
            x->ampstep1[index] = ampstep;

            nextamp = x->ampstep2[index] + (ampscale * ampstep);
            nextamp = mirror(nextamp, -1.0, 1.0);
            x->ampstep2[index] = nextamp;

            // dur
            gendy_distro durdist = x->durdist;
            double durstep = x->durstep1[index] + gendy_distribution(durdist, x->durparam);
            durstep = mirror(durstep, -1.0, 1.0);
            x->durstep1[index] = durstep;

            rate = x->durstep2[index] + (durscale * durstep);
            rate = mirror(rate, 0.0, 1.0);
            x->durstep2[index] = rate;

            speed = (x->minfreq + ((x->maxfreq - x->minfreq) * rate)) * x->isamplerate;
            speed *= knum;
        }

        z = ((1.0 - phase) * amp) + (phase * nextamp);
        phase += speed;

        out[i] = z;
    }

    x->phase = phase;
    x->amp = amp;
    x->nextamp = nextamp;
    x->speed = speed;
    x->dur = rate;

    return (w + 4);
}

static void gendy_dsp(t_gendy* x, t_signal** sp)
{
    x->isamplerate = 1 / sys_getsr();
    dsp_add(gendy_perform, 3, x, sp[0]->s_vec, (t_int)sp[0]->s_n);
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

    class_addmethod(gendy_class, (t_method)gendy_debug, gensym("debug"), 0);
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
