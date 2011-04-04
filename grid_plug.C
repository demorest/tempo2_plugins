/* grid_plug.C
 *
 * Do a simple M2-sini Chi2 grid
 *
 * Paul Demorest, 2011/03
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <signal.h>
#include "tempo2.h"

static int run=1;
void cc(int sig) { run=0; }

double psr_chi2(pulsar *p);
void store_params(pulsar *p);
void reset_params(pulsar *p);

double get_time() {
    struct timeval tt;
    gettimeofday(&tt, NULL);
    return (double)tt.tv_sec + 1e-6*(double)tt.tv_usec;
}

#include <fitsio.h>
struct mcmc_output_file {
    fitsfile *fptr;
    int nparam;
    int nmodel;
    int njump;
    int *map_p;
    int *map_t;
    int *map_j;
    int status;
};
void init_mcmc_output(struct mcmc_output_file *f, const char *fname, 
        int npar, int *map_p, int *map_t, 
        int njump, int *map_j,
        pulsar *p);
void mcmc_output(struct mcmc_output_file *f, pulsar *p, double cur_chi2);
void close_mcmc_output(struct mcmc_output_file *f);

/* The main function has to be called this, even if no
 * graphics are involved ;)
 */
extern "C" int graphicalInterface(int argc, char *argv[], 
        pulsar *psr, int *npsr) 
{
    char parFile[MAX_PSR][MAX_FILELEN];
    char timFile[MAX_PSR][MAX_FILELEN];
    *npsr = 1; // Necessary, I guess?

    printf("Chi2 grid plugin for tempo2\n");
    printf("Author:  Paul Demorest\n");
    printf("Version:           0.0\n");

    parFile[0][0]='\0';
    timFile[0][0]='\0';

    // Parse command line.. there has to be a better way to do this..
    double m2_min = 0.0, m2_max = 2.0;
    double sini_min = 0.5, sini_max = 1.0;
    int m2_steps = 128, sini_steps = 128, mtot_steps = 128;
    for (unsigned i=1; i<argc; i++) {
        if (strcmp(argv[i], "-f")==0) {
            strcpy(parFile[0], argv[i+1]);
            strcpy(timFile[0], argv[i+2]);
        } else if (strcmp(argv[i], "-sini_min")==0) {
            sini_min = atof(argv[i+1]);
        } else if (strcmp(argv[i], "-sini_max")==0) {
            sini_max = atof(argv[i+1]);
        } else if (strcmp(argv[i], "-sini_steps")==0) {
            sini_steps = atoi(argv[i+1]);
        } else if (strcmp(argv[i], "-m2_min")==0) {
            m2_min = atof(argv[i+1]);
        } else if (strcmp(argv[i], "-m2_max")==0) {
            m2_max = atof(argv[i+1]);
        } else if (strcmp(argv[i], "-m2_steps")==0) {
            m2_steps = atoi(argv[i+1]);
        }
    }

    if (parFile[0][0]=='\0' || timFile[0][0]=='\0') {
        fprintf(stderr, "Please provide .par and .tim files using -f\n");
#if 1 
        fprintf(stderr, "Optional arguments:\n");
        fprintf(stderr, "  -sini_min N    Minimum sin(i) value (0.5)\n");
        fprintf(stderr, "  -sini_max N    Maximum sin(i) value (1.0)\n");
        fprintf(stderr, "  -sini_steps N  Number of sin(i) steps (128)\n");
        fprintf(stderr, "  -m2_min N      Minimum m2 value (0.0)\n");
        fprintf(stderr, "  -m2_max N      Maximum m2 value (2.0)\n");
        fprintf(stderr, "  -m2_steps N    Number of m2 steps (128)\n");
#endif
        exit(1);
    }

    //strcpy(parFile[0], "pulsar.par");
    //strcpy(timFile[0], "pulsar.tim");

    // Read the stuff
    readParfile(psr, parFile, timFile, *npsr);
    readTimfile(psr, timFile, *npsr);
    preProcess(psr, *npsr, argc, argv);

    // Do an initial fit to get best-fit params to use as a starting
    // place.
    formBatsAll(psr, *npsr);
    formResiduals(psr, *npsr, 0); // Final arg is "removeMean"
    doFit(psr, *npsr, 0);         // Final arg is "writeModel"
    formBatsAll(psr, *npsr);      // Run these again to get post-fit
    formResiduals(psr, *npsr, 0); // residuals.

    // Get initial chi2, etc.
    double init_chi2 = psr_chi2(psr);
    printf("# fit chi2=%.5e\n", psr[0].fitChisq);
    printf("# my  chi2=%.5e\n", init_chi2);

    // step size, etc
    double mtot_min = 0.0, mtot_max = 3.0;
    double delta_m2 = (m2_max - m2_min)/(double)m2_steps;
    double delta_sini = (sini_max - sini_min)/(double)sini_steps;
    double delta_mtot = (mtot_max - mtot_min)/(double)mtot_steps;

    // Store current parameters in psr[1], trial params in psr[0]
    // Who knows if this has been allocated/etc...
    store_params(psr);

    // Set up output file
    // TODO: store initial params and cov matrix in file
    //struct mcmc_output_file out;
    //char output_fname[256];
    //sprintf(output_fname, "%s.mcmc", psr[0].name);
    //init_mcmc_output(&out, output_fname, 
    //        npmodel, map_p, map_t, npjump, map_j, psr);

    // Main iteration loop
    double t0 = get_time();
    double recent_t0=t0, recent_t1=0;
    signal(SIGINT, cc);
    double cur_m2, cur_sini, cur_mtot, cur_chi2;
    for (cur_sini=sini_min; cur_sini<sini_max; cur_sini+=delta_sini) {
    //for (cur_mtot=mtot_min; cur_mtot<mtot_max; cur_mtot+=delta_mtot) {
        for (cur_m2=m2_min; cur_m2<m2_max; cur_m2+=delta_m2) {

            reset_params(psr);

            psr[0].param[param_sini].val[0] = cur_sini;
            psr[0].param[param_sini].paramSet[0] = 1;
            //psr[0].param[param_mtot].val[0] = cur_mtot;
            psr[0].param[param_m2].val[0] = cur_m2;
            psr[0].param[param_m2].paramSet[0] = 1;

            // updateBats or FormBats?
            int ok=1;
            updateBatsAll(psr, *npsr);
            try {
                formResiduals(psr, *npsr, 0);
                doFit(psr, *npsr, 0);
            } catch (int e) {
                if (e==5) { ok=0; }
            }
            if (ok) {

                updateBatsAll(psr, *npsr);
                formResiduals(psr, *npsr, 0);
                cur_chi2 = psr_chi2(psr);

            } else {
                cur_chi2 = 1e6;
            }

            printf("%.9f %.5f %.6e\n", cur_sini, cur_m2, cur_chi2);

            if (run==0) break;
        }
        if (run==0) break;
    }

    //close_mcmc_output(&out);
}

/* Accept the current trial params by moving them from 
 * psr[0] to psr[1].
 */
void store_params(pulsar *p) {
    copyPSR(p, 0, 1);
    for (unsigned ip=0; ip<MAX_PARAMS; ip++)
        copyParam(p[0].param[ip], &(p[1].param[ip]));
}

/* Accept the current trial params by moving them from 
 * psr[0] to psr[1].
 */
void reset_params(pulsar *p) {
    copyPSR(p, 1, 0);
    for (unsigned ip=0; ip<MAX_PARAMS; ip++)
        copyParam(p[1].param[ip], &(p[0].param[ip]));
}

/* Compute the chi2 value of residuals in the struct.
 * This removes a weighted mean.
 */
double psr_chi2(pulsar *p) {
    double sum=0.0, sum2=0.0, wsum=0.0;
    for (unsigned i=0; i<p->nobs; i++) {
        double r = 1e6*p->obsn[i].residual;
        double e = p->obsn[i].toaErr;
        sum2 += r*r/(e*e);
        sum += r/(e*e);
        wsum += 1.0/(e*e);
    }
    return (sum2 - sum*sum/wsum);
}

/* Return the log-likelihood for the current fit, assuming
 * gaussian stats.  "toaErr" values already have efac, equad, 
 * etc incorporated.
 */
double psr_logL(pulsar *p) {
    double r_sum=0.0, r2_sum=0.0, w_sum=0.0, le_sum=0.0;
    for (unsigned i=0; i<p->nobs; i++) {
        double r = 1e6*p->obsn[i].residual;
        double e = p->obsn[i].toaErr;
        r_sum  += r/(e*e);
        r2_sum += r*r/(e*e);
        w_sum  += 1.0/(e*e);
        le_sum += log(e);
    }
    double chi2 = r2_sum - r_sum*r_sum/w_sum;
    return -1.0*le_sum - 0.5*chi2;
}

/* Return the log-PDF of any priors on the pulsar params. */
double psr_logPr(pulsar *p) {

    /* Flat cos(i) prior */
    double sini = p->param[param_sini].val[0];
    double pr = sini / sqrt(1.0 - sini*sini);
    return log(pr);

    /* Uniform prior on all params */
    //return 0.0;
}

void init_mcmc_output(struct mcmc_output_file *f, const char *fname, 
        int nmodel, int *map_p, int *map_t, 
        int njump, int *map_j,
        pulsar *p) {

    f->status = 0;
    int *s = &f->status;
    fits_create_file(&f->fptr, fname, s);
    if (*s) {
        fprintf(stderr, "Error opening %s for output", fname);
        fits_report_error(stderr, *s);
        exit(1);
    }

    f->nparam = nmodel + njump;
    f->nmodel = nmodel;
    f->njump = njump;
    f->map_p = map_p;
    f->map_t = map_t;
    f->map_j = map_j;

    // Main output is a fits table
    fits_create_tbl(f->fptr, BINARY_TBL, 0, 0, NULL, NULL, NULL, "MCMC", s);

    // Create columns
    int ncol=0;
    for (unsigned i=0; i<nmodel; i++) {
        parameter *pp = &(p->param[map_p[i]]);
        ncol++;
        fits_insert_col(f->fptr, ncol, pp->shortlabel[map_t[i]], "D", s);
    }
    for (unsigned i=0; i<njump; i++) {
        ncol++;
        char tmps[256];
        sprintf(tmps, "JUMP_%d", map_j[i]);
        fits_insert_col(f->fptr, ncol, tmps, "D", s);
    }
    ncol++;
    fits_insert_col(f->fptr, ncol, "CHI2", "D", s);
    if (*s) {
        fits_report_error(stderr, *s);
        exit(1);
    }
}

void mcmc_output(struct mcmc_output_file *f, pulsar *p, double chi2) {
    long nrow;
    fits_get_num_rows(f->fptr, &nrow, &f->status);
    fits_insert_rows(f->fptr, nrow, 1, &f->status);
    for (unsigned i=0; i<f->nmodel; i++) {
        parameter *pp = &(p->param[f->map_p[i]]);
        double val = pp->val[f->map_t[i]];
        fits_write_col(f->fptr, TDOUBLE, i+1, nrow+1, 
                1, 1, &val, &f->status);
    }
    for (unsigned i=0; i<f->njump; i++) {
        double val = p->jumpVal[f->map_j[i]];
        fits_write_col(f->fptr, TDOUBLE, i+f->nmodel+1, nrow+1, 
                1, 1, &val, &f->status);
    }
    fits_write_col(f->fptr, TDOUBLE, f->nparam+1, nrow+1, 
            1, 1, &chi2, &f->status);
    if (f->status) 
        fits_report_error(stderr, f->status);
}

void close_mcmc_output(struct mcmc_output_file *f) {
    fits_close_file(f->fptr, &f->status);
}
