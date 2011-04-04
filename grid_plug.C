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
    int m2_steps = 128, sini_steps = 128;
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
        fprintf(stderr, "Optional arguments:\n");
        fprintf(stderr, "  -sini_min N    Minimum sin(i) value (0.5)\n");
        fprintf(stderr, "  -sini_max N    Maximum sin(i) value (1.0)\n");
        fprintf(stderr, "  -sini_steps N  Number of sin(i) steps (128)\n");
        fprintf(stderr, "  -m2_min N      Minimum m2 value (0.0)\n");
        fprintf(stderr, "  -m2_max N      Maximum m2 value (2.0)\n");
        fprintf(stderr, "  -m2_steps N    Number of m2 steps (128)\n");
        exit(1);
    }

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

    // step size, etc
    double delta_m2 = (m2_max - m2_min)/(double)m2_steps;
    double delta_sini = (sini_max - sini_min)/(double)sini_steps;

    // Store current parameters in psr[1], trial params in psr[0]
    store_params(psr);

    // Print header lines
    printf("# chi2_0=%.5e\n", init_chi2);
    printf("# SINI M2 chi2-chi2_0\n");

    // Main iteration loop
    signal(SIGINT, cc);
    double cur_m2, cur_sini, cur_chi2;
    for (cur_sini=sini_min; cur_sini<sini_max; cur_sini+=delta_sini) {
        for (cur_m2=m2_min; cur_m2<m2_max; cur_m2+=delta_m2) {

            reset_params(psr);

            psr[0].param[param_sini].val[0] = cur_sini;
            psr[0].param[param_sini].paramSet[0] = 1;
            psr[0].param[param_m2].val[0] = cur_m2;
            psr[0].param[param_m2].paramSet[0] = 1;

            // updateBats or FormBats?
            int ok=1;
            updateBatsAll(psr, *npsr);
            try {
                formResiduals(psr, *npsr, 0);
                doFit(psr, *npsr, 0);
            } catch (int e) {
                // Note, this is a hack to ignore spots in 
                // M2-sini parameter space where the fit does not
                // converge.  Plain tempo2 exit()'s here, killing
                // the whole process.  I've modified my copy to
                // throw an error instead so we can continue on
                // with the gridding.
                if (e==5) { ok=0; }
            }
            if (ok) {

                updateBatsAll(psr, *npsr);
                formResiduals(psr, *npsr, 0);
                cur_chi2 = psr_chi2(psr);

            } else {
                cur_chi2 = 1e6;
            }

            printf("%.9f %.5f %.10e\n", cur_sini, cur_m2, cur_chi2-init_chi2);

            if (run==0) break;
        }
        if (run==0) break;
    }

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
