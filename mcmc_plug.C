/* mcmc_plug.C
 *
 * Run Markov Chain Monte Carlo parameter estimation
 * for Tempo2
 *
 * Paul Demorest, 2010/04
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <signal.h>
#include "tempo2.h"
#include "T2toolkit.h"
#include "TKfit.h"

static int run=1;
void cc(int sig) { run=0; }

void TKsingularValueDecomposition_lsq(double **A, int n, int nf,
        double **v,double *w,double **u);

double psr_chi2(pulsar *p);
double psr_logL(pulsar *p);
double psr_logPr(pulsar *p);
bool accept_step(double cur, double trial);
void store_params(pulsar *p);
void update_params(pulsar *p, int np, int *mp, int *mt,
        int nj, int *map_j);
void multi_gauss_rand(int n, double **cov, double *x);

double **matrix_malloc(int n) {
    double **res = (double **)malloc(sizeof(double *) * n);
    for (unsigned i=0; i<n; i++) 
        res[i] = (double *)malloc(sizeof(double) * n);
    return res;
}

void matrix_free(int n, double **m) {
    for (unsigned i=0; i<n; i++) free(m[i]);
    free(m);
}

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

static long seed;

// Use a flat independent distribution of cos(i) for trial
// steps.  This is meant to be better for very unconstrained
// Shapiro delay fits.
static bool indep_cosi = false;
static int sini_idx = -1;

// Jump scaling factor for multivar gaussian jumps
static double jfac = 1.3;

/* The main function has to be called this, even if no
 * graphics are involved ;)
 */
extern "C" int graphicalInterface(int argc, char *argv[], 
        pulsar *psr, int *npsr) 
{
    char parFile[MAX_PSR][MAX_FILELEN];
    char timFile[MAX_PSR][MAX_FILELEN];
    *npsr = 1; // Necessary, I guess?

    printf("Markov chain Monte Carlo plugin for tempo2\n");
    printf("Author:    Paul Demorest\n");
    printf("Version:   1.0\n");
    printf("Reference: Demorest et al. 2010, Nature 467, 1081\n");

    parFile[0][0]='\0';
    timFile[0][0]='\0';

    // Parse command line.. there has to be a better way to do this..
    int max_its=-1;
    int thin_fac = 1;
    int show_help = 0;
    for (unsigned i=1; i<argc; i++) {
        if (strcmp(argv[i], "-f")==0) {
            strcpy(parFile[0], argv[i+1]);
            strcpy(timFile[0], argv[i+2]);
        } else if (strcmp(argv[i], "-jfac")==0) {
            jfac = atof(argv[i+1]);
        } else if (strcmp(argv[i], "-cosi")==0) {
            indep_cosi = true;
        } else if (strcmp(argv[i], "-nits")==0) {
            max_its = atoi(argv[i+1]);
        } else if (strcmp(argv[i], "-thin")==0) {
            thin_fac = atoi(argv[i+1]);
        } else if (strcmp(argv[i], "-h")==0) {
            show_help = 1;
        }
    }

    if (parFile[0][0]=='\0' || timFile[0][0]=='\0' || show_help) {
        fprintf(stderr, "Please provide .par and .tim files using -f\n");
        fprintf(stderr, "Output goes to (pulsar_name).mcmc\n");
        fprintf(stderr, "Optional arguments:\n");
        fprintf(stderr, "  -nits N   Run for N iterations (default infinity)\n");
        fprintf(stderr, "  -jfac N   Set MCMC jump scaling factor (default 1.3)\n");
        fprintf(stderr, "  -cosi     Use independent cos(i) in each iter\n");
        fprintf(stderr, "  -thin N   Save only every Nth step\n");
        exit(1);
    }

    //strcpy(parFile[0], "pulsar.par");
    //strcpy(timFile[0], "pulsar.tim");

    // Read the stuff
    readParfile(psr, parFile, timFile, *npsr);
    readTimfile(psr, timFile, *npsr);
    preProcess(psr, *npsr, argc, argv);

    // Do an initial fit to get best-fit params to use as a starting
    // place for MCMC, and cov matrix to use for generating trial
    // jumps.
    formBatsAll(psr, *npsr);
    formResiduals(psr, *npsr, 0); // Final arg is "removeMean"
    doFit(psr, *npsr, 0);         // Final arg is "writeModel"
    formBatsAll(psr, *npsr);      // Run these again to get post-fit
    formResiduals(psr, *npsr, 0); // residuals.

    // Get initial chi2, etc.
    double init_chi2 = psr_chi2(psr);
    printf("fit chi2=%.5e\n", psr[0].fitChisq);
    printf("my  chi2=%.5e\n", init_chi2);

    // Compute map of param array to cov matrix
    char par_names[MAX_PARAMS][32];
    int map_p[MAX_PARAMS], map_t[MAX_PARAMS], map_j[MAX_JUMPS];
    int npfit=0, npmodel=0, npjump=0;
    printf("Fit params:\n");
    for (unsigned i=0; i<MAX_PARAMS; i++) {
        parameter *pp = &(psr[0].param[i]);
        for (unsigned j=0; j<pp->aSize; j++) {
            if (pp->fitFlag[j]>0) {
                map_p[npfit] = i;
                map_t[npfit] = j;
                strcpy(par_names[npfit], pp->shortlabel[j]);
                if (strcmp(par_names[npfit], "SINI")==0) sini_idx=npfit;
                printf("  %5s\n", pp->shortlabel[j]);
                npfit++;
            }
        }
    }
    npmodel = npfit; // Number of non-jump params
    for (unsigned i=1; i<=psr[0].nJumps; i++) {
        if (psr[0].fitJump[i]) {
            printf("  %5s %d %s\n", "JUMP", i, psr[0].jumpStr[i]);
            map_j[npjump] = i;
            npfit++;
            npjump++;
        }
    }
    printf("Total %d params (%d model, %d jumps)\n", 
            npfit, npmodel, npjump);
    if (sini_idx>0) printf("sini_idx=%d\n", sini_idx);

    // Check out covariance, make normalized cov matrix
    double **cov = matrix_malloc(npfit);
    double *errs = (double *)malloc(sizeof(double)*npfit);
    for (unsigned i=0; i<npfit; i++) {
        if (i<npmodel) {
            // Timing model param
            parameter *pi = &(psr[0].param[map_p[i]]);
            errs[i] = pi->err[map_t[i]];
        } else {
            // Jump param
            errs[i] = psr[0].jumpValErr[map_j[i-npmodel]];
        }
        for (unsigned j=0; j<npfit; j++) {
            cov[i][j] = psr[0].covar[i+1][j+1];
            cov[i][j] /= sqrt(psr[0].covar[i+1][i+1]*psr[0].covar[j+1][j+1]);
            // Blank out sini if needed
            if (indep_cosi && (i==sini_idx || j==sini_idx)) {
                if (i==j) cov[i][j] = 1.0;
                else cov[i][j] = 0.0;
            }
        }
    }

    // Rescale cov matrix using param errors. This fixes any 
    // unit discrepancies (PB, etc).
    for (unsigned i=0; i<npfit; i++) 
        for (unsigned j=0; j<npfit; j++)
            cov[i][j] *= errs[i]*errs[j];

    // Set up multivariate gaussian random generator
    multi_gauss_rand(npfit, cov, NULL);

    // Store current parameters in psr[1], trial params in psr[0]
    // Who knows if this has been allocated/etc...
    store_params(psr);

    // Set up output file
    // TODO: store initial params and cov matrix in file
    struct mcmc_output_file out;
    char output_fname[256];
    sprintf(output_fname, "%s.mcmc", psr[0].name);
    init_mcmc_output(&out, output_fname, 
            npmodel, map_p, map_t, npjump, map_j, psr);

    // Main iteration loop
    //int max_its=1000000;
    int it=0, naccept=0;
    int recent_it=0, recent_accept=0, it_stat = 10;
    double cur_logL = psr_logL(psr);
    double cur_chi2 = psr_chi2(psr);
    double cur_logPr = psr_logPr(psr);
    double t0 = get_time();
    double recent_t0=t0, recent_t1=0;
    signal(SIGINT, cc);
    while (run && (it<max_its || max_its<=0)) {

        // Generate trial params
        update_params(psr, npmodel, map_p, map_t, npjump, map_j);

        // Enforce any hard limits (sini<1.0 etc)
        bool a = true;
        if (psr[0].param[param_sini].val[0] > 1.0) a = false;
        if (psr[0].param[param_m2].val[0] < 0.0)   a = false;
        if (psr[0].param[param_mtot].val[0] < 0.0) a = false;

        if (a) {

            // Compute new log-likelihood
            //formBatsAll(psr, *npsr); // Necessary?
            updateBatsAll(psr, *npsr); // Necessary?
            formResiduals(psr, *npsr, 0);
            double trial_logL = psr_logL(psr);
            double trial_logPr = psr_logPr(psr);
            double trial_chi2 = psr_chi2(psr);

            // Check whether or not to accept jump
            a = accept_step(cur_logL + cur_logPr, trial_logL + trial_logPr);

            // Update current state if appropriate
            if (a) {
                store_params(psr);
                cur_logL = trial_logL;
                cur_logPr = trial_logPr;
                cur_chi2 = trial_chi2;
            }

        }

        // Output current state
        if ((it%thin_fac)==0) 
            mcmc_output(&out, &psr[1], cur_chi2);

        // Update stuff
        it++; recent_it++;
        if (a) { naccept++; recent_accept++; }

        // Print status
        if (recent_it == it_stat) {
            recent_t1 = get_time();
            double it_rate = (double)recent_it / (recent_t1-recent_t0);
            printf("%d iters, %.1f it/sec, a=%.3f\n", 
                    it, (float)recent_it/(recent_t1-recent_t0), 
                    (float)recent_accept/(float)recent_it);

            // Tune so ~1 line per xx sec outputs
            it_stat = 2.0*it_rate;

            // Clear counters
            recent_it = recent_accept = 0;
            recent_t0 = recent_t1;

            // Flush output
            fits_flush_file(out.fptr, &out.status);
        }

    }
    double t1 = get_time();

    close_mcmc_output(&out);

    // Print some general diagnostics?
    printf("a=%.3f\n", (double)naccept/(double)it);
    printf("%.1f it/s\n", (double)it/(t1-t0));
}

/* Accept the current trial params by moving them from 
 * psr[0] to psr[1].
 */
void store_params(pulsar *p) {
    copyPSR(p, 0, 1);
    for (unsigned ip=0; ip<MAX_PARAMS; ip++)
        copyParam(p[0].param[ip], &(p[1].param[ip]));
}

void update_params(pulsar *p, int np, int *mp, int *mt, int nj, int *mj) {

    // Copy current good params (1) to trial params (0) 
    copyPSR(p, 1, 0);
    for (unsigned ip=0; ip<MAX_PARAMS; ip++)
        copyParam(p[1].param[ip], &(p[0].param[ip]));

    // Generate random step
    static double *x=NULL;
    int n = np + nj;
    if (x==NULL) x = (double *)malloc(sizeof(double)*n);
    multi_gauss_rand(n, NULL, x);

    // Add it to trial params
    //const double fac = 2.0/sqrt((double)n);
    const double fac = jfac/sqrt((double)n);
    for (unsigned i=0; i<np; i++) {
        parameter *pp = &(p->param[mp[i]]);
        pp->val[mt[i]] += (longdouble)(fac * x[i]);
    }
    for (unsigned i=0; i<nj; i++) {
        p->jumpVal[mj[i]] += fac * x[i+np];
    }

    // Independent cos(i) distribution
    if (indep_cosi) {
        //double max_cosi = 0.33;
        double max_cosi = 1.0;
        double cosi = max_cosi*TKranDev(&seed);
        double sini = sqrt(1.0 - cosi*cosi);
        p->param[mp[sini_idx]].val[0] = sini;
    }
    
}

/* Compare log likelihood ratio and decide whether to accept
 * the step.
 */
bool accept_step(double cur, double trial) {
    if (trial > cur) return true;
    double log_pr = trial - cur;
    double log_rr = log(TKranDev(&seed));
    if (log_rr < log_pr) return true;
    return false;
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

/* Call with x==NULL and cov!=NULL to set up, then vice-versa
 * to get random values.
 */
void multi_gauss_rand(int n, double **cov, double *x) {

    static double **A = NULL;
    static double *e = NULL;
    static double *tmp = NULL;
    unsigned i, j, k;

    // Setup step
    if (cov!=NULL) {

        A = matrix_malloc(n);
        e = (double *)malloc(sizeof(double) * n);
        tmp = (double *)malloc(sizeof(double) * n);

        // Scale cov matrix
        for (i=0; i<n; i++) {
            e[i] = sqrt(cov[i][i]);
            for (j=0; j<n; j++) {
                A[i][j] = cov[i][j] / sqrt(cov[i][i]*cov[j][j]);
            }
        }

        // Use SVD to get the matrix sqrt..
        double **v = matrix_malloc(n);
        double **u = matrix_malloc(n);
        double *w = (double *)malloc(sizeof(double) * n);
        TKsingularValueDecomposition_lsq(A, n, n, v, w, u);
        for (i=0; i<n; i++)
            for (j=0; j<n; j++)
                A[i][j] = 0.0;
        for (i=0; i<n; i++)
            for (j=0; j<n; j++)
                for (k=0; k<n; k++)
                    A[i][j] += v[i][k]*sqrt(w[k])*v[j][k];
                    //A[i][j] += u[i][k]*sqrt(w[k])*v[j][k];
                    /* output "u" seems to be messed up.. that's ok since for
                     * our symmetric matrix, u=v so we can just use v
                     * for this calc.
                     */

        // Free temps
        //for (i=0; i<n; i++) 
        //    for (j=0; j<n; j++) 
        //        printf("%+.6e %+.6e %+.6e\n", u[i][j], v[i][j], A[i][j]);
        matrix_free(n,v);
        matrix_free(n,u);
        free(w);

        seed = TKsetSeed();
    }

    // Generate
    if (x!=NULL) {

        for (i=0; i<n; i++) {
            tmp[i] = TKgaussDev(&seed);
            x[i] = 0.0;
        }

        for (i=0; i<n; i++)
            for (j=0; j<n; j++)
                x[i] += A[i][j]*tmp[j];

        for (i=0; i<n; i++)
            x[i] *= e[i];
    }
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
