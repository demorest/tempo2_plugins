/* resid_plug.C
 *
 * Outputs tempo2 residuals and other useful info
 * in a simple ASCII format.  By default output 
 * goes to resid2.dat.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"


extern "C" int tempoOutput(int argc,char *argv[],pulsar *psr,int npsr) 
{  
  int i,j,npol=0;
  double d,err,det;
  int k;

  const char fname[256] = "resid2.dat";
  FILE *f = fopen(fname, "w");
  if (f==NULL) { 
      fprintf(stderr, "Error opening %s\n", fname); 
      return 0;
  }

  // Print residuals in "print_resid" format.
  // NOTE: I'm not sure what 'torb' actually is, but 
  // it's definitely not 'orbital phase' as was
  // computed by tempo.
  for (int i=0; i<psr[0].nobs; i++) {
    if (psr[0].obsn[i].deleted==0) {
      fprintf(f, "%15.9Lf %9.4f %+.8Le %6.3e %.8Lf\n",
          psr[0].obsn[i].bat,
          psr[0].obsn[i].freqSSB*1e-6,
          psr[0].obsn[i].residual*1e6,
          psr[0].obsn[i].toaErr,
          psr[0].obsn[i].torb);
    }
  }

  fclose(f);

}

