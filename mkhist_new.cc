#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define debug 0

#define R_win 5.0
#define bin_width 0.05
#define bin_num (int)(R_win / bin_width * 2.0 + 1.0)

#define deg_rad 0.017453292519943296 /*   pi/180  */
#define rad_deg 57.295779513082322   /*   180/pi  */
#define pi 3.14159265358979323846    /* pi */

int main(int argc, char *argv[]) {

  FILE *fpon, *fpoff;
  double on[bin_num][bin_num];
  double off[bin_num][bin_num];
  int i, j;
  int id, xd, yd;
  double maxon = 0;
  double maxoff = 0;
  double corr = 1.0;

  int nhit, ndetc, nfitc;
  double mjd, npea, zen, azi, xc, yc, npec, zenc, azic, omega, chi2p, chi2c,
      compactness, pincness, theta, phi;

  memset(&on[0][0], '\0', sizeof on);
  memset(&off[0][0], '\0', sizeof off);

  while (scanf("%lf %d %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf "
               "%lf %d %lf %lf",
               &mjd, &nhit, &npea, &zen, &azi, &xc, &yc, &ndetc, &nfitc, &npec,
               &zenc, &azic, &omega, &chi2p, &chi2c, &compactness, &pincness,
               &id, &theta, &phi) > 0) {
    //    if(zen>=50) continue;
    if (nfitc > 100 && zenc < 50. && zenc > 0. ) {
  //     double g = gmst( mjd ) + tibet_lo;
  // double lst = ( g - floor( g / 360.0 ) * 360.0 );
if (mjd-floor(mjd)>0.32)
 continue;  
      xd = floor((theta * sin(phi * deg_rad) + R_win) / bin_width + 0.5);
      yd = floor((theta * cos(phi * deg_rad) + R_win) / bin_width + 0.5);
      // yd = floor( ( r * sin( theta * deg_rad ) + 5.0 ) * DIFm +0.5 );
      // xd = floor( ( r * cos( theta * deg_rad ) + 5.0 ) * DIFm +0.5 );
      if (id == 0) {
        maxon++;
        on[xd][yd] += corr;
        //    fprintf(stderr,"any4:%d\n",any4);
      }
      if (id != 0) {
        maxoff++;
        off[xd][yd] += corr;
      }
    }
  }
  fpon = stdout;
  for (i = 0; i < bin_num; i++) {
    for (j = 0; j < bin_num; j++) {
      // if (on[ i][j ]>0 || off[ i][j ]>0)
      fprintf(fpon, "%d %d %lf %lf\n", i, j, on[i][j], off[i][j]);
      if (debug)
        fprintf(stderr, "%d %d %lf %lf\n", i, j, on[i][j], off[i][j]);
    }
  }
  fclose(fpon);

  fprintf(stderr, "maxon = %f, maxoff = %f \n", maxon, maxoff);
  return 0;
}
