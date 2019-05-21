#include "Astro.c"
#include "Astro.h"
#include "Moon.c"
#include "Sun.c"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define debug_moon_pos 0
#define R_win 5  /* radius of the window */
#define id_max 4 /* the number of off-source windows in each side */

#define deg_rad 0.017453292519943296 /*   pi/180  */
#define rad_deg 57.295779513082322   /*   180/pi  */
#define pi 3.14159265358979323846    /* pi */

int main() {

  /*	for selection	*/
  double THETA, THETA_min, PHI, PHI_off, theta, phi, RA, DEC, phi_off, ra, dec;

  /*	data buffer	*/
  int nhit, ndetc, nfitc;
  double mjd, npea, zen, azi, xc, yc, npec, zenc, azic, omega, chi2p, chi2c,
      compactness, pincness;

  THETA_min = asin(sin(R_win * deg_rad) / sin(pi / (2 * id_max + 1))) * rad_deg;
  /*	sin c = sin a / sin A  , C=90deg, and C is the point of contact of two
   * adjacent circles	*/
  // printf("The windows will overlap when THETA_Moon < %lf ! \n",THETA_min);

  /*	read data from experiment data	*/
  while (
      (scanf("%lf %d %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf",
             &mjd, &nhit, &npea, &zen, &azi, &xc, &yc, &ndetc, &nfitc, &npec,
             &zenc, &azic, &omega, &chi2p, &chi2c, &compactness, &pincness)) >
      0) {

    /********** Moon Orbit **********/
    moon_orbit(mjd, &RA, &DEC);

    /**********  Theta & Phi of Moon *********/
    equator_horizon(mjd, RA, DEC, &THETA, &PHI);

    azic = 60. - azic;

    if (azic > 180.)
      azic -= 360.;
    if (azic < -180.)
      azic += 360.;

    if (debug_moon_pos) {
      printf("The windows will overlap when THETA_Moon < %lf ! %lf %lf %lf %lf "
             "%lf \r",
             THETA_min, mjd, RA, DEC, THETA,
             PHI); /*	the boundary for non-overlap and check Moon's position
                    */
      continue;
    }

    /********** avoid windows overlap **********/
    if (THETA > THETA_min && THETA < 45. && THETA > 0.)
    // if(THETA > 7.8 && THETA < 70)
    {
      for (int id = -id_max; id <= id_max; id++) {
        /*	to calculate the coordinate PHI after equal-theta shift	*/
        PHI_off = PHI + 2 * asin(sin(R_win * deg_rad) / sin(THETA * deg_rad)) *
                            rad_deg * id;
        /*	sin c = sin a / sin A  , C=90deg, and C is the point of contact
         * of two adjacent circles	*/

        if (PHI_off > 180.0)
          PHI_off -= 360.0;

        if (PHI_off < -180.0)
          PHI_off += 360.0;

        theta = distance_horizontal(THETA, PHI_off, zenc, azic);

        if (theta < R_win) {
          //						horizon_equator( mjd,
          //THETA, PHI_off, &RA, &DEC );
          // azi=direction_equatorial(RA,DEC,ra,dec);
          phi_off =
              azic - 2 * asin(sin(R_win * deg_rad) / sin(THETA * deg_rad)) *
                         rad_deg * id;
          if (phi_off > 180.)
            phi_off -= 360.;
          if (phi_off < -180.)
            phi_off += 360.;

          horizon_equator(mjd, zenc, phi_off, &ra, &dec);
          phi = direction_equatorial(RA, DEC, ra, dec);
          printf("%lf %d %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf "
                 "%lf %d %lf %lf\n",
                 mjd, nhit, npea, zen, azi, xc, yc, ndetc, nfitc, npec, zenc,
                 azic, omega, chi2p, chi2c, compactness, pincness, id, theta,
                 phi);
        }
      }
    }
  }
  return 0;
}
