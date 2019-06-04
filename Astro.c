#include "Astro.h"
#include <math.h>
#include <stdio.h>

/********** GMST **********/

double gmst(double mjd) {

  double ut1, Tu, alpham;

  ut1 = mjd - floor(mjd);
  Tu = (mjd - 51544.5) / 36525.0;
  alpham = 280.460618370 + 36000.770053608 * Tu + 3.8793333331e-4 * Tu * Tu -
           2.5833333331e-8 * Tu * Tu * Tu;

  return (180.0 + ut1 * 360.0 + alpham);
}

double azirange(double azi) {
  if (azi > 180.0)
    azi -= 360.0;
  if (azi < -180.)
    azi += 360.0;
  return azi;
}

/********** transfomation of (zenith & azimuth)
            from (decrination & right asension)  **********/

void equator_horizon(double mjd, double ras, double dec, double *zen,
                     double *azi)
// double mjd, ras, dec, *zen, *azi;
{

  double lst, H, g;
  double sH, sdec, sphi;
  double cH, cdec, cphi;
  double zenith, azimuth;
  double sazi, cazi;

  g = gmst(mjd) + tibet_lo;
  lst = (g - floor(g / 360.0) * 360.0);
  H = lst - ras;

  sH = sin(H * deg_rad);
  sdec = sin(dec * deg_rad);
  sphi = sin(tibet_la * deg_rad);
  cH = cos(H * deg_rad);
  cdec = cos(dec * deg_rad);
  cphi = cos(tibet_la * deg_rad);

  zenith = acos(sphi * sdec + cphi * cdec * cH) * rad_deg;

  sazi = (-cdec * sH) / sin(zenith * deg_rad);
  cazi = (cphi * sdec - sphi * cdec * cH) / sin(zenith * deg_rad);

  azimuth = atan2(sazi, cazi) * rad_deg;
  if (azimuth > 180.0)
    azimuth -= 360.0;
  if (azimuth < -180.0)
    azimuth += 360.0;

  *zen = zenith;
  *azi = azimuth;
}

/********** transfomation of (decrination & right asension)
                                    from (zenith & azimuth) **********/

void horizon_equator(double mjd, double zen, double azi, double *ras,
                     double *dec)
// double mjd, *ras, *dec, zen, azi;
{

  double lst, H, g;
  double sazi, szen, sphi;
  double cazi, czen, cphi;
  double DEC, RAS;
  double sH, cH;

  g = gmst(mjd) + tibet_lo;
  lst = (g - floor(g / 360.0) * 360.0);

  sazi = sin(azi * deg_rad);
  szen = sin(zen * deg_rad);
  sphi = sin(tibet_la * deg_rad);
  cazi = cos(azi * deg_rad);
  czen = cos(zen * deg_rad);
  cphi = cos(tibet_la * deg_rad);

  DEC = asin(sphi * czen + cphi * szen * cazi) * rad_deg;

  sH = (-szen * sazi) / cos(DEC * deg_rad);
  cH = (cphi * czen - sphi * szen * cazi) / cos(DEC * deg_rad);

  H = atan2(sH, cH) * rad_deg;

  RAS = lst - H;
  if (RAS > 360.0)
    RAS -= 360.0;
  if (RAS < 0.0)
    RAS += 360.0;

  *ras = RAS;
  *dec = DEC;
}

/********** caluculation of distance
 *                          by (altitude & longitude)
 *                             altitude=0 at equator!!
 *                                                   *********/

double distance_equatorial(double ra1, double dec1, double ra2, double dec2) {
  return (acos(cos((ra2 - ra1) * deg_rad) * cos(dec1 * deg_rad) *
                   cos(dec2 * deg_rad) +
               sin(dec1 * deg_rad) * sin(dec2 * deg_rad)) *
          rad_deg);
}

/********** caluculation of distance
 *                          by (zenith & azimuth)
 *                             zenith=0 at north pole!!
 *                                                   *********/

double distance_horizontal(double zen1, double azi1, double zen2, double azi2) {
  return (acos(cos((azi1 - azi2) * deg_rad) * sin(zen1 * deg_rad) *
                   sin(zen2 * deg_rad) +
               cos(zen1 * deg_rad) * cos(zen2 * deg_rad)) *
          rad_deg);
}

/********** caluculation of position angle
 *
 *		the direction is point(ra2,dec2) relative to point(ra1,dec1)
 *		angle=0 when the direction is to the north pole, and clockwise.
 *
 * **********/

double direction_equatorial(double ra1, double dec1, double ra2, double dec2) {
  return (atan2(cos(dec2 * deg_rad) * sin((ra2 - ra1) * deg_rad),
                cos(dec1 * deg_rad) * sin(dec2 * deg_rad) -
                    sin(dec1 * deg_rad) * cos(dec2 * deg_rad) *
                        cos((ra2 - ra1) * deg_rad)) *
          rad_deg);
}

/********** caluculation of position angle
 *  *
 *   *              the direction is point(ra2,dec2) relative to point(ra1,dec1)
 *    *              angle=0 when the direction is to the north pole, and
 * clockwise.
 *     *
 *      * **********/

double direction_horizontal(double zen1, double azi1, double zen2,
                            double azi2) {
  return (atan2(sin(zen2 * deg_rad) * sin((azi1 - azi2) * deg_rad),
                sin(zen1 * deg_rad) * cos(zen2 * deg_rad) -
                    cos(zen1 * deg_rad) * sin(zen2 * deg_rad) *
                        cos((azi1 - azi2) * deg_rad)) *
          rad_deg);
}

/********** randam function **********/

double randam_fraction() {

  double randam;
  static long X = 53402397;

  X = X * 65539 + 125654;
  if (X < 0) {
    X += 2147483647;
    X += 1;
  }

  randam = X / 2147483648.0;
  return (0.001 * floor(randam * 10.0) - 0.005);
}

/**********  Bit Swap **********/

short swap2(short dd) { return (((dd & 0xff00) >> 8) | ((dd & 0x00ff) << 8)); }

unsigned short uns_swap2(unsigned short dd) {
  return (((dd & 0xff00) >> 8) | ((dd & 0x00ff) << 8));
}

unsigned int swap4(unsigned int dd) {
  return (((dd & 0xff000000) >> 24) | ((dd & 0x00ff0000) >> 8) |
          ((dd & 0x0000ff00) << 8) | ((dd & 0x000000ff) << 24));
}

/********** transfomation of (galactic Dec. & galactic RA.)
                        from (decrination & right asension)  **********/

void eguator_galaxcy(double ras, double dec, double *gras, double *gdec)
// double ras, dec, *gras, *gdec;
{

  double sr, sd, si;
  double cr, cd, ci;
  double sin_ras, cos_ras;

  sr = sin((ras - 282.85) * deg_rad);
  sd = sin(dec * deg_rad);
  si = sin(62.87 * deg_rad);
  cr = cos((ras - 282.85) * deg_rad);
  cd = cos(dec * deg_rad);
  ci = cos(62.87 * deg_rad);

  *gdec = asin((sd * ci) - (cd * sr * si)) * rad_deg;

  cos_ras = (cd * cr);
  sin_ras = (sd * si) + (cd * sr * ci);

  *gras = atan2(sin_ras, cos_ras) * rad_deg + 32.93;
  if (*gras > 180.0)
    *gras -= 360.0;
  if (*gras < -180.0)
    *gras += 360.0;
}

/********** transfomation of (decrination & right asension)
                        from (galactic Dec. & galactic RA.) **********/
void galaxcy_equator(double gras, double gdec, double *ras, double *dec)
// double *ras, *dec, gras, gdec;
{

  double sr, sd, si;
  double cr, cd, ci;
  double sin_gras, cos_gras;

  sr = sin((gras - 32.93) * deg_rad);
  sd = sin(gdec * deg_rad);
  si = sin(62.87 * deg_rad);
  cr = cos((gras - 32.93) * deg_rad);
  cd = cos(gdec * deg_rad);
  ci = cos(62.87 * deg_rad);

  *dec = asin((sd * ci) + (cd * sr * si)) * rad_deg;

  cos_gras = (cd * cr);
  sin_gras = (-(sd * si) + (cd * sr * ci));

  *ras = atan2(sin_gras, cos_gras) * rad_deg + 282.85;
  if (*ras > 180.0)
    *ras -= 360.0;
  if (*ras < -180.0)
    *ras += 360.0;
}
