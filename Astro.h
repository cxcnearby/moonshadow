
/********** Prottotype Diclaration **********/

void               moon_orbit(double mjd,double *MRAS,double *MDEC);
void               sun_orbit(double mjd,double *SRAS,double *SDEC);
double             gmst( double mjd );
double             azirange(double azi);
void               equator_horizon( double mjd, double ras, double dec, double *zen, double *azi );
void               horizon_equator( double mjd, double zen, double azi, double *ras, double *dec );
void               galaxcy_equator( double gras, double gdec, double *ras, double *dec );
void               equator_galaxcy( double ras, double dec, double *gras, double *gdec );
double             distance_equatorial(double ra1, double dec1, double ra2, double dec2);
double             distance_horizontal(double zen1, double azi1, double zen2, double azi2);
double             direction_equatorial( double ra1, double dec1 ,double ra2, double dec2 );
double             direction_horizontal( double zen1, double azi1 ,double zen2, double azi2 );
double             randam_fraction();
short              swap2(short dd);
unsigned short     uns_swap2(unsigned short dd);
unsigned int       swap4(unsigned int dd);


/********** Constant **********/

#define deg_rad 0.017453292519943296      /*   pi/180  */
#define rad_deg 57.295779513082322      /*   180/pi  */
#define pi 3.14159265358979323846       /* pi */


/********** Location of Tibet **********/

#define  tibet_la 29.3585          /*  latitude  */
#define  tibet_lo 100.1374          /* longtitude */  

