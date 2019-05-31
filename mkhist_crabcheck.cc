#include "TTree.h"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
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

typedef struct {
  double mjd;
  int nhit;
  double npea;
  double zen;
  double azi;
  double xc;
  double yc;
  int ndetc;
  int nfitc;
  double npec;
  double zenc;
  double azic;
  double omega;
  double chi2p;
  double chi2c;
  double compactness;
  double pincness;
  int id;
  double theta;
  double phi;
  double xd;
  double yd;
} result_t;

int main(int argc, char *argv[]) {

  if (argc != 2) {
    printf("Usage :zcat input.gz | %s out.root\n", argv[0]);
    exit(1);
  }

  result_t result;

  int i, j;
  int id, xd, yd;

  TFile *f = new TFile(argv[1], "recreate");
  TTree *t = new TTree("t", "events within windows");
  t->Branch("mjd",&result.mjd);
  t->Branch("nhit",&result.nhit);
  t->Branch("npea",&result.npea);
  t->Branch("zen",&result.zen);
  t->Branch("azi",&result.azi);
  t->Branch("xc",&result.xc);
  t->Branch("yc",&result.yc);
  t->Branch("ndetc",&result.ndetc);
  t->Branch("nfitc",&result.nfitc);
  t->Branch("npec",&result.npec);
  t->Branch("zenc",&result.zenc);
  t->Branch("azic",&result.azic);
  t->Branch("omega",&result.omega);
  t->Branch("chi2p",&result.chi2p);
  t->Branch("chi2c",&result.chi2c);
  t->Branch("compactness",&result.compactness);
  t->Branch("pincness",&result.pincness);
  t->Branch("id",&result.id);
  t->Branch("theta",&result.theta);
  t->Branch("phi",&result.phi);
  t->Branch("xd",&result.xd);
  t->Branch("yd",&result.yd);

  while (scanf("%lf %d %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf "
               "%lf %d %lf %lf",
               &result.mjd, &result.nhit, &result.npea, &result.zen,
               &result.azi, &result.xc, &result.yc, &result.ndetc,
               &result.nfitc, &result.npec, &result.zenc, &result.azic,
               &result.omega, &result.chi2p, &result.chi2c, &result.compactness,
               &result.pincness, &result.id, &result.theta, &result.phi) > 0) {
    //    if(zen>=50) continue;

    result.xd = result.theta * sin(result.phi * deg_rad);
    result.yd = result.theta * cos(result.phi * deg_rad);
    // yd = floor( ( r * sin( theta * deg_rad ) + 5.0 ) * DIFm +0.5 );
    // xd = floor( ( r * cos( theta * deg_rad ) + 5.0 ) * DIFm +0.5 );
    t->Fill();
  }
  f->Write();
  f->Close();
  return 0;
}
