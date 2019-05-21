/** @file rstread.cc
 *  @author  zwang
 *  @date    2017.7
/* ========================================================= */
#include "TBranch.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TTree.h"
#include <stdio.h>
#include <stdlib.h>

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

} result_t;

int main(int argc, char *argv[]) {
  FILE *fp;
  result_t result;

  TFile *hfile = new TFile(argv[1]);
  TTree *rst = (TTree *)hfile->Get("trec");
  rst->SetBranchAddress("mjd", &result.mjd);
  rst->SetBranchAddress("nhit", &result.nhit);
  rst->SetBranchAddress("npea", &result.npea);
  rst->SetBranchAddress("zen", &result.zen);
  rst->SetBranchAddress("azi", &result.azi);
  rst->SetBranchAddress("xc", &result.xc);
  rst->SetBranchAddress("yc", &result.yc);
  rst->SetBranchAddress("ndetc", &result.ndetc);
  rst->SetBranchAddress("nfitc", &result.nfitc);
  rst->SetBranchAddress("npec", &result.npec);
  rst->SetBranchAddress("zenc", &result.zenc);
  rst->SetBranchAddress("azic", &result.azic);
  rst->SetBranchAddress("omega", &result.omega);
  rst->SetBranchAddress("chi2p", &result.chi2p);
  rst->SetBranchAddress("chi2c", &result.chi2c);
  rst->SetBranchAddress("compactness", &result.compactness);
  rst->SetBranchAddress("pincness", &result.pincness);

  int Ev = rst->GetEntries();
  for (int i = 0; i < Ev; i++) {
    rst->GetEntry(i);
    // if (result.inout==1&&result.sigma<1&&result.sumpd>30)
    printf("%lf %d %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",
           result.mjd, result.nhit, result.npea, result.zen, result.azi,
           result.xc, result.yc, result.ndetc, result.nfitc, result.npec,
           result.zenc, result.azic, result.omega, result.chi2p, result.chi2c,
           result.compactness, result.pincness);
  }
  //  delete hfile;
  delete rst;
  return 0;
}
