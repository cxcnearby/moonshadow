/** @file rstread.cc
 *  @author  zwang
 *  @date    2017.7
/* ========================================================= */
#include "TBranch.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TTree.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;

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

  string dir = argv[1];
  string file = dir + "/*.root";

  TChain *chain = new TChain("trec");
  chain->Add(file.c_str());

  chain->SetBranchAddress("mjd", &result.mjd);
  chain->SetBranchAddress("nhit", &result.nhit);
  chain->SetBranchAddress("npea", &result.npea);
  chain->SetBranchAddress("zen", &result.zen);
  chain->SetBranchAddress("azi", &result.azi);
  chain->SetBranchAddress("xc", &result.xc);
  chain->SetBranchAddress("yc", &result.yc);
  chain->SetBranchAddress("ndetc", &result.ndetc);
  chain->SetBranchAddress("nfitc", &result.nfitc);
  chain->SetBranchAddress("npec", &result.npec);
  chain->SetBranchAddress("zenc", &result.zenc);
  chain->SetBranchAddress("azic", &result.azic);
  chain->SetBranchAddress("omega", &result.omega);
  chain->SetBranchAddress("chi2p", &result.chi2p);
  chain->SetBranchAddress("chi2c", &result.chi2c);
  chain->SetBranchAddress("compactness", &result.compactness);
  chain->SetBranchAddress("pincness", &result.pincness);

  int Ev = chain->GetEntries();
  for (int i = 0; i < Ev; i++) {
    chain->GetEntry(i);
    // if (result.inout==1&&result.sigma<1&&result.sumpd>30)
    printf("%lf %d %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",
           result.mjd, result.nhit, result.npea, result.zen, result.azi,
           result.xc, result.yc, result.ndetc, result.nfitc, result.npec,
           result.zenc, result.azic, result.omega, result.chi2p, result.chi2c,
           result.compactness, result.pincness);
  }
  //  delete hfile;
  delete chain;
  return 0;
}
