/** @file evtselect.cc
 *  @author  changxc
 *  @date    2019.06.01
/* ========================================================= */
#include "Astro.c"
#include "Astro.h"
#include "Moon.c"
#include "Sun.c"
#include "TBranch.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TTree.h"
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

struct s_wcdaeventsSel {
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
};

int main(int argc, char *argv[]) {

  if (argc < 4) {
    printf("%s  InputPath  OutputFile.root  nfitc_min\nUse \\* instead of * if appearing in InputFile.root\n",
           argv[0]);
    exit(0);
  }
  string sInputFile = argv[1];
  string sOutputFile = argv[2];
  int iNFitcMin = atoi(argv[3]);

  s_wcdaeventsSel s_EventSel;

  TFile *fSelected = new TFile(sOutputFile.c_str(), "recreate");
  TTree *tsel = new TTree("tsel", "nfitc selected events");
  tsel->Branch("mjd", &s_EventSel.mjd);
  tsel->Branch("nhit", &s_EventSel.nhit);
  tsel->Branch("npea", &s_EventSel.npea);
  tsel->Branch("zen", &s_EventSel.zen);
  tsel->Branch("azi", &s_EventSel.azi);
  tsel->Branch("xc", &s_EventSel.xc);
  tsel->Branch("yc", &s_EventSel.yc);
  tsel->Branch("ndetc", &s_EventSel.ndetc);
  tsel->Branch("nfitc", &s_EventSel.nfitc);
  tsel->Branch("npec", &s_EventSel.npec);
  tsel->Branch("zenc", &s_EventSel.zenc);
  tsel->Branch("azic", &s_EventSel.azic);
  tsel->Branch("omega", &s_EventSel.omega);
  tsel->Branch("chi2p", &s_EventSel.chi2p);
  tsel->Branch("chi2c", &s_EventSel.chi2c);
  tsel->Branch("compactness", &s_EventSel.compactness);
  tsel->Branch("pincness", &s_EventSel.pincness);
  tsel->Branch("id", &s_EventSel.id);
  tsel->Branch("theta", &s_EventSel.theta);
  tsel->Branch("phi", &s_EventSel.phi);
  tsel->Branch("xd", &s_EventSel.xd);
  tsel->Branch("yd", &s_EventSel.yd);

  TChain *cInput = new TChain("tsel");
  cInput->Add(sInputFile.c_str());

  cInput->SetBranchAddress("mjd", &s_EventSel.mjd);
  cInput->SetBranchAddress("nhit", &s_EventSel.nhit);
  cInput->SetBranchAddress("npea", &s_EventSel.npea);
  cInput->SetBranchAddress("zen", &s_EventSel.zen);
  cInput->SetBranchAddress("azi", &s_EventSel.azi);
  cInput->SetBranchAddress("xc", &s_EventSel.xc);
  cInput->SetBranchAddress("yc", &s_EventSel.yc);
  cInput->SetBranchAddress("ndetc", &s_EventSel.ndetc);
  cInput->SetBranchAddress("nfitc", &s_EventSel.nfitc);
  cInput->SetBranchAddress("npec", &s_EventSel.npec);
  cInput->SetBranchAddress("zenc", &s_EventSel.zenc);
  cInput->SetBranchAddress("azic", &s_EventSel.azic);
  cInput->SetBranchAddress("omega", &s_EventSel.omega);
  cInput->SetBranchAddress("chi2p", &s_EventSel.chi2p);
  cInput->SetBranchAddress("chi2c", &s_EventSel.chi2c);
  cInput->SetBranchAddress("compactness", &s_EventSel.compactness);
  cInput->SetBranchAddress("pincness", &s_EventSel.pincness);
  cInput->SetBranchAddress("id", &s_EventSel.id);
  cInput->SetBranchAddress("theta", &s_EventSel.theta);
  cInput->SetBranchAddress("phi", &s_EventSel.phi);
  cInput->SetBranchAddress("xd", &s_EventSel.xd);
  cInput->SetBranchAddress("yd", &s_EventSel.yd);


  Long64_t nentries = cInput->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
    cInput->GetEntry(i);
    if (s_EventSel.nfitc >= iNFitcMin)
      tsel->Fill();
  }
  delete cInput;
  tsel->Write("", TObject::kOverwrite);
  fSelected->Close();
  return 0;
}
