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

struct s_wcdaevents {
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
};

int main(int argc, char *argv[]) {

  if (argc < 7) {
    printf("%s  InputPath  OutputFile.root  WindowRadius  bgWindowNumber  "
           "ZenithMax  SourceType  [source RA]  [source DEC]\n(SourceType: 1 "
           "for Moon, 2 for Sun, 0 need RA & DEC)\n",
           argv[0]);
    exit(0);
  }
  string sInputPath = argv[1];
  string sInputFile = sInputPath + "/*.root";
  string sOutputFile = argv[2];

  double dWindowRadius = atof(argv[3]);
  int iWindowNumberOneSide = atof(argv[4]);
  int iSourceType = atoi(argv[5]);
  double dZenithMax = atof(argv[6]);
  double dZenithMin = asin(sin(dWindowRadius * deg_rad) /
                           sin(PI / (2 * iWindowNumberOneSide + 1))) *
                      rad_deg;

  double dRA, dDEC, dZENITH, dAZI;

  if ((argc == 8) && (iSourceType == 0)) {
    dRA = atof(argv[6]);
    dDEC = atof(argv[7]);
  }

  s_wcdaevents s_Event;

  int id;
  double theta, phi, xd, yd;

  TFile *fSelected = new TFile(sOutputFile.c_str(), "recreate");
  TTree *tsel = new TTree("tsel", "selected events");
  tsel->Branch("mjd", &s_Event.mjd);
  tsel->Branch("nhit", &s_Event.nhit);
  tsel->Branch("npea", &s_Event.npea);
  tsel->Branch("zen", &s_Event.zen);
  tsel->Branch("azi", &s_Event.azi);
  tsel->Branch("xc", &s_Event.xc);
  tsel->Branch("yc", &s_Event.yc);
  tsel->Branch("ndetc", &s_Event.ndetc);
  tsel->Branch("nfitc", &s_Event.nfitc);
  tsel->Branch("npec", &s_Event.npec);
  tsel->Branch("zenc", &s_Event.zenc);
  tsel->Branch("azic", &s_Event.azic);
  tsel->Branch("omega", &s_Event.omega);
  tsel->Branch("chi2p", &s_Event.chi2p);
  tsel->Branch("chi2c", &s_Event.chi2c);
  tsel->Branch("compactness", &s_Event.compactness);
  tsel->Branch("pincness", &s_Event.pincness);
  tsel->Branch("id", &id);
  tsel->Branch("theta", &theta);
  tsel->Branch("phi", &phi);
  tsel->Branch("xd", &xd);
  tsel->Branch("yd", &yd);

  TChain *cInput = new TChain("trec");
  cInput->Add(sInputFile.c_str());

  cInput->SetBranchAddress("mjd", &s_Event.mjd);
  cInput->SetBranchAddress("nhit", &s_Event.nhit);
  cInput->SetBranchAddress("npea", &s_Event.npea);
  cInput->SetBranchAddress("zen", &s_Event.zen);
  cInput->SetBranchAddress("azi", &s_Event.azi);
  cInput->SetBranchAddress("xc", &s_Event.xc);
  cInput->SetBranchAddress("yc", &s_Event.yc);
  cInput->SetBranchAddress("ndetc", &s_Event.ndetc);
  cInput->SetBranchAddress("nfitc", &s_Event.nfitc);
  cInput->SetBranchAddress("npec", &s_Event.npec);
  cInput->SetBranchAddress("zenc", &s_Event.zenc);
  cInput->SetBranchAddress("azic", &s_Event.azic);
  cInput->SetBranchAddress("omega", &s_Event.omega);
  cInput->SetBranchAddress("chi2p", &s_Event.chi2p);
  cInput->SetBranchAddress("chi2c", &s_Event.chi2c);
  cInput->SetBranchAddress("compactness", &s_Event.compactness);
  cInput->SetBranchAddress("pincness", &s_Event.pincness);

  int nentries = cInput->GetEntriesFast();
  for (int i = 0; i < nentries; i++) {
    cInput->GetEntry(i);
    switch (iSourceType) {
    case 1:
      moon_orbit(s_Event.mjd, &dRA, &dDEC);
      break;
    case 2:
      sun_orbit(s_Event.mjd, &dRA, &dDEC);
      break;
    default:
      break;
    }
    equator_horizon(s_Event.mjd, dRA, dDEC, &dZENITH, &dAZI);
    if (dZENITH > dZenithMin && dZENITH < dZenithMax) {
      s_Event.azic = 60. - s_Event.azic;
      s_Event.azic = azirange(s_Event.azic);
      for (id = -iWindowNumberOneSide; id <= iWindowNumberOneSide; id++) {
        double dAziDiff =
            2 * asin(sin(dWindowRadius * deg_rad) / sin(dZENITH * deg_rad)) *
            rad_deg * id;
        double dAZIid = dAZI + dAziDiff;
        dAZIid = azirange(dAZIid);
        theta =
            distance_horizontal(dZENITH, dAZIid, s_Event.zenc, s_Event.azic);
        if (theta < dWindowRadius) {
          double daziid = s_Event.azic - dAziDiff;
          daziid = azirange(daziid);
          double dra, ddec;
          horizon_equator(s_Event.mjd, s_Event.zenc, daziid, &dra, &ddec);
          phi = direction_equatorial(dRA, dDEC, dra, ddec);
          xd = theta * sin(phi * deg_rad);
          yd = theta * cos(phi * deg_rad);
          tsel->Fill();
          break;
        }
      }
    }
  }
  //  delete hfile;
  delete cInput;
  fSelected->Write();
  fSelected->Close();
  return 0;
}
