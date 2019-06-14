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

  FILE *fpLog;
  string::size_type pos = sOutputFile.rfind(".root");
  string sFpLog = sOutputFile.substr(0, pos) + ".log";
  if ((fpLog = fopen(sFpLog.c_str(), "w")) == NULL)
    printf("cannot open log file\n");
  double dWindowRadius = atof(argv[3]);
  int iWindowNumberOneSide = atof(argv[4]);
  double dZenithMax = atof(argv[5]);
  double dZenithMin = asin(sin(dWindowRadius * deg_rad) /
                           sin(PI / (2 * iWindowNumberOneSide + 1))) *
                      rad_deg;
  int iSourceType = atoi(argv[6]);

  double dRA, dDEC, dZENITH, dAZI;

  if ((argc == 9) && (iSourceType == 0)) {
    dRA = atof(argv[7]);
    dDEC = atof(argv[8]);
  }

  s_wcdaevents s_Event;

  int id;
  double theta, phi, xd, yd;
  clock_t ctStart, ctFinish;
  //  double dLiveTime = 0.;
  double dLiveTime1 = 0.;
  //  double dLT0 = 0.;
  //  double dLT1 = 0.;
  //  double dExposureTime = 0.;
  double dExposureTime1 = 0.;
  //  double dET0 = 0.;
  //  double dET1 = 0.;

  ctStart = clock();

  TFile *fSelected = new TFile(sOutputFile.c_str(), "recreate");
  TH1F *hlive =
      new TH1F("hlive", "events time dist. 1e-5 day per bin", 200000, 0., 2.);
  TH1F *hexpo =
      new TH1F("hexpo", "events time dist. 1e-5 day per bin", 200000, 0., 2.);
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

  Long64_t nentries = cInput->GetEntries();
  int iTimeInte;
  for (Long64_t i = 0; i < nentries; i++) {
    cInput->GetEntry(i);
    if (i == 0)
      iTimeInte = int(s_Event.mjd);
    double dTimeDeci = s_Event.mjd - iTimeInte;
    hlive->Fill(dTimeDeci);
    //    if (s_Event.mjd - dLT1 > 0 && s_Event.mjd - dLT1 < 2e-7) {
    //      dLT1 = s_Event.mjd;
    //    } else {
    //      dLiveTime += dLT1 - dLT0;
    //      dLT0 = s_Event.mjd;
    //      dLT1 = dLT0;
    //    }
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
      hexpo->Fill(dTimeDeci);
      //      if (s_Event.mjd - dET1 > 0 && s_Event.mjd - dET1 < 2e-7) {
      //        dET1 = s_Event.mjd;
      //      } else {
      //        dExposureTime += dET1 - dET0;
      //        dET0 = s_Event.mjd;
      //        dET1 = dET0;
      //      }
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
  for (int i = 1; i <= 200000; i++) {
    double dLive = hlive->GetBinContent(i);
    double dExpo = hexpo->GetBinContent(i);
    if (dLive > 0.1)
      dLiveTime1 += 1e-5;
    if (dExpo > 0.1)
      dExposureTime1 += 1e-5;
  }
  delete hlive;
  delete hexpo;
  delete cInput;
  tsel->Write("", TObject::kOverwrite);
  fSelected->Close();
  ctFinish = clock();
  fprintf(fpLog, "%f %f\n", dLiveTime1, dExposureTime1);
  cout << "LiveTime: " << dLiveTime1 << " d; ExpoTime: " << dExposureTime1
       << " d. total events: " << nentries << ". selecting use "
       << double((ctFinish - ctStart) / CLOCKS_PER_SEC) << " s." << endl;
  return 0;
}
