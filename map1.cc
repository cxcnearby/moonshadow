#include "TBranch.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include <TH1.h>
#include <TH2.h>
#include "TROOT.h"
#include "TTree.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Astro.h"
#include "ener_reso.h"


using namespace std;

double sig_smooth(vector<vector<float>> &vOn, vector<vector<float>> &vOff, int x, int y, int id, double binwidth, int nbin,
                  double r_smooth, double *Non, double *Noff);
double gaus32(int kk, double r);

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

int nk = 100;

int main(int argc, char *argv[]) {

  if (argc < 12) {
    printf("%s  InputFile.root  OutputFile.root  WindowRadius  bgWindowNumber  BinWidth  nfitc_min  nfitc_max  Compactness1  Compactness2  Compactness3  SmoothRadius1  SmoothRadius2  SmoothRadius3\nUse \\* instead of * if appearing in InputFile.root\n",
           argv[0]);
    exit(0);
  }
  string sInputFile = argv[1];
  string sOutputFile = argv[2];
  double dWindowRadius = atof(argv[3]);
  int iWindowNumberOneSide = atof(argv[4]);
  double dBinWidth = atof(argv[5]);
  double dBinWidth1 = 0.1;
  double dNFitcMin = atof(argv[6]);
  double dNFitcMax = atof(argv[7]);
  double dCompactness[3] = {atof(argv[8]), atof(argv[9]), atof(argv[10])};
  double dSmoothRadius[3] = {atof(argv[11]), atof(argv[12]), atof(argv[13])};
  int iBinNumber = int(dWindowRadius / dBinWidth * 2.0 + 1.0);

  vector<vector<float>> vOn1(iBinNumber + 1, vector<float>(iBinNumber + 1));
  vector<vector<float>> vOff1(iBinNumber + 1, vector<float>(iBinNumber + 1));
  vector<vector<float>> vOn2(iBinNumber + 1, vector<float>(iBinNumber + 1));
  vector<vector<float>> vOff2(iBinNumber + 1, vector<float>(iBinNumber + 1));
  vector<vector<float>> vOn3(iBinNumber + 1, vector<float>(iBinNumber + 1));
  vector<vector<float>> vOff3(iBinNumber + 1, vector<float>(iBinNumber + 1));

  s_wcdaeventsSel s_EventSel;
  double dMapRange = 3.0;
  int iMapBinNumber = int(dMapRange / dBinWidth * 2.0 + 1.0);
  int iBinShift = (iBinNumber - iMapBinNumber) / 2;
  double dSig, dSigSmooth;
  double dNOn, dNOff;

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

  TFile *fMap = new TFile(sOutputFile.c_str(), "recreate");
  TH2F *hOn1 = new TH2F(
      "hon1", "on events", iBinNumber, -(dWindowRadius + dBinWidth / 2.0),
      dWindowRadius + dBinWidth / 2.0, iBinNumber,
      -(dWindowRadius + dBinWidth / 2.0), dWindowRadius + dBinWidth / 2.0);
  TH2F *hOff1 = new TH2F(
      "hoff1", "off events", iBinNumber, -(dWindowRadius + dBinWidth / 2.0),
      dWindowRadius + dBinWidth / 2.0, iBinNumber,
      -(dWindowRadius + dBinWidth / 2.0), dWindowRadius + dBinWidth / 2.0);
  TH2F *hExc1 = new TH2F(
      "hexc1", "on-off events", iBinNumber, -(dWindowRadius + dBinWidth / 2.0),
      dWindowRadius + dBinWidth / 2.0, iBinNumber,
      -(dWindowRadius + dBinWidth / 2.0), dWindowRadius + dBinWidth / 2.0);
  TH2F *hOn2 = new TH2F(
      "hon2", "on events", iBinNumber, -(dWindowRadius + dBinWidth / 2.0),
      dWindowRadius + dBinWidth / 2.0, iBinNumber,
      -(dWindowRadius + dBinWidth / 2.0), dWindowRadius + dBinWidth / 2.0);
  TH2F *hOff2 = new TH2F(
      "hoff2", "off events", iBinNumber, -(dWindowRadius + dBinWidth / 2.0),
      dWindowRadius + dBinWidth / 2.0, iBinNumber,
      -(dWindowRadius + dBinWidth / 2.0), dWindowRadius + dBinWidth / 2.0);
  TH2F *hExc2 = new TH2F(
      "hexc2", "on-off events", iBinNumber, -(dWindowRadius + dBinWidth / 2.0),
      dWindowRadius + dBinWidth / 2.0, iBinNumber,
      -(dWindowRadius + dBinWidth / 2.0), dWindowRadius + dBinWidth / 2.0);
  TH2F *hOn3 = new TH2F(
      "hon3", "on events", iBinNumber, -(dWindowRadius + dBinWidth / 2.0),
      dWindowRadius + dBinWidth / 2.0, iBinNumber,
      -(dWindowRadius + dBinWidth / 2.0), dWindowRadius + dBinWidth / 2.0);
  TH2F *hOff3 = new TH2F(
      "hoff3", "off events", iBinNumber, -(dWindowRadius + dBinWidth / 2.0),
      dWindowRadius + dBinWidth / 2.0, iBinNumber,
      -(dWindowRadius + dBinWidth / 2.0), dWindowRadius + dBinWidth / 2.0);
  TH2F *hExc3 = new TH2F(
      "hexc3", "on-off events", iBinNumber, -(dWindowRadius + dBinWidth / 2.0),
      dWindowRadius + dBinWidth / 2.0, iBinNumber,
      -(dWindowRadius + dBinWidth / 2.0), dWindowRadius + dBinWidth / 2.0);
  // TH2F *hExc = new TH2F(
  //     "hexc", "on-off events", iMapBinNumber, -(dMapRange + dBinWidth / 2.0),
  //     dMapRange + dBinWidth / 2.0, iMapBinNumber,
  //     -(dMapRange + dBinWidth / 2.0), dMapRange + dBinWidth / 2.0);
  TH2F *hSigS2Dc1r1 = new TH2F("hsigs2dc1r1", "2-D significance (smoothed)", iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0, iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0);
  TH2F *hSigS2Dc1r2 = new TH2F("hsigs2dc1r2", "2-D significance (smoothed)", iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0, iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0);
  TH2F *hSigS2Dc1r3 = new TH2F("hsigs2dc1r3", "2-D significance (smoothed)", iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0, iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0);
  TH2F *hSigS2Dc2r1 = new TH2F("hsigs2dc2r1", "2-D significance (smoothed)", iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0, iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0);
  TH2F *hSigS2Dc2r2 = new TH2F("hsigs2dc2r2", "2-D significance (smoothed)", iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0, iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0);
  TH2F *hSigS2Dc2r3 = new TH2F("hsigs2dc2r3", "2-D significance (smoothed)", iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0, iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0);
  TH2F *hSigS2Dc3r1 = new TH2F("hsigs2dc3r1", "2-D significance (smoothed)", iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0, iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0);
  TH2F *hSigS2Dc3r2 = new TH2F("hsigs2dc3r2", "2-D significance (smoothed)", iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0, iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0);
  TH2F *hSigS2Dc3r3 = new TH2F("hsigs2dc3r3", "2-D significance (smoothed)", iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0, iMapBinNumber,
                          -(dMapRange + dBinWidth / 2.0),
                          dMapRange + dBinWidth / 2.0);

  Long64_t nentries = cInput->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
    cInput->GetEntry(i);
    if (s_EventSel.nfitc >= dNFitcMin && s_EventSel.nfitc < dNFitcMax && s_EventSel.compactness < 500) {
      if (s_EventSel.compactness < dCompactness[0]) {
        if (s_EventSel.id == 0) {
          hOn1->Fill(s_EventSel.xd, s_EventSel.yd);
        } else {
          hOff1->Fill(s_EventSel.xd, s_EventSel.yd);
        }
      }
      if (s_EventSel.compactness > dCompactness[0] && s_EventSel.compactness < dCompactness[1]) {
        if (s_EventSel.id == 0) {
          hOn2->Fill(s_EventSel.xd, s_EventSel.yd);
        } else {
          hOff2->Fill(s_EventSel.xd, s_EventSel.yd);
        }
      }
      if (s_EventSel.compactness > dCompactness[1]) {
        if (s_EventSel.id == 0) {
          hOn3->Fill(s_EventSel.xd, s_EventSel.yd);
        } else {
          hOff3->Fill(s_EventSel.xd, s_EventSel.yd);
        }
      }
    }
  }
  delete cInput;
  hOn1->Write();
  hOff1->Write();
  hOn2->Write();
  hOff2->Write();
  hOn3->Write();
  hOff3->Write();

  for (int i = 1; i <= iBinNumber; i++) {
    for (int j = 1; j <= iBinNumber; j++) {
      vOn1[i][j] = hOn1->GetBinContent(i, j);
      vOff1[i][j] = hOff1->GetBinContent(i, j);
      hExc1->SetBinContent(i, j, vOn1[i][j] - vOff1[i][j] / (2 * iWindowNumberOneSide));
      vOn2[i][j] = hOn2->GetBinContent(i, j);
      vOff2[i][j] = hOff2->GetBinContent(i, j);
      hExc2->SetBinContent(i, j, vOn2[i][j] - vOff2[i][j] / (2 * iWindowNumberOneSide));
      vOn3[i][j] = hOn3->GetBinContent(i, j);
      vOff3[i][j] = hOff3->GetBinContent(i, j);
      hExc3->SetBinContent(i, j, vOn3[i][j] - vOff3[i][j] / (2 * iWindowNumberOneSide));
    }
  }
  hExc1->Write();
  hExc2->Write();
  hExc3->Write();

  for (int k = 0; k < 3 ; k++) {
    for (int i = 1; i <= iMapBinNumber; i++) {
      for (int j = 1; j <= iMapBinNumber; j++) {
        dSigSmooth = sig_smooth(vOn1, vOff1, i + iBinShift, j + iBinShift, iWindowNumberOneSide, dBinWidth,
                                 iBinNumber, dSmoothRadius[k], &dNOn, &dNOff);
        // hExc->SetBinContent(i, j,
        //                     vOn[i + iBinShift][j + iBinShift] - vOff[i + iBinShift][j + iBinShift] / (2 * iWindowNumberOneSide));
        if (vOn1[i + iBinShift][j + iBinShift] != 0)
          dSig = (vOn1[i + iBinShift][j + iBinShift] - vOff1[i + iBinShift][j + iBinShift] / (2 * iWindowNumberOneSide)) /
                 sqrt(vOn1[i + iBinShift][j + iBinShift] + vOff1[i + iBinShift][j + iBinShift] / (2 * iWindowNumberOneSide * 2 *
                                                iWindowNumberOneSide));
        else
          dSig = -100;
        switch (k) {
        case 0:
          hSigS2Dc1r1->SetBinContent(i, j, dSigSmooth);
                if (i == (iMapBinNumber + 1) / 2 && j == (iMapBinNumber + 1) / 2)
        printf("center %f %f %f %f\n", dSmoothRadius[k], dSigSmooth, dNOn, dNOff);

          break;
        case 1:
          hSigS2Dc1r2->SetBinContent(i, j, dSigSmooth);
                if (i == (iMapBinNumber + 1) / 2 && j == (iMapBinNumber + 1) / 2)
        printf("center %f %f %f %f\n", dSmoothRadius[k], dSigSmooth, dNOn, dNOff);

          break;
        case 2:
          hSigS2Dc1r3->SetBinContent(i, j, dSigSmooth);
                if (i == (iMapBinNumber + 1) / 2 && j == (iMapBinNumber + 1) / 2)
        printf("center %f %f %f %f\n", dSmoothRadius[k], dSigSmooth, dNOn, dNOff);

          break;
        default:
          break;
        }
      }
    }
  }
  hSigS2Dc1r1->Write();
  hSigS2Dc1r2->Write();
  hSigS2Dc1r3->Write();
  for (int k = 0; k < 3 ; k++) {
    for (int i = 1; i <= iMapBinNumber; i++) {
      for (int j = 1; j <= iMapBinNumber; j++) {
        dSigSmooth = sig_smooth(vOn2, vOff2, i + iBinShift, j + iBinShift, iWindowNumberOneSide, dBinWidth,
                                 iBinNumber, dSmoothRadius[k], &dNOn, &dNOff);
        // hExc->SetBinContent(i, j,
        //                     vOn[i + iBinShift][j + iBinShift] - vOff[i + iBinShift][j + iBinShift] / (2 * iWindowNumberOneSide));
        if (vOn2[i + iBinShift][j + iBinShift] != 0)
          dSig = (vOn2[i + iBinShift][j + iBinShift] - vOff2[i + iBinShift][j + iBinShift] / (2 * iWindowNumberOneSide)) /
                 sqrt(vOn2[i + iBinShift][j + iBinShift] + vOff2[i + iBinShift][j + iBinShift] / (2 * iWindowNumberOneSide * 2 *
                                                iWindowNumberOneSide));
        else
          dSig = -100;
        switch (k) {
        case 0:
          hSigS2Dc2r1->SetBinContent(i, j, dSigSmooth);
                if (i == (iMapBinNumber + 1) / 2 && j == (iMapBinNumber + 1) / 2)
        printf("center %f %f %f %f\n", dSmoothRadius[k], dSigSmooth, dNOn, dNOff);

          break;
        case 1:
          hSigS2Dc2r2->SetBinContent(i, j, dSigSmooth);
                if (i == (iMapBinNumber + 1) / 2 && j == (iMapBinNumber + 1) / 2)
        printf("center %f %f %f %f\n", dSmoothRadius[k], dSigSmooth, dNOn, dNOff);

          break;
        case 2:
          hSigS2Dc2r3->SetBinContent(i, j, dSigSmooth);
                if (i == (iMapBinNumber + 1) / 2 && j == (iMapBinNumber + 1) / 2)
        printf("center %f %f %f %f\n", dSmoothRadius[k], dSigSmooth, dNOn, dNOff);

          break;
        default:
          break;
        }
      }
    }
  }
  hSigS2Dc2r1->Write();
  hSigS2Dc2r2->Write();
  hSigS2Dc2r3->Write();

  for (int k = 0; k < 3 ; k++) {
    for (int i = 1; i <= iMapBinNumber; i++) {
      for (int j = 1; j <= iMapBinNumber; j++) {
        dSigSmooth = sig_smooth(vOn3, vOff3, i + iBinShift, j + iBinShift, iWindowNumberOneSide, dBinWidth,
                                 iBinNumber, dSmoothRadius[k], &dNOn, &dNOff);
        // hExc->SetBinContent(i, j,
        //                     vOn[i + iBinShift][j + iBinShift] - vOff[i + iBinShift][j + iBinShift] / (2 * iWindowNumberOneSide));
        if (vOn3[i + iBinShift][j + iBinShift] != 0)
          dSig = (vOn3[i + iBinShift][j + iBinShift] - vOff3[i + iBinShift][j + iBinShift] / (2 * iWindowNumberOneSide)) /
                 sqrt(vOn3[i + iBinShift][j + iBinShift] + vOff3[i + iBinShift][j + iBinShift] / (2 * iWindowNumberOneSide * 2 *
                                                iWindowNumberOneSide));
        else
          dSig = -100;
        switch (k) {
        case 0:
          hSigS2Dc3r1->SetBinContent(i, j, dSigSmooth);
                if (i == (iMapBinNumber + 1) / 2 && j == (iMapBinNumber + 1) / 2)
        printf("center %f %f %f %f\n", dSmoothRadius[k], dSigSmooth, dNOn, dNOff);

          break;
        case 1:
          hSigS2Dc3r2->SetBinContent(i, j, dSigSmooth);
                if (i == (iMapBinNumber + 1) / 2 && j == (iMapBinNumber + 1) / 2)
        printf("center %f %f %f %f\n", dSmoothRadius[k], dSigSmooth, dNOn, dNOff);

          break;
        case 2:
          hSigS2Dc3r3->SetBinContent(i, j, dSigSmooth);
                if (i == (iMapBinNumber + 1) / 2 && j == (iMapBinNumber + 1) / 2)
        printf("center %f %f %f %f\n", dSmoothRadius[k], dSigSmooth, dNOn, dNOff);

          break;
        default:
          break;
        }
      }
    }
  }
  hSigS2Dc3r1->Write();
  hSigS2Dc3r2->Write();
  hSigS2Dc3r3->Write();

  fMap->Close();
  return 0;
}

double sig_smooth(vector<vector<float>> &vOn, vector<vector<float>> &vOff, int x, int y, int id, double binwidth, int nbin,
                  double r_smooth, double *Non, double *Noff) {
  int n_r;
  int i, j, i0, j0, i1, j1;
  double w = 0, r = 0;
  double dN = 0, dNoff = 0;
  double dis_square;
  double r_smooth_square = r_smooth * r_smooth;
  n_r = int(r_smooth / binwidth);

  *Non = 0.;
  *Noff = 0.;

  i0 = x - n_r;
  i1 = x + n_r;
  j0 = y - n_r;
  j1 = y + n_r;

  if (i0 < 1)
    i0 = 1;
  if (i1 > nbin)
    i1 = nbin;

  if (j0 < 1)
    j0 = 1;
  if (j1 > nbin)
    j1 = nbin;

  for (i = i0; i <= i1; i++) {
    for (j = j0; j <= j1; j++) {

      dis_square =
          ((i - x) * (i - x) + (j - y) * (j - y)) * binwidth * binwidth;
      if (dis_square > r_smooth_square)
        continue;

      if (nk > 12)
        w = 1;
      else
        w = gaus32(nk, sqrt(dis_square) * deg_rad);
      *Non += vOn[i][j] * w;
      *Noff += vOff[i][j] * w;
      dN += vOn[i][j] * w * w;
      dNoff += vOff[i][j] * w * w;
    }
  }

  if (dN != 0)
    return (*Non - *Noff / (2 * id)) / sqrt(dN + dNoff / (2 * 2 * id * id));
  else
    return -100;
}

double gaus32(int kk, double r) {
  int i;
  double result = 0;
  double reso2 = 0;

  for (i = 0; i < 3; i++) {
    reso2 = par[kk][i] * par[kk][i] * deg_rad * deg_rad;
    result += par[kk][i + 3] / (2 * PI * reso2) * exp(-r * r / (2 * reso2));
  }
  return result;
}
