#include "ener_reso.h"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define deg_rad 0.017453292519943296 /*   pi/180  */
#define rad_deg 57.295779513082322   /*   180/pi  */
#define pi 3.14159265358979323846    /* pi */

#define R_win 5.0
#define bin_width 0.05
#define bin_num (int)(R_win / bin_width * 2.0 + 1.0)

#define fig_x 3.0
#define fig_y 3.0
#define fig_x_bin_num (int)(2.0 * fig_x / bin_width + 1.0)
#define fig_y_bin_num (int)(2.0 * fig_y / bin_width + 1.0)

#define id_max 4

double on[bin_num][bin_num], off[bin_num][bin_num];
double sig_smooth3(int x, int y, double r_smooth, double *Non, double *Noff);
double gaus32(int kk, double r);

int Flag = 0;
int nk = 100;
int main(int argc, char *argv[]) {
  double reso = 1.0;
  double Non, Noff, sig;
  double max = 0;
  int i, j;

  if (argc < 2) {
    fprintf(stderr, "Usage :cat input | %s out.put [reso] [nk] \n", argv[0]);
    exit(1);
  }
  if (argc > 2)
    reso = atof(argv[2]);
  if (argc > 3)
    nk = atoi(argv[3]);
  fprintf(stderr, "reso= %f nk= %d \n", reso, nk);
  if (argc > 4)
    Flag = atoi(argv[4]);
  fprintf(stderr, "Flag= %d \n", Flag);
  gaus_init(Flag);

  for (i = 0; i < bin_num; i++) {
    for (j = 0; j < bin_num; j++) {
      fscanf(stdin, "%*d %*d %lf %lf", &on[i][j], &off[i][j]);
      max += on[i][j] + off[i][j];
    }
  }
  fprintf(stderr, "max= %g\n", max);

  TFile *hfile = new TFile(argv[1], "recreate");
  TH1D *h0 = new TH1D("h0", "signi.", 140, -7.0, 7.0);
  TH2D *h1 = new TH2D("h1", "on source events", fig_x_bin_num,
                      -(fig_x + bin_width / 2.0), fig_x + bin_width / 2.0,
                      fig_y_bin_num, -(fig_y + bin_width / 2.0),
                      fig_y + bin_width / 2.0);
  TH2D *h2 = new TH2D("h2", "off source events", fig_x_bin_num,
                      -(fig_x + bin_width / 2.0), fig_x + bin_width / 2.0,
                      fig_y_bin_num, -(fig_y + bin_width / 2.0),
                      fig_y + bin_width / 2.0);
  TH2D *h3 = new TH2D("h3", "on-off source events", fig_x_bin_num,
                      -(fig_x + bin_width / 2.0), fig_x + bin_width / 2.0,
                      fig_y_bin_num, -(fig_y + bin_width / 2.0),
                      fig_y + bin_width / 2.0);
  TH2D *h4 =
      new TH2D("h4", "on source", fig_x_bin_num, -(fig_x + bin_width / 2.0),
               fig_x + bin_width / 2.0, fig_y_bin_num,
               -(fig_y + bin_width / 2.0), fig_y + bin_width / 2.0);
  TH2D *h5 =
      new TH2D("h5", "off source", fig_x_bin_num, -(fig_x + bin_width / 2.0),
               fig_x + bin_width / 2.0, fig_y_bin_num,
               -(fig_y + bin_width / 2.0), fig_y + bin_width / 2.0);
  TH2D *h6 =
      new TH2D("h6", "on-off source", fig_x_bin_num, -(fig_x + bin_width / 2.0),
               fig_x + bin_width / 2.0, fig_y_bin_num,
               -(fig_y + bin_width / 2.0), fig_y + bin_width / 2.0);
  TH2D *h7 =
      new TH2D("h7", "signi.  ", fig_x_bin_num, -(fig_x + bin_width / 2.0),
               fig_x + bin_width / 2.0, fig_y_bin_num,
               -(fig_y + bin_width / 2.0), fig_y + bin_width / 2.0);
  TH1D *hra = new TH1D("hra", "on-off events", 11, -2.75, 2.75);
  TH1D *hra1 = new TH1D("hra1", "on-off events", fig_x_bin_num, -(fig_x + bin_width / 2.0), fig_x + bin_width / 2.0);
  TH1D *hdec = new TH1D("hdec", "on-off events", 11, -2.75, 2.75);
  TH1D *hdec1 = new TH1D("hdec1", "on-off events", fig_x_bin_num, -(fig_x + bin_width / 2.0), fig_x + bin_width / 2.0);

  for (i = (bin_num - fig_x_bin_num) / 2; i < (bin_num + fig_x_bin_num) / 2;
       i++) {
    for (j = (bin_num - fig_y_bin_num) / 2; j < (bin_num + fig_y_bin_num) / 2;
         j++) {
      //    for(i=40;i<161;i++) {
      //      for(j=40;j<161;j++) {
      if (off[i][j] == 0.)
        continue;
      sig = sig_smooth3(i, j, reso, &Non, &Noff);
      h0->Fill(sig);
      h1->Fill(i * bin_width - R_win, j * bin_width - R_win, on[i][j]);
      h2->Fill(i * bin_width - R_win, j * bin_width - R_win,
               off[i][j] / (2 * id_max));
      h3->Fill(i * bin_width - R_win, j * bin_width - R_win,
               on[i][j] - off[i][j] / (2 * id_max));
      h4->Fill(i * bin_width - R_win, j * bin_width - R_win, Non);
      h5->Fill(i * bin_width - R_win, j * bin_width - R_win,
               Noff / (2 * id_max));
      h6->Fill(i * bin_width - R_win, j * bin_width - R_win,
               Non - Noff / (2 * id_max));
      h7->Fill(i * bin_width - R_win, j * bin_width - R_win, sig);
      if (j >= (bin_num - 41) / 2 && j <= (bin_num + 41) / 2)
        hra->Fill(i * bin_width - R_win, on[i][j] - off[i][j] / (2 * id_max));
      hra1->Fill(i * bin_width - R_win, on[i][j] - off[i][j] / (2 * id_max));
      if (i >= (bin_num - 41) / 2 && i <= (bin_num + 41) / 2)
        hdec->Fill(j * bin_width - R_win, on[i][j] - off[i][j] / (2 * id_max));
        hdec1->Fill(j * bin_width - R_win, on[i][j] - off[i][j] / (2 * id_max));
      if (i == (bin_num - 1) / 2 && j == (bin_num - 1) / 2)
        printf("result %d %f %f %f %f %f\n", nk, reso, sig, Non, Noff, max);
    }
  }

  hfile->Write();
  hfile->Close();
  return 0;
}

double sig_smooth3(int x, int y, double r_smooth, double *Non, double *Noff) {
  int n_r;
  int i, j, i0, j0, i1, j1;
  double w = 0, r = 0;
  double dN = 0, dNoff = 0;
  double dis_square;
  double r_smooth_square = r_smooth * r_smooth;
  n_r = int(r_smooth / bin_width);

  *Non = 0.;
  *Noff = 0.;

  i0 = x - n_r;
  i1 = x + n_r;
  j0 = y - n_r;
  j1 = y + n_r;

  if (i0 < 0)
    i0 = 0;
  if (i1 >= bin_num)
    i1 = bin_num - 1;

  if (j0 < 0)
    j0 = 0;
  if (j1 >= bin_num)
    j1 = bin_num - 1;

  for (i = i0; i <= i1; i++) {
    for (j = j0; j <= j1; j++) {

      dis_square =
          ((i - x) * (i - x) + (j - y) * (j - y)) * bin_width * bin_width;
      if (dis_square > r_smooth_square)
        continue;

      if (nk > 12)
        w = 1;
      else
        w = gaus32(nk, sqrt(dis_square) * deg_rad);
      *Non += on[i][j] * w;
      *Noff += off[i][j] * w;
      dN += on[i][j] * w * w;
      dNoff += off[i][j] * w * w;
    }
  }

  if (dN != 0)
    return (*Non - *Noff / (2 * id_max)) /
           sqrt(dN + dNoff / (2 * 2 * id_max * id_max));
  else
    return -100;
}

double gaus32(int kk, double r) {
  int i;
  double result = 0;
  double reso2 = 0;

  for (i = 0; i < 3; i++) {
    reso2 = par[kk][i] * par[kk][i] * deg_rad * deg_rad;
    result += par[kk][i + 3] / (2 * pi * reso2) * exp(-r * r / (2 * reso2));
  }
  return result;
}
