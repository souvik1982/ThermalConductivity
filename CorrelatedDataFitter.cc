#include <iostream>
#include <vector>
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TH1D.h"
#include "TCanvas.h"

using namespace std;

unsigned int global_n;
vector<double> global_v_x, global_v_xErr, global_v_y, global_v_yErr;
Double_t* global_array_Cij_inv;
Double_t fitFunction(double x, Double_t *par);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

double quad(double a=0, double b=0, double c=0, double d=0, double e=0, double f=0, double g=0, double h=0, double i=0, double j=0, double k=0);

class CorrelatedDataFitter
{
public:
  CorrelatedDataFitter(double correlation, vector<double> v_x, vector<double> v_y, vector<double> v_xErr, vector<double> v_yErr);
  ~CorrelatedDataFitter();
  TCanvas* getMinuitFit(string title, double &m, double &dm, double &c, double &dc);
  TGraphErrors* getAnalyticalFit(double &m, double &dm, double &c, double &dc);
};

CorrelatedDataFitter::CorrelatedDataFitter(double correlation, vector<double> v_x, vector<double> v_y, vector<double> v_xErr, vector<double> v_yErr)
{
  global_n = v_x.size();
  global_v_x = v_x;
  global_v_y = v_y;
  global_v_xErr = v_xErr;
  global_v_yErr = v_yErr;

  // Add arbitrary correlations between data points before fitting them.
  // First, construct the covariance matrix Cij and invert it.
  Double_t array_Cij[global_n*global_n];
  cout<<"Cij = "<<endl;
  for (unsigned int i=0; i<global_n; ++i)
  {
    for (unsigned int j=0; j<global_n; ++j)
    {
      if (i == j)
        array_Cij[i*global_n+j] = pow(v_yErr.at(i), 2);
      else
        array_Cij[i*global_n+j] = correlation*v_yErr.at(i)*v_yErr.at(j);
      cout<<array_Cij[i*global_n+j]<<" ";
    }
    cout<<endl;
  }
  cout<<"Inv(Cij) =  "<<endl;
  TMatrixD Cij_inv(global_n, global_n);
  Cij_inv.SetMatrixArray(array_Cij);
  Cij_inv.Invert();
  global_array_Cij_inv = new Double_t[global_n*global_n+1];
  memcpy(global_array_Cij_inv, Cij_inv.GetMatrixArray(), global_n*global_n*sizeof(Double_t));
  for (unsigned int i = 0; i<global_n; ++i)
  {
    for (unsigned int j = 0; j<global_n; ++j)
      cout<<global_array_Cij_inv[i*global_n+j]<<" ";
    cout<<endl;
  }
}

TGraphErrors* CorrelatedDataFitter::getAnalyticalFit(double &m, double &dm, double &c, double &dc)
{
  double A=0, B=0, C=0, D=0, E=0;
  double dE=0, dC=0;

  for (unsigned int i=0; i<global_n; ++i)
  {
    for (unsigned int j=0; j<global_n; ++j)
    {
      A += global_v_x.at(i) * global_v_x.at(j) * global_array_Cij_inv[i*global_n+j];
      B += global_v_x.at(i) * global_array_Cij_inv[i*global_n+j];
      C += global_v_x.at(i) * global_v_y.at(j) * global_array_Cij_inv[i*global_n+j];
      D += global_array_Cij_inv[i*global_n+j];
      E += global_v_y.at(j) * global_array_Cij_inv[i*global_n+j];

      // For uncertainties
      dE += global_v_yErr.at(i) * global_array_Cij_inv[i*global_n+j];
      dC += global_v_x.at(i) * global_v_yErr.at(j) * global_array_Cij_inv[i*global_n+j];
    }
  }

  m = (E*B - C*D)/(B*B - A*D);
  c = (C*B - E*A)/(B*B - A*D);
  dm = 2.*(B*dE - D*dC)/(B*B - A*D);

  TF1 *f_linear = new TF1("f_linear", "[0] + [1]*x", global_v_x.at(0) - 2, global_v_x.at(global_n - 1) + 2);
  f_linear->SetParameter(0, c); f_linear->FixParameter(0, c);
  f_linear->SetParameter(1, m); f_linear->FixParameter(1, m);
  f_linear->SetLineColor(1);

  TGraphErrors *g = new TGraphErrors(global_n, &global_v_x[0], &global_v_y[0], &global_v_xErr[0], &global_v_yErr[0]);
  g->Fit(f_linear, "R");

  return g;
}

TCanvas* CorrelatedDataFitter::getMinuitFit(string title, double &slope, double &slopeErr, double &intercept, double &interceptErr)
{
  TMinuit *gMinuit = new TMinuit(2);
  gMinuit->SetFCN(fcn);

  Double_t arglist[10];
  Int_t ierflag = 0;

  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflag);

  // Set starting values and step sizes for parameters
  double par_init[2] = {5e-5, 5e-7};
  double par_step[2] = {1e-6, 1e-8};
  gMinuit->mnparm(0, "intercept", par_init[0], par_step[0], 0, 0, ierflag);
  gMinuit->mnparm(1, "slope", par_init[1], par_step[1], 0, 0, ierflag);

  // Minimize
  arglist[0] = 500;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflag);

  // Get the straight line fit
  TString string_0, string_1;
  Double_t xlo_0, xlo_1;
  Double_t xhi_0, xhi_1;
  Int_t i_0, i_1;
  gMinuit->mnpout(0, string_0, intercept, interceptErr, xlo_0, xhi_0, i_0);
  gMinuit->mnpout(1, string_1, slope, slopeErr, xlo_1, xhi_1, i_1);
  cout<<"Minuit fit, slope = "<<slope<<", slopeErr = "<<slopeErr<<endl;
  cout<<"Minuit fit, intercept = "<<intercept<<", interceptErr = "<<interceptErr<<endl;
  TF1 *f_linear = new TF1("f_linear", "[0] + [1]*x", global_v_x.at(0) - 2, global_v_x.at(global_n - 1) + 2);
  f_linear->SetParameter(0, intercept); f_linear->FixParameter(0, intercept);
  f_linear->SetParameter(1, slope); f_linear->FixParameter(1, slope);
  f_linear->SetLineColor(1);
  TGraphErrors *g = new TGraphErrors(global_n, &global_v_x[0], &global_v_y[0], &global_v_xErr[0], &global_v_yErr[0]);
  g->SetTitle(title.c_str());
  g->SetMarkerStyle(8);
  g->Fit(f_linear, "R");

  // Compute the band
  Double_t *errorMatrix = new Double_t[4];
  gMinuit->mnemat(errorMatrix, 2);
  cout<<"errorMatrix[0] = "<<errorMatrix[0]<<endl;
  cout<<"errorMatrix[1] = "<<errorMatrix[1]<<endl;
  cout<<"errorMatrix[2] = "<<errorMatrix[2]<<endl;
  cout<<"errorMatrix[3] = "<<errorMatrix[3]<<endl;
  Double_t correlationParams = errorMatrix[1]/sqrt(errorMatrix[0]*errorMatrix[3]);
  cout<<"correlationParams = "<<correlationParams<<endl;
  double N = 100;
  TH1D *h_band = new TH1D("h_band", "h_band", N, global_v_x[0], global_v_x[global_n-1]);
  for (int i=0; i<N; ++i)
  {
    double x = ((i+0.5)*(global_v_x[global_n-1]-global_v_x[0])/N + global_v_x[0]);
    double y = slope*x + intercept;
    double dy = 2.*sqrt(pow(x*slopeErr, 2) + pow(interceptErr, 2) + 2.*correlationParams*x*slopeErr*interceptErr);
    h_band->SetBinContent(i+1, y);
    h_band->SetBinError(i+1, dy);
  }
  h_band->SetFillColorAlpha(1, 0.3); h_band->SetFillStyle(3209);

  // Put the plot together with data points
  TCanvas *c = new TCanvas();
  g->Draw("AP");
  h_band->Draw("E3 SAME");

  return c;
}

CorrelatedDataFitter::~CorrelatedDataFitter()
{
  // delete[] global_array_Cij_inv;
}

Double_t fitFunction(Double_t x, Double_t *par)
{
  Double_t value = par[0] + par[1] * x;
  return value;
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  double chi2 = 0;
  for (unsigned int i=0; i<global_n; ++i)
  {
    for (unsigned int j=0; j<global_n; ++j)
    {
      chi2 += (global_v_y.at(i) - fitFunction(global_v_x.at(i), par))*(global_v_y.at(j) - fitFunction(global_v_x.at(j), par))*global_array_Cij_inv[i*global_n+j];
    }
  }
  f = chi2;
}

double quad(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j, double k){
  return pow(a*a+b*b+c*c+d*d+e*e+f*f+g*g+h*h+i*i+j*j+k*k, 0.5);
}
