#include <iostream>
#include <vector>
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TMatrixD.h"

using namespace std;

unsigned int global_n;
vector<double> global_v_x, global_v_xErr, global_v_y, global_v_yErr;
Double_t *array_Cij_inv;
Double_t *global_par;
Double_t fitFunction(double x, Double_t *par);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

class CorrelatedDataFitter
{
public:
  TGraphErrors *g_;

  CorrelatedDataFitter(double correlation, vector<double> v_x, vector<double> v_y, vector<double> v_xErr, vector<double> v_yErr);
  TGraphErrors* getFit();
};

CorrelatedDataFitter::CorrelatedDataFitter(double correlation, vector<double> v_x, vector<double> v_y, vector<double> v_xErr, vector<double> v_yErr)
{
  global_n = v_x.size();
  global_v_x = v_x;
  global_v_y = v_y;
  global_v_xErr = v_xErr;
  global_v_yErr = v_yErr;

  g_ = new TGraphErrors(global_n, &global_v_x[0], &global_v_y[0], &global_v_xErr[0], &global_v_yErr[0]);

  // Add arbitrary correlations between data points before fitting them.
  // First, construct the covariance matrix Cij and invert it.
  double array_Cij[global_n*global_n];
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
  array_Cij_inv = Cij_inv.GetMatrixArray();
  for (unsigned int i = 0; i<global_n; ++i)
  {
    for (unsigned int j = 0; j<global_n; ++j)
      cout<<array_Cij_inv[i*global_n+j]<<" ";
    cout<<endl;
  }
}

TGraphErrors* CorrelatedDataFitter::getFit()
{
  // TMinuit *gMinuit = new TMinuit();
  // gMinuit->SetFCN(fcn);

  return g_;
}

Double_t fitFunction(double x, Double_t *par)
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
      chi2 += (global_v_y.at(i) - fitFunction(global_v_x.at(i), global_par))*(global_v_y.at(j) - fitFunction(global_v_x.at(j), global_par))*array_Cij_inv[i*global_n+j];
    }
  }
  f = chi2;
}
