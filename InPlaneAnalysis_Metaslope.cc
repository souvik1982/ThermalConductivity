/*
==========================================
InPlaneAnalysis_Metaslope

Data analysis program to process MetaSlope for
Thermal Conductivity measurement with In-plane Apparatus

April 2023, Souvik Das
==========================================
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <ctime>
#include "TCanvas.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TVirtualFitter.h"

#include "CorrelatedDataFitter.cc"

using namespace std;

double corrUnc = 0.5; // Correlation of uncertainties between data points.

// Forward declared functions
void splitLine(string &line, vector<string> &words, char separator); // Splits a line into a vector of strings
double quad(double a=0, double b=0, double c=0, double d=0, double e=0, double f=0, double g=0, double h=0, double i=0, double j=0, double k=0);
// extern void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

int main()
{
  cout<<"==========================================="<<endl;
  cout<<"Thermal Conductivity In-Plane Analysis"<<endl;
  cout<<"Metaslope program"<<endl;
  cout<<"Author: Souvik Das, 2023"<<endl;
  cout<<"=== Meta-data ============================="<<endl;
  cout<<"Enter name of analyst: ";
  std::string meta_analyst;
  getline(cin, meta_analyst);
  time_t t = time(0); tm* now = localtime(&t);
  string meta_analysisDate = to_string(now->tm_mon+1)+"/"+to_string(now->tm_mday)+"/"+to_string(now->tm_year+1900);
  cout<<"Analysis data: "<<meta_analysisDate<<endl;
  cout<<"Enter material codename: ";
  std::string meta_material;
  getline(cin, meta_material);
  cout<<"Which sample numbers should I analyze? These folders must exist with completed and inspected analyses resulting in SampleResults.html within them. \n Enter integers separated by commas: ";
  std::string line;
  getline(cin, line);
  vector<string> v_samples;
  splitLine(line, v_samples, ',');

  // Search for SampleResults.html in subfolders containing samples of this material analyzed
  // Fill the vectors for the metaslope graph
  vector<double> v_length, v_length_err,
                 v_width, v_width_err,
                 v_thickness, v_thickness_err,
                 v_RA, v_RA_err;
  for (unsigned int i = 0; i < v_samples.size(); ++i)
  {
    string filename = meta_material+"_"+v_samples.at(i)+"/SampleResults.html";
    ifstream ifs_sampleResults(filename);
    if (ifs_sampleResults.good())
    {
      string line;
      while (getline(ifs_sampleResults, line))
      {
        vector<string> v_words;
        splitLine(line, v_words, ',');
        string key = v_words.at(0);
        if (key == "Length (mm)")        {v_length.push_back(stod(v_words.at(1))); v_length_err.push_back(stod(v_words.at(2)));}
        else if (key == "RA (Km^2/W)")   {v_RA.push_back(stod(v_words.at(1))); v_RA_err.push_back(stod(v_words.at(2)));}
      }
      cout<<"LOG: Processed "<<filename<<endl;
    }
    else
    {
      std::cout<<"ERROR: "<<filename<<" does not exist. Skipping it.";
    }
  }

  TF1 *f_linear = new TF1("f_linear", "[0] + [1]*x");
  f_linear->SetLineColor(1);

  TGraphErrors *g_metaSlope = new TGraphErrors(v_length.size(), &v_length[0], &v_RA[0], &v_length_err[0], &v_RA_err[0]);
  g_metaSlope->SetTitle("; Sample length (#mm); Thermal insulance (Km^{2}/W)");

  g_metaSlope->Fit(f_linear);
  TGraphErrors *g_metaSlope_confidence = new TGraphErrors(v_length.size(), &v_length[0], &v_RA[0]);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(g_metaSlope_confidence);
  g_metaSlope_confidence->SetFillColorAlpha(1, 0.3); g_metaSlope_confidence->SetFillStyle(3209);

  TCanvas *c_metaSlope = new TCanvas();
  g_metaSlope_confidence->Draw("AP3");
  g_metaSlope->Draw("PSAME");
  c_metaSlope->SaveAs(("c_metaSlope_"+meta_material+".png").c_str());
  c_metaSlope->SaveAs(("c_metaSlope_"+meta_material+".pdf").c_str());

  double slope = f_linear->GetParameter(1);
  double slope_err = f_linear->GetParError(1);
  double intercept = f_linear->GetParameter(0);
  double intercept_err = f_linear->GetParError(0);

  double k = 1./(slope*1e3);
  double k_err = k * quad(slope_err/slope);

  ofstream ofs_result("Result.html");
  ofs_result<<"<pre>"<<endl;
  ofs_result<<"Through-plane thermal conductivity measurement"<<endl;
  ofs_result<<"Analyst: "<<meta_analyst<<endl;
  ofs_result<<"Analysis date: "<<meta_analysisDate<<endl;
  ofs_result<<"Material: "<<meta_material<<endl;
  ofs_result<<"k = "<<k<<" +/- "<<k_err<<" W/mK"<<endl;
  ofs_result<<"</pre>"<<endl;
  ofs_result<<"Underlying sample results: <br/>"<<endl;
  for (unsigned int i = 0; i < v_samples.size(); ++i)
  {
    string href = meta_material+"_"+v_samples.at(i)+"/SampleResults.html";
    ofs_result<<"<a href = '"<<href<<"' target='_blank'>"<<meta_material<<"_"+v_samples.at(i)<<"</a> <br/>"<<endl;
  }
  ofs_result<<"<img src='c_metaSlope_"<<meta_material<<".png'/> <br/>"<<endl;
  ofs_result.close();
  cout<<"Result in Result.html"<<endl;
  system("open Result.html");

  // Test inversion of correlation function
  CorrelatedDataFitter correlatedDataFitter(0.0, v_length, v_RA, v_length_err, v_RA_err);

}

string removeLeadingSpaces(string &s)
{
  size_t start = s.find_first_not_of(' ');
  return (start == std::string::npos) ? "" : s.substr(start);
}

void splitLine(string &line, vector<string> &words, char delimiter)
{
  stringstream ss(line);
  string word;
  while(ss.good())
  {
    getline(ss, word, delimiter);
    words.push_back(removeLeadingSpaces(word));
  }
}

double quad(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j, double k){
  return pow(a*a+b*b+c*c+d*d+e*e+f*f+g*g+h*h+i*i+j*j+k*k, 0.5);
}
