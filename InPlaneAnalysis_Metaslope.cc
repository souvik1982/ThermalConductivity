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
#include "TMinuit.h"

using namespace std;

double correlation = 0.5; // Correlation of uncertainties between data points.

// Forward declared functions
void splitLine(string &line, vector<string> &words, char separator); // Splits a line into a vector of strings

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

  // Read FEAHeatLoss.html file and incorporate uncertainties into v_RA_err
  ifstream ifs_FEA;
  ifs_FEA.open("FEAHeatLoss.html");

  // Start writing Result.html
  ofstream ofs_result("Result.html");
  ofs_result<<"<pre>"<<endl;
  ofs_result<<"====================================================="<<endl;
  ofs_result<<"In-plane thermal conductivity measurement"<<endl;
  ofs_result<<"Output of InPlaneAnalysis_Metaslope. Souvik Das, 2023"<<endl;
  ofs_result<<"====================================================="<<endl;
  ofs_result<<"Analyst: "<<meta_analyst<<endl;
  ofs_result<<"Analysis date: "<<meta_analysisDate<<endl;
  ofs_result<<"Material: "<<meta_material<<endl;
  ofs_result<<"</pre>"<<endl;
  ofs_result<<"<hr/>"<<endl;
  ofs_result<<"Underlying sample results: <br/>"<<endl;

  // Search for SampleResults.html in subfolders containing samples of this material analyzed
  // Fill the vectors for the metaslope graph
  // Add FEA heat loss information if available
  vector<Double_t> v_length, v_length_err,
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
      ofs_result<<"<a href = '"<<filename<<"' target='_blank'>"<<meta_material<<"_"+v_samples.at(i)<<"</a>";
    }
    else
      std::cout<<"ERROR: "<<filename<<" does not exist. Skipping it.";

    while (getline(ifs_FEA, line))
    {
      vector<string> v_words;
      splitLine(line, v_words, ',');
      string key = v_words.at(0);
      if (key == meta_material+"_"+v_samples.at(i))
      {
        double error_FEA = stod(v_words.at(1));
        v_RA_err.at(i) += v_RA.at(i) * error_FEA;
        cout<<"LOG: FEA Heat Loss uncertainty. Added fractional uncertainty of "<<error_FEA<<" to "<<key<<endl;
        ofs_result<<", FEA Heat Loss uncertainty = "<<error_FEA;
      }
    }
    ifs_FEA.clear();
    ifs_FEA.seekg(0, ios::beg);
    ofs_result<<"</br>"<<endl;
  }

  ofs_result<<"<hr/>"<<endl;

  // Naive fit with no correlation between data points
  TF1 *f_linear = new TF1("f_linear", "[0] + [1]*x", v_length.at(0) - 2, v_length.at(v_length.size() - 1) + 2);
  f_linear->SetLineColor(1);

  TGraphErrors *g_metaslope = new TGraphErrors(v_length.size(), &v_length[0], &v_RA[0], &v_length_err[0], &v_RA_err[0]);
  g_metaslope->SetMarkerStyle(8);
  g_metaslope->Fit(f_linear, "R");
  TGraphErrors *g_metaslope_confidence = new TGraphErrors(v_length.size(), &v_length[0], &v_RA[0]);
  g_metaslope_confidence->SetTitle("; Sample length (mm); Thermal insulance (Km^{2}/W)");
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(g_metaslope_confidence);
  g_metaslope_confidence->SetFillColorAlpha(1, 0.3); g_metaslope_confidence->SetFillStyle(3209);

  TCanvas *c_metaSlope = new TCanvas();
  g_metaslope_confidence->Draw("AP3");
  g_metaslope->Draw("PSAME");
  c_metaSlope->SaveAs(("c_metaSlope_"+meta_material+".png").c_str());
  c_metaSlope->SaveAs(("c_metaSlope_"+meta_material+".pdf").c_str());

  double slope = f_linear->GetParameter(1);
  double slope_err = f_linear->GetParError(1);
  double intercept = f_linear->GetParameter(0);
  double intercept_err = f_linear->GetParError(0);

  double k = 1./(slope*1e3);
  double k_err = k * slope_err/slope;
  ofs_result<<"<img src='c_metaSlope_"<<meta_material<<".png'/> <br/>"<<endl;
  ofs_result<<"Iterative fit with no correlations between data points. k = "<<k<<" +/- "<<k_err<<" W/mK"<<endl;
  ofs_result<<"<hr/>"<<endl;

  // Analytical fit with correlation
  CorrelatedDataFitter correlatedDataFitter(correlation, v_length, v_RA, v_length_err, v_RA_err);
  double slope_analytical, slope_err_analytical, intercept_analytical, intercept_err_analytical;
  TGraphErrors *g_metaslope_analytical = correlatedDataFitter.getAnalyticalFit(slope_analytical, slope_err_analytical, intercept_analytical, intercept_err_analytical);
  g_metaslope_analytical->SetTitle("; Sample length (mm); Thermal insulance (Km^{2}/W)");
  g_metaslope_analytical->SetMarkerStyle(8);

  double k_analytical = 1./(slope_analytical*1e3);
  double k_err_analytical = k_analytical * slope_err_analytical/slope_analytical;

  TCanvas *c_metaslope_analytical = new TCanvas();
  g_metaslope_analytical->Draw("AP");
  c_metaslope_analytical->SaveAs(("c_metaSlope_analytical_"+meta_material+".png").c_str());
  ofs_result<<"<img src='c_metaSlope_analytical_"<<meta_material<<".png'/> <br/>"<<endl;
  ofs_result<<"Analytical fit with "<<correlation<<" correlation between data point uncertainties. k = "<<k_analytical<<" +/- "<<k_err_analytical<<" W/mK"<<endl;
  ofs_result<<"<hr/>"<<endl;

  // Minuit Fit with correlation
  double slope_minuit, slope_err_minuit, intercept_minuit, intercept_err_minuit;
  TCanvas *c_metaSlope_minuit = correlatedDataFitter.getMinuitFit("; Sample length (mm); Thermal insulance (Km^{2}/W)", slope_minuit, slope_err_minuit, intercept_minuit, intercept_err_minuit);

  double k_minuit = 1./(slope_minuit*1e3);
  double k_err_minuit = k_minuit * slope_err_minuit / slope_minuit;

  c_metaSlope_minuit->SaveAs(("c_metaSlope_minuit_"+meta_material+".png").c_str());
  c_metaSlope_minuit->SaveAs(("c_metaSlope_minuit_"+meta_material+".pdf").c_str());
  ofs_result<<"<img src='c_metaSlope_minuit_"<<meta_material<<".png'/> <br/>"<<endl;
  ofs_result<<"Iterative Minuit fit with "<<correlation<<" correlation between data point uncertainties. k = "<<k_minuit<<" +/- "<<k_err_minuit<<" W/mK"<<endl;

  ofs_result.close();
  cout<<"Result in Result.html"<<endl;
  system("open Result.html");
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
