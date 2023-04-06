/*
==========================================
ThroughPlaneAnalysis_Metaslope

Data analysis program to process MetaSlope for
Thermal Conductivity measurement with Through-plane Apparatus

January 2023, Souvik Das
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

using namespace std;

// Forward declared functions
void splitLine(string &line, vector<string> &words, char separator); // Splits a line into a vector of strings

int main()
{
  cout<<"==========================================="<<endl;
  cout<<"Thermal Conductivity Through-Plane Analysis"<<endl;
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
  double pi = 3.14159265358979;
  double sample_diameter;
  vector<double> v_thickness,
                 v_thickness_err,
                 v_R,
                 v_R_err;
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
        if (key == "Thickness (um)") {v_thickness.push_back(stod(v_words.at(1))); v_thickness_err.push_back(stod(v_words.at(2)));}
        else if (key == "R (K/W)")   {v_R.push_back(stod(v_words.at(1))); v_R_err.push_back(stod(v_words.at(2)));}
        else if (key == "Apparatus fluxmeter diameter (mm)") sample_diameter = stod(v_words.at(1));
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

  TGraphErrors *g_metaSlope = new TGraphErrors(v_thickness.size(), &v_thickness[0], &v_R[0], &v_thickness_err[0], &v_R_err[0]);
  g_metaSlope->SetTitle("; Sample thickness (#mum); Thermal resistance (K/W)");
  g_metaSlope->Fit(f_linear, "Q");

  TCanvas *c_metaSlope = new TCanvas();
  g_metaSlope->Draw("AP");
  c_metaSlope->SaveAs(("c_metaSlope_"+meta_material+".png").c_str());
  c_metaSlope->SaveAs(("c_metaSlope_"+meta_material+".pdf").c_str());

  double slope = f_linear->GetParameter(1);
  double slope_err = f_linear->GetParError(1);
  double intercept = f_linear->GetParameter(0);
  double intercept_err = f_linear->GetParError(0);

  double sample_area = pi * pow(sample_diameter/2000., 2);

  double k = 1./(slope*sample_area*1e6);
  double k_err = k * slope_err/slope;

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
