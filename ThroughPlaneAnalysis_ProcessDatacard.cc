/*
==========================================
ThroughPlaneAnalysis_ProcessDatacard

Data analysis program to process datacard for
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
double quad(double a=0, double b=0, double c=0, double d=0, double e=0, double f=0, double g=0, double h=0, double i=0, double j=0, double k=0);

// Main program
int main()
{
  // Read data from Datacard.html
  ifstream ifs_datacard("Datacard.html");
  if (!ifs_datacard.good())
  {
    cout<<"ERROR: Datacard.html not found. Please run ThroughPlaneAnalysisReadRaw on the CSV file from the experiment first."<<endl;
    return 0;
  }
  string meta_experimentalist,
         meta_date,
         meta_material,
         meta_sample;
  vector<double> v_heaterRegionAvg, v_coolerRegionAvg,
                 v_heaterRegionErr, v_coolerRegionErr,
                 v_heaterRegionDiffAvg, v_coolerRegionDiffAvg,
                 v_heaterRegionDiffErr, v_coolerRegionDiffErr;
  string line;
  bool b_readData = true;
  while (b_readData)
  {
    vector<string> v_data;
    getline(ifs_datacard, line, '\n');
    splitLine(line, v_data, ',');
    if (v_data.at(0) == "Experimentalist")     meta_experimentalist = v_data.at(1);
    if (v_data.at(0) == "Date of measurement") meta_date = v_data.at(1);
    if (v_data.at(0) == "Material")            meta_material = v_data.at(1);
    if (v_data.at(0) == "Sample")              meta_sample = v_data.at(1);
    if (v_data.at(0) == "Hot Fluxmeter Temperatures (C)")
      for (unsigned int i = 1; i < v_data.size(); ++i)
        v_heaterRegionAvg.push_back(stod(v_data.at(i)));
    if (v_data.at(0) == "Cold Fluxmeter Temperatures (C)")
      for (unsigned int i = 1; i < v_data.size(); ++i)
        v_coolerRegionAvg.push_back(stod(v_data.at(i)));
    if (v_data.at(0) == "Hot Fluxmeter Temperature Uncertainties (C)")
      for (unsigned int i = 1; i < v_data.size(); ++i)
        v_heaterRegionErr.push_back(stod(v_data.at(i)));
    if (v_data.at(0) == "Cold Fluxmeter Temperature Uncertainties (C)")
      for (unsigned int i = 1; i < v_data.size(); ++i)
        v_coolerRegionErr.push_back(stod(v_data.at(i)));
    if (v_data.at(0) == "Hot Fluxmeter Temperatures Relative to Hottest Thermistor (C)")
      for (unsigned int i = 1; i < v_data.size(); ++i)
        v_heaterRegionDiffAvg.push_back(stod(v_data.at(i)));
    if (v_data.at(0) == "Cold Fluxmeter Temperatures Relative to Hottest Thermistor (C)")
      for (unsigned int i = 1; i < v_data.size(); ++i)
        v_coolerRegionDiffAvg.push_back(stod(v_data.at(i)));
    if (v_data.at(0) == "Hot Fluxmeter Temperature Relative to Hottest Thermistor Uncertainties (C)")
      for (unsigned int i = 1; i < v_data.size(); ++i)
        v_heaterRegionDiffErr.push_back(stod(v_data.at(i)));
    if (v_data.at(0) == "Cold Fluxmeter Temperature Relative to Hottest Thermistor Uncertainties (C)")
      for (unsigned int i = 1; i < v_data.size(); ++i)
        v_coolerRegionDiffErr.push_back(stod(v_data.at(i)));


    if (v_data.at(0) == "</pre>") b_readData = false;
  }
  ifs_datacard.close();

  // Read machine information from Apparatus.html
  ifstream ifs_apparatus("Apparatus.html");
  if (!ifs_apparatus.good())
  {
    cout<<"ERROR: Apparatus.html containing apparatus-specific information not found. Please softlink to it from the current directory. Aborting."<<endl;
    return 0;
  }
  string apparatus_name;
  double apparatus_thermistorDistance,
         apparatus_thermistorBoreDiameter,
         apparatus_fluxmeterDiameter,
         apparatus_fluxmeterConductivity;
  while (getline(ifs_apparatus, line, '\n'))
  {
    vector<string> v_data;
    splitLine(line, v_data, ':');
    if (v_data.at(0) == "Apparatus Name")                 apparatus_name = v_data.at(1);
    if (v_data.at(0) == "Inter-thermistor distance (mm)") apparatus_thermistorDistance = stod(v_data.at(1));
    if (v_data.at(0) == "Thermistor bore diameter (mm)")  apparatus_thermistorBoreDiameter = stod(v_data.at(1));
    if (v_data.at(0) == "Fluxmeter diameter (mm)")        apparatus_fluxmeterDiameter = stod(v_data.at(1));
    if (v_data.at(0) == "Fluxmeter conductivity (W/mK)")  apparatus_fluxmeterConductivity = stod(v_data.at(1));
  }
  ifs_apparatus.close();

  // Read bias information from BiasCalibration.html
  ifstream ifs_biasCalibration("BiasCalibration.html");
  if (!ifs_biasCalibration.good())
  {
    cout<<"ERROR: BiasCalibration.html containing bias calibration information not found. Please softlink to it from the current directory. Aborting."<<endl;
    return 0;
  }
  vector<double> v_heaterBias, v_coolerBias;
  while (getline(ifs_biasCalibration, line, '\n'))
  {
    vector<string> v_data;
    splitLine(line, v_data, ',');
    if (v_data.at(0) == "Heater Biases")
      for (unsigned int i = 1; i < v_data.size(); ++i)
        v_heaterBias.push_back(stod(v_data.at(i)));
    if (v_data.at(0) == "Cooler Biases")
      for (unsigned int i = 1; i < v_data.size(); ++i)
        v_coolerBias.push_back(stod(v_data.at(i)));
  }
  ifs_biasCalibration.close();
  // cout<<"v_heaterBias.size() = "<<v_heaterBias.size()<<endl;
  // cout<<"v_coolerBias.size() = "<<v_coolerBias.size()<<endl;

  const unsigned int n_heaterThermistors = v_heaterRegionAvg.size();
  const unsigned int n_coolerThermistors = v_coolerRegionAvg.size();
  // cout<<"n_heaterThermistors = "<<n_heaterThermistors<<endl;
  // cout<<"n_coolerThermistors = "<<n_coolerThermistors<<endl;

  cout<<"==========================================="<<endl;
  cout<<"Thermal Conductivity Through-Plane Analysis"<<endl;
  cout<<"Datacard processing program"<<endl;
  cout<<"Author: Souvik Das, 2023"<<endl;
  cout<<"=== Meta-data ============================="<<endl;
  cout<<"Enter name of analyst: ";
  std::string meta_analyst;
  getline(cin, meta_analyst);
  time_t t = time(0); tm* now = localtime(&t);
  string meta_analysisDate = to_string(now->tm_mon+1)+"/"+to_string(now->tm_mday)+"/"+to_string(now->tm_year+1900);
  cout<<"Analysis data: "<<meta_analysisDate<<endl;
  cout<<"Experimentalist: "<<meta_experimentalist<<endl;
  cout<<"Measurement Date: "<<meta_date<<endl;
  cout<<"Material: "<<meta_material<<endl;
  cout<<"Sample: "<<meta_sample<<endl;
  cout<<"Apparatus: "<<apparatus_name<<endl;

  // Calculate positions of thermistors in hot and cold fluxmeters
  // Add bias corrections to thermistor measurements
  // Calculate positions of the sample boundaries
  vector<double> v_heaterPosition, v_coolerPosition,
                 v_heaterPositionErr, v_coolerPositionErr;
  for (unsigned int i = 0; i < n_heaterThermistors; ++i)
  {
    v_heaterPosition.push_back(i*apparatus_thermistorDistance);
    v_heaterPositionErr.push_back(apparatus_thermistorBoreDiameter);
    v_heaterRegionAvg.at(i) += v_heaterBias.at(i);
    v_heaterRegionDiffAvg.at(i) += v_heaterBias.at(i);
  }
  for (unsigned int i = 0; i < n_coolerThermistors; ++i)
  {
    v_coolerPosition.push_back(i*apparatus_thermistorDistance);
    v_coolerPositionErr.push_back(apparatus_thermistorBoreDiameter);
    v_coolerRegionAvg.at(i) += v_coolerBias.at(i);
    v_coolerRegionDiffAvg.at(i) += (v_coolerBias.at(i) - v_coolerBias.at(0));
  }
  double x_sampleHotEnd = n_heaterThermistors * apparatus_thermistorDistance; // wrt hottest thermistor of hot fluxmeter
  double x_sampleColdEnd = -apparatus_thermistorDistance;                     // wrt hottest thermistor of cold fluxmeter
  cout<<"x_sampleHotEnd = "<<x_sampleHotEnd<<endl;
  cout<<"x_sampleColdEnd = "<<x_sampleColdEnd<<endl;

  // Fit the basic Heater Region plot
  TGraphErrors *g_HeaterFlux = new TGraphErrors(n_heaterThermistors, &v_heaterPosition[0], &v_heaterRegionAvg[0], &v_heaterPositionErr[0], &v_heaterRegionErr[0]);
  g_HeaterFlux->SetTitle("; Thermistor Position (mm); Average Thermistor Temperature (^{#circ}C)");
  g_HeaterFlux->GetXaxis()->SetLimits(-5., 45.);
  g_HeaterFlux->SetMarkerStyle(8);
  g_HeaterFlux->SetMarkerSize(1);
  TF1 *f_linear_heater = new TF1("f_linear_heater", "[0]+[1]*x", -2., 40.);
  f_linear_heater->SetParLimits(0, 5.1, 50.1);
  f_linear_heater->SetParLimits(1, -1., 0.);
  f_linear_heater->SetLineColor(kBlack);
  g_HeaterFlux->Fit(f_linear_heater, "QR+");
  TCanvas* c_HeaterFlux = new TCanvas();
  g_HeaterFlux->Draw("AP");
  c_HeaterFlux->SaveAs("c_HeaterFlux.png");
  c_HeaterFlux->SaveAs("c_HeaterFlux.pdf");

  // Fit the basic Cooler Region plot
  TGraphErrors *g_CoolerFlux = new TGraphErrors(n_coolerThermistors, &v_coolerPosition[0], &v_coolerRegionAvg[0], &v_coolerPositionErr[0], &v_coolerRegionErr[0]);
  g_CoolerFlux->SetTitle("; Thermistor Position (mm); Average Thermistor Temperature (^{#circ}C)");
  g_CoolerFlux->GetXaxis()->SetLimits(-5., 45.);
  g_CoolerFlux->SetMarkerStyle(8);
  g_CoolerFlux->SetMarkerSize(1);
  TF1 *f_linear_cooler = new TF1("f_linear_cooler", "[0]+[1]*x", -2., 40.);
  f_linear_cooler->SetParLimits(0, -5.1, 15.1);
  f_linear_cooler->SetParLimits(1, -0.10, 0.);
  f_linear_cooler->SetLineColor(kBlack);
  g_CoolerFlux->Fit(f_linear_cooler, "R+");
  TCanvas* c_CoolerFlux = new TCanvas();
  g_CoolerFlux->Draw("AP");
  c_CoolerFlux->SaveAs("c_CoolerFlux.png");
  c_CoolerFlux->SaveAs("c_CoolerFlux.pdf");

  // Fit the Heater Region Diff plot
  TGraphErrors *g_HeaterDiffFlux = new TGraphErrors(n_heaterThermistors, &v_heaterPosition[0], &v_heaterRegionDiffAvg[0], &v_heaterPositionErr[0], &v_heaterRegionDiffErr[0]);
  g_HeaterDiffFlux->SetTitle("; Thermistor Position (mm); Average Relative Temperature (^{#circ}C)");
  g_HeaterDiffFlux->GetXaxis()->SetLimits(-5., 45.);
  g_HeaterDiffFlux->SetMarkerStyle(8);
  g_HeaterDiffFlux->SetMarkerSize(1);
  TF1 *f_linear_heater_diff = new TF1("f_linear_heater_diff", "[0]+[1]*x", -2., 40.);
  f_linear_heater_diff->SetParLimits(0, -0.01, 0.01);
  f_linear_heater_diff->SetParLimits(1, -1., 0.);
  f_linear_heater_diff->SetLineColor(kBlack);
  g_HeaterDiffFlux->Fit(f_linear_heater_diff, "QR+");
  TCanvas* c_HeaterDiffFlux = new TCanvas();
  g_HeaterDiffFlux->Draw("AP");
  c_HeaterDiffFlux->SaveAs("c_HeaterDiffFlux.png");
  c_HeaterDiffFlux->SaveAs("c_HeaterDiffFlux.pdf");

  // Fit the basic Cooler Region Diff plot
  TGraphErrors *g_CoolerDiffFlux = new TGraphErrors(n_coolerThermistors, &v_coolerPosition[0], &v_coolerRegionDiffAvg[0], &v_coolerPositionErr[0], &v_coolerRegionDiffErr[0]);
  g_CoolerDiffFlux->SetTitle("; Thermistor Position (mm); Average Relative Temperature (^{#circ}C)");
  g_CoolerDiffFlux->GetXaxis()->SetLimits(-5., 45.);
  g_CoolerDiffFlux->SetMarkerStyle(8);
  g_CoolerDiffFlux->SetMarkerSize(1);
  TF1 *f_linear_cooler_diff = new TF1("f_linear_cooler_diff", "[0]+[1]*x", -2., 40.);
  f_linear_cooler_diff->SetParLimits(0, -0.01, 0.01);
  f_linear_cooler_diff->SetParLimits(1, -1., 0.);
  f_linear_cooler_diff->SetLineColor(kBlack);
  g_CoolerDiffFlux->Fit(f_linear_cooler_diff, "QR+");
  TCanvas* c_CoolerDiffFlux = new TCanvas();
  g_CoolerDiffFlux->Draw("AP");
  c_CoolerDiffFlux->SaveAs("c_CoolerDiffFlux.png");
  c_CoolerDiffFlux->SaveAs("c_CoolerDiffFlux.pdf");

  // Calculate results
  double pi = 3.14159265358979;
  double area = pi * pow(apparatus_fluxmeterDiameter/2000., 2);

  double j_in = -f_linear_heater_diff->GetParameter(1)*1000.*apparatus_fluxmeterConductivity;
  double j_in_err = f_linear_heater_diff->GetParError(1)*1000.*apparatus_fluxmeterConductivity;

  double j_out = -f_linear_cooler_diff->GetParameter(1)*1000.*apparatus_fluxmeterConductivity;
  double j_out_err = f_linear_cooler_diff->GetParError(1)*1000.*apparatus_fluxmeterConductivity;

  double j_avg = (j_in + j_out)/2.;
  double j_err = quad((j_in - j_out)/2., j_in_err/(2.*j_in), j_out_err/(2.*j_out));

  double i_avg = j_avg * area;
  double i_err = j_err * area;

  double T_hot = f_linear_heater->GetParameter(0) + f_linear_heater->GetParameter(1)*x_sampleHotEnd;
  double T_hot_err = quad(f_linear_heater->GetParError(0), f_linear_heater->GetParError(1)*x_sampleHotEnd);

  double T_cold = f_linear_cooler->GetParameter(0) + f_linear_cooler->GetParameter(1)*x_sampleColdEnd;
  double T_cold_err = quad(f_linear_cooler->GetParError(0), f_linear_cooler->GetParError(1)*x_sampleColdEnd);

  double deltaT = T_hot - T_cold;
  double deltaT_err = quad(T_hot_err, T_cold_err);

  double R = deltaT/i_avg;
  double R_err = R * quad((deltaT_err/deltaT), (i_err/i_avg));

  // Results which should be readable by CSV and HTML readers
  system("cp Datacard.html SampleResults.html");
  ofstream ofs_results("SampleResults.html", ios::app);
  ofs_results<<"<pre>"<<endl;
  ofs_results<<"Analyst, "<<meta_analyst<<endl;
  ofs_results<<"Analysis date, "<<meta_analysisDate<<endl;
  ofs_results<<"Apparatus, "<<apparatus_name<<endl;
  ofs_results<<"Apparatus inter-thermistor distance (mm), "<<apparatus_thermistorDistance<<endl;
  ofs_results<<"Apparatus thermistor bore diameter (mm), "<<apparatus_thermistorBoreDiameter<<endl;
  ofs_results<<"Apparatus fluxmeter diameter (mm), "<<apparatus_fluxmeterDiameter<<endl;
  ofs_results<<"Apparatus fluxmeter conductivity (W/mK), "<<apparatus_fluxmeterConductivity<<endl;
  ofs_results<<"Heater Biases: "<<v_heaterBias.at(0);
  for (unsigned int i = 1; i < v_heaterBias.size(); ++i)
    ofs_results<<", "<<v_heaterBias.at(i);
  ofs_results<<endl;
  ofs_results<<"Cooler Biases: "<<v_coolerBias.at(0);
  for (unsigned int i = 1; i < v_coolerBias.size(); ++i)
    ofs_results<<", "<<v_coolerBias.at(i);
  ofs_results<<endl;
  ofs_results<<"Hot fluxmeter flux (W/m^2), "<<j_in<<" +/- "<<j_in_err<<endl;
  ofs_results<<"Cold fluxmeter flux (W/m^2), "<<j_out<<" +/- "<<j_out_err<<endl;
  ofs_results<<"Heat flux through sample (W/m^2), "<<j_avg<<" +/- "<<j_err<<endl;
  ofs_results<<"Heat flow through sample (W), "<<i_avg<<" +/- "<<i_err<<endl;
  ofs_results<<"Temperature at hot end of sample (C), "<<T_hot<<" +/- "<<T_hot_err<<endl;
  ofs_results<<"Temperature at cold end of sample (C), "<<T_cold<<" +/- "<<T_cold_err<<endl;
  ofs_results<<"deltaT (C), "<<deltaT<<" +/- "<<deltaT_err<<endl;
  ofs_results<<"R (K/W), "<<R<<", "<<R_err<<endl;
  ofs_results<<"</pre>"<<endl;
  ofs_results<<"<table border='1'>"<<endl;
  ofs_results<<" <tr>"<<endl;
  ofs_results<<"  <td>"<<endl;
  ofs_results<<"   Hot flux meter fit <br/>"<<endl;
  ofs_results<<"   <img src='c_HeaterFlux.png'/>"<<endl;
  ofs_results<<"  </td>"<<endl;
  ofs_results<<"  <td>"<<endl;
  ofs_results<<"   Cold flux meter fit <br/>"<<endl;
  ofs_results<<"   <img src='c_CoolerFlux.png'/> <br/>"<<endl;
  ofs_results<<"  </td>"<<endl;
  ofs_results<<" </tr>"<<endl;
  ofs_results<<" <tr>"<<endl;
  ofs_results<<"  <td>"<<endl;
  ofs_results<<"   Hot flux meter difference fit <br/>"<<endl;
  ofs_results<<"   <img src='c_HeaterDiffFlux.png'/>"<<endl;
  ofs_results<<"  </td>"<<endl;
  ofs_results<<"  <td>"<<endl;
  ofs_results<<"   Cold flux meter difference fit <br/>"<<endl;
  ofs_results<<"   <img src='c_CoolerDiffFlux.png'/> <br/>"<<endl;
  ofs_results<<"  </td>"<<endl;
  ofs_results<<" </tr>"<<endl;
  ofs_results<<" <tr>"<<endl;
  ofs_results<<"</table>"<<endl;
  ofs_results.close();
  cout<<"Results in SampleResults.html"<<endl;
  system("open SampleResults.html");
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
