/*
==========================================
ThroughPlaneAnalysis_ReadRaw

Data analysis program to create datacard for
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
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TArrow.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"

using namespace std;

const int minIntegrationTime = 1800; // seconds

// Forward declared functions
void splitLine(string &line, vector<string> &words, char separator); // Splits a line into a vector of strings
void selectRegion(vector<vector<double>> &heater_sample, vector<vector<double>> &peltier_sample, unsigned int time_size, unsigned int &min, unsigned int &max, unsigned int &size); // Select integration window with least combined standard deviation

// Main program
int main(int argc, char **argv)
{
  if (argc != 2)
  {
    cout<<"ERROR: Filename of the CSV file containing meta and raw data was not given. Aborting."<<endl;
    return 0;
  }
  std::string filename = argv[1];

  // Read and enforce metadata in raw data
  ifstream ifs_rawData(filename.c_str());
  string meta_experimentalist,
         meta_date,
         meta_apparatus,
         meta_material,
         meta_sample,
         meta_thickness,
         meta_thicknessErr,
         meta_heaterVoltage,
         meta_coolerVoltage;
  string line;
  bool b_readMeta = true;
  while (b_readMeta)
  {
    vector<string> v_metaData;
    getline(ifs_rawData, line, '\n');
    splitLine(line, v_metaData, ',');
    string key = v_metaData.at(0);
         if (key == "Experimentalist") meta_experimentalist = v_metaData.at(1);
    else if (key == "Date")            meta_date = v_metaData.at(1);
    else if (key == "Apparatus")       meta_apparatus = v_metaData.at(1);
    else if (key == "Material")        meta_material = v_metaData.at(1);
    else if (key == "Sample")          meta_sample = v_metaData.at(1);
    else if (key == "Thickness (um)")  {meta_thickness = v_metaData.at(1); meta_thicknessErr = v_metaData.at(2);}
    else if (key == "Heater Voltage (V)")  meta_heaterVoltage = v_metaData.at(1);
    else if (key == "Cooler Voltage (V)")  meta_coolerVoltage = v_metaData.at(1);
    else if (key == "Item")            b_readMeta = false;
  }
  cout<<"==========================================="<<endl;
  cout<<"Thermal Conductivity Through-Plane Analysis"<<endl;
  cout<<"Data ingestion program"<<endl;
  cout<<"Author: Souvik Das, 2023"<<endl;
  cout<<"=== Meta-data ============================="<<endl;
  cout<<"Experimentalist: "<<meta_experimentalist<<endl;
  cout<<"Date: "<<meta_date<<endl;
  cout<<"Apparatus: "<<meta_apparatus<<endl;
  cout<<"Material: "<<meta_material<<endl;
  cout<<"Sample: "<<meta_sample<<endl;
  cout<<"Heater Voltage: "<<meta_heaterVoltage<<" V"<<endl;
  cout<<"Cooler Voltage: "<<meta_coolerVoltage<<" V"<<endl;
  // End reading metadata

  // Read machine information from Apparatus.html
  string env_thermalBase = string(getenv("THERMALBASE"));
  if (env_thermalBase == "")
  {
    cout<<"ERROR: Environmental variable THERMALBASE not set. Aborting."<<endl;
    return 0;
  }
  ifstream ifs_apparatus((env_thermalBase+"/Apparatus.html").c_str());
  if (!ifs_apparatus.good())
  {
    cout<<"ERROR: Apparatus.html containing apparatus-specific information not found. Please softlink to it from the current directory. Aborting."<<endl;
    return 0;
  }
  vector<string> v_heaterChannel, v_coolerChannel;
  string line_apparatus;
  while (getline(ifs_apparatus, line_apparatus, '\n'))
  while (getline(ifs_apparatus, line_apparatus, '\n'))
  {
    vector<string> v_data;
    splitLine(line_apparatus, v_data, ',');
    if (v_data.at(0) == "<pre>")
    {
      getline(ifs_apparatus, line_apparatus, '\n');
      v_data.clear();
      splitLine(line_apparatus, v_data, ',');
      if (v_data.at(0) == "Apparatus Name" && v_data.at(1) == meta_apparatus)
      {
        while (v_data.at(0) != "</pre>")
        {
          getline(ifs_apparatus, line_apparatus, '\n');
          v_data.clear();
          splitLine(line_apparatus, v_data, ',');
          if (v_data.at(0) == "Heater Channels")
          {
            for (unsigned int i = 1; i < v_data.size(); ++i)
            {
              string content = v_data.at(i);
              if (content.size() == 3) v_heaterChannel.push_back(v_data.at(i));
            }
          }
          else if (v_data.at(0) == "Cooler Channels")
          {
            for (unsigned int i = 1; i < v_data.size(); ++i)
            {
              string content = v_data.at(i);
              if (content.size() == 3) v_coolerChannel.push_back(v_data.at(i));
            }
          }
        }
      }
    }
  }
  const unsigned int n_heaterThermistors = v_heaterChannel.size();
  const unsigned int n_coolerThermistors = v_coolerChannel.size();
  cout<<" Heater Channels: ";
  for (unsigned int i = 0; i < v_heaterChannel.size(); ++i)
    cout<<v_heaterChannel.at(i)<<" ";
  cout<<endl;
  cout<<" Cooler Channels: ";
  for (unsigned int i = 0; i < v_coolerChannel.size(); ++i)
    cout<<v_coolerChannel.at(i)<<" ";
  cout<<endl;
  cout<<"If any of these entries are missing of incorrect, please alert the experimentalist and report the incident to the supervisor. - Souvik Das"<<endl;
  cout<<"==========================================="<<endl;
  ifs_apparatus.close();


  // Find channel indices based on metadata and data
  vector<unsigned int> v_heaterIndex,
                       v_coolerIndex;
  vector<string> v_label;
  splitLine(line, v_label, ',');
  if (v_label.at(0) != "Item")
  {
    cout<<"ERROR: Starting of data labels with 'Item' after metadata was not found. Aborting"<<endl;
    return 0;
  }
  for (unsigned int i = 0; i < n_heaterThermistors; ++i)
  {
    for (unsigned int j = 1; j < v_label.size(); ++j)
    {
      if (v_heaterChannel.at(i) + "_through" == v_label.at(j).substr(2, 11))
      {
        cout<<"LOG: Found channel "<<v_heaterChannel.at(i)<<endl;
        v_heaterIndex.push_back(j);
      }
    }
  }
  for (unsigned int i = 0; i < n_coolerThermistors; ++i)
  {
    for (unsigned int j = 1; j < v_label.size(); ++j)
    {
      if (v_coolerChannel.at(i) + "_through" == v_label.at(j).substr(2, 11))
      {
        cout<<"LOG: Found channel "<<v_coolerChannel.at(i)<<endl;
        v_coolerIndex.push_back(j);
      }
    }
  }
  if (v_heaterIndex.size() == 0 || v_coolerIndex.size() == 0)
  {
    cout<<"ERROR: Number of heater channels found = "<<v_heaterIndex.size()<<", cooler channels found = "<<v_coolerIndex.size()<<". This is insufficient for data analysis. Aborting."<<endl;
    return 0;
  }

  // Extract data
  vector<vector<double> > v_v_heaterTemp(n_heaterThermistors),
                          v_v_coolerTemp(n_coolerThermistors),
                          v_v_heaterTime(n_heaterThermistors),
                          v_v_coolerTime(n_coolerThermistors);
  cout<<"LOG: Starting data extraction."<<endl;
  while (getline(ifs_rawData, line))
  {
    vector<string> v_line;
    splitLine(line, v_line, ',');

    for (unsigned int i = 0; i < n_heaterThermistors; ++i)
    {
      v_v_heaterTemp.at(i).push_back(stod(v_line.at(v_heaterIndex.at(i))));
      v_v_heaterTime.at(i).push_back(stod(v_line.at(v_heaterIndex.at(i)+1)));
    }
    for (unsigned int i = 0; i < n_coolerThermistors; ++i)
    {
      v_v_coolerTemp.at(i).push_back(stod(v_line.at(v_coolerIndex.at(i))));
      v_v_coolerTime.at(i).push_back(stod(v_line.at(v_coolerIndex.at(i)+1)));
    }
  }
  ifs_rawData.close();
  cout<<"LOG: Completed data extraction."<<endl;

  // Compute time between two data points of the same thermistor
  float dt = 0;
  for (unsigned int i = 0; i < n_heaterThermistors; ++i)
    dt += v_v_heaterTime.at(i).at(1) - v_v_heaterTime.at(i).at(0);
  for (unsigned int i = 0; i < n_coolerThermistors; ++i)
    dt += v_v_coolerTime.at(i).at(1) - v_v_coolerTime.at(i).at(0);
  dt /= float(n_heaterThermistors + n_coolerThermistors);
  cout<<"LOG: Time between two consecutive data points on the same thermistor = "<<dt<<" sec."<<endl;

  // Select region with lowest standard deviation greater than 1,800 seconds
  unsigned int i_min, i_max, i_size;
  unsigned int time_region_size=minIntegrationTime/dt;
  cout<<"LOG: Starting finding stable temperature region."<<endl;
  selectRegion(v_v_heaterTemp, v_v_coolerTemp, time_region_size, i_min, i_max, i_size);
  cout<<"LOG: Completed finding stable temperature region."<<endl;
  float time_min = v_v_heaterTime.at(0).at(i_min);
  float time_max = v_v_heaterTime.at(0).at(i_max-1);

  // Plot the whole graph of the raw data
  vector<int> v_heater_color = {kRed+4, kRed+3, kRed+2, kRed+1, kRed, kRed-4};
  vector<int> v_cooler_color = {kBlue+4, kBlue+3, kBlue+2, kBlue+1, kBlue, kBlue-4};
  vector<int> v_plotting_style = {20, 21, 22, 23, 33, 34};
  TMultiGraph *g_WholeGraph = new TMultiGraph();
  g_WholeGraph->SetTitle("; Time (s); Temperature (^{#circ}C)");
  for (unsigned int i = 0; i < n_heaterThermistors; ++i)
  {
    TGraph *aux = new TGraph(v_v_heaterTemp.at(i).size(), &v_v_heaterTime.at(i)[0], &v_v_heaterTemp.at(i)[0]);
    aux->SetLineColor(v_heater_color.at(i));
    g_WholeGraph->Add(aux);
  }
  for (unsigned int i = 0; i < n_coolerThermistors; ++i)
  {
    TGraph *aux = new TGraph(v_v_coolerTemp.at(i).size(), &v_v_coolerTime.at(i)[0], &v_v_coolerTemp.at(i)[0]);
    aux->SetLineColor(v_cooler_color.at(i));
    g_WholeGraph->Add(aux);
  }
  TCanvas *c_WholeGraph = new TCanvas("c_WholeGraph");
  g_WholeGraph->Draw("A");
  float y_max = g_WholeGraph->GetYaxis()->GetXmax();
  float y_min = g_WholeGraph->GetYaxis()->GetXmin();
  TArrow *arrow_min = new TArrow(time_min, y_max, time_min, y_min, 0.01, "|>");
  TArrow *arrow_max = new TArrow(time_max, y_max, time_max, y_min, 0.01, "|>");
  arrow_min->Draw();
  arrow_max->Draw();
  c_WholeGraph->SaveAs("c_WholeGraph.png");
  c_WholeGraph->SaveAs("c_WholeGraph.pdf");

  // Fill 1 min average of the selected integration region
  // Fill difference plots of the integration region
  double nStep = int(60./dt);
  vector<vector<double> > v_v_heaterRegionTemp(n_heaterThermistors),
                          v_v_coolerRegionTemp(n_coolerThermistors),
                          v_v_heaterRegionTime(n_heaterThermistors),
                          v_v_coolerRegionTime(n_coolerThermistors),
                          v_v_heaterRegionTemp_Err(n_heaterThermistors),
                          v_v_coolerRegionTemp_Err(n_heaterThermistors);
  vector<vector<double> > v_v_heaterRegionDiffTemp(n_heaterThermistors),
                          v_v_coolerRegionDiffTemp(n_coolerThermistors),
                          v_v_heaterRegionDiffTemp_Err(n_heaterThermistors),
                          v_v_coolerRegionDiffTemp_Err(n_coolerThermistors);
  for (unsigned int i = i_min; i < i_max - nStep; i += nStep)
  {
    for (unsigned int j = 0; j < n_heaterThermistors; ++j)
    {
      double temp0 = 0, temp1 = 0, temp2 = 0, temp3 = 0;
      for (unsigned int k = 0; k < nStep; ++k)
      {
        temp0 += v_v_heaterTemp.at(0).at(i+k);
        temp1 += v_v_heaterTemp.at(j).at(i+k);
        temp2 += pow(v_v_heaterTemp.at(j).at(i+k), 2);
        temp3 += pow(v_v_heaterTemp.at(j).at(i+k) - v_v_heaterTemp.at(0).at(i+k), 2);
      }
      temp0 /= nStep;
      temp1 /= nStep;
      double stdDev = pow(temp2 - nStep*temp1*temp1, 0.5);
      double diffStdDev = pow(temp3 - nStep*pow(temp1 - temp0, 2), 0.5);
      v_v_heaterRegionTemp.at(j).push_back(temp1);
      v_v_heaterRegionTemp_Err.at(j).push_back(stdDev);
      v_v_heaterRegionTime.at(j).push_back(v_v_heaterTime.at(j).at(i));
      v_v_heaterRegionDiffTemp.at(j).push_back(temp1 - temp0);
      v_v_heaterRegionDiffTemp_Err.at(j).push_back(diffStdDev);
    }
    for (unsigned int j = 0; j < n_coolerThermistors; ++j)
    {
      double temp0 = 0, temp1 = 0, temp2 = 0, temp3 = 0;
      for (unsigned int k = 0; k < nStep; ++k)
      {
        temp0 += v_v_coolerTemp.at(0).at(i+k);
        temp1 += v_v_coolerTemp.at(j).at(i+k);
        temp2 += pow(v_v_coolerTemp.at(j).at(i+k), 2);
        temp3 += pow(v_v_coolerTemp.at(j).at(i+k) - v_v_coolerTemp.at(0).at(i+k), 2);
      }
      temp0 /= nStep;
      temp1 /= nStep;
      double stdDev = pow(temp2 - nStep*temp1*temp1, 0.5);
      double diffStdDev = pow(temp3 - nStep*pow(temp1 - temp0, 2), 0.5);
      v_v_coolerRegionTemp.at(j).push_back(temp1);
      v_v_coolerRegionTemp_Err.at(j).push_back(stdDev);
      v_v_coolerRegionTime.at(j).push_back(v_v_coolerTime.at(j).at(i));
      v_v_coolerRegionDiffTemp.at(j).push_back(temp1 - temp0);
      v_v_coolerRegionDiffTemp_Err.at(j).push_back(diffStdDev);
    }
  }

  // Fit the basic Region plots
  vector<double> v_heaterRegionAvg(n_heaterThermistors),
                 v_heaterRegionErr(n_heaterThermistors),
                 v_coolerRegionAvg(n_coolerThermistors),
                 v_coolerRegionErr(n_coolerThermistors);
  TMultiGraph *g_HeaterGraph = new TMultiGraph();
  g_HeaterGraph->SetTitle("; Time (s); Temperature (^{#circ}C)");
  TLegend *leg_HeaterGraph = new TLegend(0.50, 0.70, 0.89, 0.89);
  leg_HeaterGraph->SetLineColor(0);
  leg_HeaterGraph->SetFillColor(0);
  for (unsigned int i = 0; i < n_heaterThermistors; ++i)
  {
    TGraphErrors *aux = new TGraphErrors(v_v_heaterRegionTemp.at(i).size(), &v_v_heaterRegionTime.at(i)[0], &v_v_heaterRegionTemp.at(i)[0], 0, &v_v_heaterRegionTemp_Err.at(i)[0]);
    aux->SetLineColor(v_heater_color.at(i));
    aux->SetMarkerColor(v_heater_color.at(i));
    aux->SetMarkerStyle(v_plotting_style.at(i));
    aux->SetMarkerSize(1);
    TF1 *f_aux = new TF1("f_aux", "[0]");
    f_aux->SetLineColor(v_heater_color.at(i));
    aux->Fit(f_aux, "Q");
    v_heaterRegionAvg.at(i)=f_aux->GetParameter(0);
    v_heaterRegionErr.at(i)=f_aux->GetParError(0);
    g_HeaterGraph->Add(aux);
    leg_HeaterGraph->AddEntry(aux, ("Hot Thermistor "+to_string(i+1)).c_str(), "P");
  }
  TMultiGraph *g_CoolerGraph = new TMultiGraph();
  g_CoolerGraph->SetTitle("; Time (s); Temperature (^{#circ}C)");
  TLegend *leg_CoolerGraph = new TLegend(0.50, 0.70, 0.89, 0.89);
  leg_CoolerGraph->SetLineColor(0);
  leg_CoolerGraph->SetFillColor(0);
  for (unsigned int i = 0; i < n_coolerThermistors; ++i)
  {
    TGraphErrors *aux = new TGraphErrors(v_v_coolerRegionTemp.at(i).size(), &v_v_coolerRegionTime.at(i)[0], &v_v_coolerRegionTemp.at(i)[0], 0, &v_v_coolerRegionTemp_Err.at(i)[0]);
    aux->SetLineColor(v_cooler_color.at(i));
    aux->SetMarkerColor(v_cooler_color.at(i));
    aux->SetMarkerStyle(v_plotting_style.at(i));
    aux->SetMarkerSize(1);
    TF1 *f_aux = new TF1("f_aux", "[0]");
    f_aux->SetLineColor(v_cooler_color.at(i));
    aux->Fit(f_aux, "Q");
    v_coolerRegionAvg.at(i)=f_aux->GetParameter(0);
    v_coolerRegionErr.at(i)=f_aux->GetParError(0);
    g_CoolerGraph->Add(aux);
    leg_CoolerGraph->AddEntry(aux, ("Cold Thermistor "+to_string(i+1)).c_str(), "P");
  }

  TCanvas *c_HeaterGraph = new TCanvas("c_HeaterGraph");
  g_HeaterGraph->SetMaximum(g_HeaterGraph->GetYaxis()->GetXmax()+0.8);
  g_HeaterGraph->SetMinimum(g_HeaterGraph->GetYaxis()->GetXmin()-0.15);
  g_HeaterGraph->Draw("AP");
  leg_HeaterGraph->Draw();
  c_HeaterGraph->SaveAs("c_HeaterGraph.png");
  c_HeaterGraph->SaveAs("c_HeaterGraph.pdf");

  TCanvas *c_CoolerGraph = new TCanvas("c_CoolerGraph");
  g_CoolerGraph->SetMaximum(g_CoolerGraph->GetYaxis()->GetXmax()+0.8);
  g_CoolerGraph->SetMinimum(g_CoolerGraph->GetYaxis()->GetXmin()-0.15);
  g_CoolerGraph->Draw("AP");
  leg_CoolerGraph->Draw();
  c_CoolerGraph->SaveAs("c_CoolerGraph.png");
  c_CoolerGraph->SaveAs("c_CoolerGraph.pdf");

  // Fit the Diff Region plots
  vector<double> v_heaterRegionDiffAvg(n_heaterThermistors),
                 v_heaterRegionDiffErr(n_heaterThermistors),
                 v_coolerRegionDiffAvg(n_coolerThermistors),
                 v_coolerRegionDiffErr(n_coolerThermistors);
  TMultiGraph *g_HeaterDiffGraph = new TMultiGraph();
  g_HeaterDiffGraph->SetTitle("; Time (s); Temperature (^{#circ}C)");
  for (unsigned int i = 0; i < n_heaterThermistors; ++i)
  {
    TGraphErrors *aux = new TGraphErrors(v_v_heaterRegionDiffTemp.at(i).size(), &v_v_heaterRegionTime.at(i)[0], &v_v_heaterRegionDiffTemp.at(i)[0], 0, &v_v_heaterRegionDiffTemp_Err.at(i)[0]);
    aux->SetLineColor(v_heater_color.at(i));
    aux->SetMarkerColor(v_heater_color.at(i));
    aux->SetMarkerStyle(v_plotting_style.at(i));
    aux->SetMarkerSize(1);
    TF1 *f_aux = new TF1("f_aux", "[0]");
    f_aux->SetLineColor(v_heater_color.at(i));
    aux->Fit(f_aux, "Q");
    v_heaterRegionDiffAvg.at(i)=f_aux->GetParameter(0);
    if (i != 0) v_heaterRegionDiffErr.at(i)=f_aux->GetParError(0);
    else v_heaterRegionDiffErr.at(i)=0.;
    g_HeaterDiffGraph->Add(aux);
  }
  TMultiGraph *g_CoolerDiffGraph = new TMultiGraph();
  g_CoolerDiffGraph->SetTitle("; Time (s); Temperature (^{#circ}C)");
  for (unsigned int i = 0; i < n_coolerThermistors; ++i)
  {
    TGraphErrors *aux = new TGraphErrors(v_v_coolerRegionDiffTemp.at(i).size(), &v_v_coolerRegionTime.at(i)[0], &v_v_coolerRegionDiffTemp.at(i)[0], 0, &v_v_coolerRegionDiffTemp_Err.at(i)[0]);
    aux->SetLineColor(v_cooler_color.at(i));
    aux->SetMarkerColor(v_cooler_color.at(i));
    aux->SetMarkerStyle(v_plotting_style.at(i));
    aux->SetMarkerSize(1);
    TF1 *f_aux = new TF1("f_aux", "[0]");
    f_aux->SetLineColor(v_cooler_color.at(i));
    aux->Fit(f_aux, "Q");
    v_coolerRegionDiffAvg.at(i)=f_aux->GetParameter(0);
    if (i != 0) v_coolerRegionDiffErr.at(i)=f_aux->GetParError(0);
    else v_coolerRegionDiffErr.at(i)=0.;
    g_CoolerDiffGraph->Add(aux);
  }

  TCanvas *c_HeaterDiffGraph = new TCanvas("c_HeaterDiffGraph");
  g_HeaterDiffGraph->SetMaximum(g_HeaterDiffGraph->GetYaxis()->GetXmax()+0.8);
  g_HeaterDiffGraph->SetMinimum(g_HeaterDiffGraph->GetYaxis()->GetXmin()-0.15);
  g_HeaterDiffGraph->Draw("AP");
  leg_HeaterGraph->Draw();
  c_HeaterDiffGraph->SaveAs("c_HeaterDiffGraph.png");
  c_HeaterDiffGraph->SaveAs("c_HeaterDiffGraph.pdf");

  TCanvas *c_CoolerDiffGraph = new TCanvas("c_CoolerDiffGraph");
  g_CoolerDiffGraph->SetMaximum(g_CoolerDiffGraph->GetYaxis()->GetXmax()+0.8);
  g_CoolerDiffGraph->SetMinimum(g_CoolerDiffGraph->GetYaxis()->GetXmin()-0.15);
  g_CoolerDiffGraph->Draw("AP");
  leg_CoolerGraph->Draw();
  c_CoolerDiffGraph->SaveAs("c_CoolerDiffGraph.png");
  c_CoolerDiffGraph->SaveAs("c_CoolerDiffGraph.pdf");

  cout<<"==========================================="<<endl;

  // Write datacard which should be readable by CSV and HTML readers
  ofstream ofs_datacard("Datacard.html");
  ofs_datacard<<"<pre>"<<endl;
  ofs_datacard<<"Through-plane thermal conductivity measurement"<<endl;
  ofs_datacard<<"Experimentalist, "<<meta_experimentalist<<endl;
  ofs_datacard<<"Date of measurement, "<<meta_date<<endl;
  ofs_datacard<<"Apparatus, "<<meta_apparatus<<endl;
  ofs_datacard<<"Material, "<<meta_material<<endl;
  ofs_datacard<<"Sample, "<<meta_sample<<endl;
  ofs_datacard<<"Thickness (um), "<<meta_thickness<<", "<<meta_thicknessErr<<endl;
  ofs_datacard<<"Heater Voltage, "<<meta_heaterVoltage<<endl;
  ofs_datacard<<"Cooler Voltage, "<<meta_coolerVoltage<<endl;
  ofs_datacard<<"Hot Fluxmeter Channels, "<<v_heaterChannel.at(0);
  for (unsigned int i = 1; i < n_heaterThermistors; ++i)
    ofs_datacard<<", "<<v_heaterChannel.at(i);
  ofs_datacard<<endl;
  ofs_datacard<<"Cold Fluxmeter Channels, "<<v_coolerChannel.at(0);
  for (unsigned int i = 1; i < n_coolerThermistors; ++i)
    ofs_datacard<<", "<<v_coolerChannel.at(i);
  ofs_datacard<<endl;
  ofs_datacard<<"Hot Fluxmeter Temperatures (C), "<<v_heaterRegionAvg.at(0);
  for (unsigned int i = 1; i < n_heaterThermistors; ++i)
    ofs_datacard<<", "<<v_heaterRegionAvg.at(i);
  ofs_datacard<<endl;
  ofs_datacard<<"Cold Fluxmeter Temperatures (C), "<<v_coolerRegionAvg.at(0);
  for (unsigned int i = 1; i < n_coolerThermistors; ++i)
    ofs_datacard<<", "<<v_coolerRegionAvg.at(i);
  ofs_datacard<<endl;
  ofs_datacard<<"Hot Fluxmeter Temperature Uncertainties (C), "<<v_heaterRegionErr.at(0);
  for (unsigned int i = 1; i < n_heaterThermistors; ++i)
    ofs_datacard<<", "<<v_heaterRegionErr.at(i);
  ofs_datacard<<endl;
  ofs_datacard<<"Cold Fluxmeter Temperature Uncertainties (C), "<<v_coolerRegionErr.at(0);
  for (unsigned int i = 1; i < n_coolerThermistors; ++i)
    ofs_datacard<<", "<<v_coolerRegionErr.at(i);
  ofs_datacard<<endl;
  ofs_datacard<<"Hot Fluxmeter Temperatures Relative to Hottest Thermistor (C), "<<v_heaterRegionDiffAvg.at(0);
  for (unsigned int i = 1; i < n_heaterThermistors; ++i)
    ofs_datacard<<", "<<v_heaterRegionDiffAvg.at(i);
  ofs_datacard<<endl;
  ofs_datacard<<"Cold Fluxmeter Temperatures Relative to Hottest Thermistor (C), "<<v_coolerRegionDiffAvg.at(0);
  for (unsigned int i = 1; i < n_coolerThermistors; ++i)
    ofs_datacard<<", "<<v_coolerRegionDiffAvg.at(i);
  ofs_datacard<<endl;
  ofs_datacard<<"Hot Fluxmeter Temperature Relative to Hottest Thermistor Uncertainties (C), "<<v_heaterRegionDiffErr.at(0);
  for (unsigned int i = 1; i < n_heaterThermistors; ++i)
    ofs_datacard<<", "<<v_heaterRegionDiffErr.at(i);
  ofs_datacard<<endl;
  ofs_datacard<<"Cold Fluxmeter Temperature Relative to Hottest Thermistor Uncertainties (C), "<<v_coolerRegionDiffErr.at(0);
  for (unsigned int i = 1; i < n_coolerThermistors; ++i)
    ofs_datacard<<", "<<v_coolerRegionDiffErr.at(i);
  ofs_datacard<<endl;
  ofs_datacard<<"</pre>"<<endl;
  ofs_datacard<<"<img src='c_WholeGraph.png'/><br/>"<<endl;
  ofs_datacard<<"<table border='1'>"<<endl;
  ofs_datacard<<" <tr>"<<endl;
  ofs_datacard<<"  <td>"<<endl;
  ofs_datacard<<"   Heater Thermistors <br/>"<<endl;
  ofs_datacard<<"   <img src='c_HeaterGraph.png'/>"<<endl;
  ofs_datacard<<"  </td>"<<endl;
  ofs_datacard<<"  <td>"<<endl;
  ofs_datacard<<"   Cooler Thermistors <br/>"<<endl;
  ofs_datacard<<"   <img src='c_CoolerGraph.png'/>"<<endl;
  ofs_datacard<<"  </td>"<<endl;
  ofs_datacard<<" </tr>"<<endl;
  ofs_datacard<<" <tr>"<<endl;
  ofs_datacard<<"  <td>"<<endl;
  ofs_datacard<<"   Heater Thermistors relative to hottest thermistor <br/>"<<endl;
  ofs_datacard<<"   <img src='c_HeaterDiffGraph.png'/>"<<endl;
  ofs_datacard<<"  </td>"<<endl;
  ofs_datacard<<"  <td>"<<endl;
  ofs_datacard<<"   Cooler Thermistors relative to the hottest thermistor <br/>"<<endl;
  ofs_datacard<<"   <img src='c_CoolerDiffGraph.png'/><br/>"<<endl;
  ofs_datacard<<"  </td>"<<endl;
  ofs_datacard<<" </tr>"<<endl;
  ofs_datacard<<"</table>"<<endl;
  ofs_datacard.close();
  cout<<"Datacard in Datacard.html"<<endl;
  system("open Datacard.html");
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

// Selects region with lowest stddev
void selectRegion(vector<vector<double>> &heater_sample,vector<vector<double>> &peltier_sample,unsigned int time_size,unsigned int &min,unsigned int &max,unsigned int &size){

    int chunk = heater_sample[0].size();

    double min_avg_stddev=99.0;

    double sum_avg_stddev=0;

    for(unsigned int ii=chunk;ii>=time_size;ii--){
        int end = 1+ chunk-ii;
        for(int kk=0;kk<end;kk++){
            vector<double> heater_stddev(heater_sample.size());
            vector<double> peltier_stddev(peltier_sample.size());
            double heater_avg_stddev=0;
            double peltier_avg_stddev=0;
            for(unsigned int ll=0;ll<heater_sample.size();ll+=1){
                heater_stddev[ll]=0;
            }
            for(unsigned int ll=0;ll<peltier_sample.size();ll+=1){
                peltier_stddev[ll]=0;
            }
            for(unsigned int ll=0;ll<heater_sample.size();ll+=1){
                double heater_avg=0;
                for(unsigned int jj=kk;jj<kk+ii;jj++){
                    heater_avg+=heater_sample[ll][jj];
                }
                heater_avg=heater_avg/ii;
                for(unsigned int jj=kk;jj<kk+ii;jj++){
                    heater_stddev[ll]+=(heater_sample[ll][jj]-heater_avg)*(heater_sample[ll][jj]-heater_avg);
                }
                heater_stddev[ll]=heater_stddev[ll]/(ii-1);
                heater_stddev[ll]=sqrt(heater_stddev[ll]);
            }

            for(unsigned int ll=0;ll<peltier_sample.size();ll+=1){
                double peltier_avg=0;
                for(unsigned int jj=kk;jj<kk+ii;jj++){
                    peltier_avg+=peltier_sample[ll][jj];
                }
                peltier_avg=peltier_avg/ii;
                for(unsigned int jj=kk;jj<kk+ii;jj++){
                    peltier_stddev[ll]+=(peltier_sample[ll][jj]-peltier_avg)*(peltier_sample[ll][jj]-peltier_avg);
                }
                peltier_stddev[ll]=peltier_stddev[ll]/(ii-1);
                peltier_stddev[ll]=sqrt(peltier_stddev[ll]);
            }
            for(unsigned int ll=0;ll<heater_sample.size();ll+=1){
                heater_avg_stddev+=heater_stddev[ll];
            }

            for(unsigned int ll=0;ll<peltier_sample.size();ll+=1){
                peltier_avg_stddev+=peltier_stddev[ll];
            }
            heater_avg_stddev=heater_avg_stddev/heater_stddev.size();
            peltier_avg_stddev=peltier_avg_stddev/peltier_stddev.size();
            sum_avg_stddev=heater_avg_stddev*heater_avg_stddev+peltier_avg_stddev*peltier_avg_stddev;
            if(sum_avg_stddev<=min_avg_stddev) min_avg_stddev=sum_avg_stddev,min=kk,max=min+ii,size=ii;
        }
    }
}
