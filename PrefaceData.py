# Preface raw data with experimentalist's observations
# Author: Souvik Das, souvik@purdue.edu
# December 2022
# Written for Python 3

import csv
import os
import tkinter
from tkinter.filedialog import askopenfile
from shutil import move

filename = tkinter.filedialog.askopenfilename(initialdir="../ThermalConductivityRawData")
csv_inFile = open(filename, "r")
reader = csv.reader(csv_inFile)
row = next(reader)
if (row[0] == "Experimentalist"):
  print("ERROR: Raw file has already been annotated. Aborting.", flush=True)
  quit()
csv_inFile.seek(0)

print("Please enter the following information.")
name = input("Experimentalist's name (First, Last): ")
date = input("Date of the experiment (mm/dd/yyyy): ")
apparatus = input("Apparatus (1 or 2): ")
material = input("Material codename (https://docs.google.com/spreadsheets/d/1PuO5Ts9dGGj1flBs73u0V22GOnD-rF94tCsjXpQOD58/edit): ")
sample = input("Sample number (1 - 10 typically): ")
thickness = input("Sample thickness (um): ")
thickness_err = input("Sample thickness uncertainty (um): ")
heater_voltage = input("Heater Voltage (V): ")
cooler_voltage = input("Cooler Voltage (V): ")

csv_outFile = open(filename+"_temp", "w", newline='')
writer = csv.writer(csv_outFile)
writer.writerow(("Experimentalist", name))
writer.writerow(("Date", date))
writer.writerow(("Apparatus", apparatus))
writer.writerow(("Material", material))
writer.writerow(("Sample", sample))
writer.writerow(("Thickness (um)", thickness, thickness_err))
writer.writerow(("Heater Voltage (V)", heater_voltage))
writer.writerow(("Cooler Voltage (V)", cooler_voltage))
writer.writerow('')
for row in reader:
  # print(row)
  writer.writerow(row)

csv_inFile.close()
csv_outFile.close()
move(filename+"_temp", filename)
