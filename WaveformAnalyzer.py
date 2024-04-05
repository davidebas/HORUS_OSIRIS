# -*- coding: utf-8 -*-
# OSIRIS WAVEFORM BASIC ANALYZER
# Author: Davide Basilico davide.basilico@mi.infn.it  2024 April 05
import sys
import uproot
import re
import matplotlib.pyplot as plt
import numpy as np
import time
import awkward as ak
from datetime import datetime
from Waveform import Waveform
from tqdm import tqdm
from utils import extract_date_and_formatted_date, read_all_PMTs_coordinates, load_PMT_coordinates, find_PMT_Coordinates, print_progress

def main():

	if len(sys.argv) != 3:
		print("Usage: python WF_Analyzer.py [Input_RootFile] [Output_TxtFile]")
		sys.exit(1)

	file = uproot.open(sys.argv[1])
	tree = file['EventTree']
	Entries = int(tree.num_entries) - 1 

	# Read all the PMTs coordinates (x,y,z)
	all_PMTs_x, all_PMTs_y, all_PMTs_z, all_PMTs_theta = read_all_PMTs_coordinates("OSIRIS_cable_map.conf")

	with open(sys.argv[2], 'w'):
		pass

	extracted_string, extracted_date = extract_date_and_formatted_date(sys.argv[1])

	print("### Welcome to WaveformAnalyzer ###")	
	print("Run: " , extracted_string, "\nDate:", extracted_date)	

	with open(sys.argv[2], 'a') as file:
		for j in tqdm(range(0, Entries), total=Entries, desc="Processing", leave=True,   ):
			events = tree.arrays(["eventId","trgNsec","IDdata.samples","IDdata.GCUID","IDdata.channelID"], entry_start=j, entry_stop=j+1)
			trgNsec = ak.to_numpy(events["trgNsec"])
			IDdata_GCUID = ak.to_numpy(events["IDdata.GCUID"])
			IDdata_channelID = ak.to_numpy(events["IDdata.channelID"])			
			IDdata_samples_vector = ak.to_numpy(events["IDdata.samples"])

			for i in range(0, IDdata_channelID.shape[1]-1):

				if IDdata_channelID[0,i] % 2 == 0:  # Only high-gains
					continue

				waveform = Waveform(IDdata_samples_vector[0, i, :], threshold_method='std_dev', baseline_entries=50)

				fired_PMTs, integrated_charge, idx = waveform.analyze_waveform()


				if fired_PMTs > 0:
					coordinate = find_PMT_Coordinates(IDdata_GCUID[0,i], IDdata_channelID[0,i])
					time_3 = time.time()
					with open(sys.argv[2], 'a') as file:
						file.write(f"{extracted_string}\t{extracted_date}\t{j}\t{integrated_charge}\t{i}\t{idx}\t{trgNsec[0]}\t{IDdata_channelID[0, i]}	\t{IDdata_GCUID[0, i]}\t{coordinate[0]}\t{coordinate[1]}\t{coordinate[2]}\t{int(IDdata_channelID.shape[1]/2)}\n")



if __name__ == "__main__":
	time_start = time.time()
	main()
	time_end = time.time()
	print(f"Completed in: {time_end-time_start:.2f} s")
