# -*- coding: utf-8 -*-
# OSIRIS WAVEFORM BASIC ANALYZER  v.0.1.0
# Author: Davide Basilico davide.basilico@mi.infn.it  2024 May 06

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
from utils import (
    extract_date_and_formatted_date,
    read_all_PMTs_coordinates,
    load_PMT_coordinates,
    find_PMT_Coordinates,
    print_progress,
    generate_synthetic_waveform
)
from numba import njit,jit

OD_ids = [(17,1), (18,1), (19,5), (20,1), (37, 3), (38,1), (38,5), (39,1), (40,1), (42, 1), (42, 5), (43, 3)]

def main():

	synth_mode = 'false'
	
	if len(sys.argv) < 3:
		print("Usage: python WF_Analyzer.py [Input_RootFile] [Output_TxtFile]")
		sys.exit(1)
		
	header_synth = "name	date	index	charge	WF_RiseTime	\n"
	header = "name\tdate\tindex\tcharge\tGCU\tWF_RiseTime\ttrgTime\tID_channel\tGCUID\tx_PMT\ty_PMT\tz_PMT\tLivePMTs\tgain\tOD\tShape_Ch\n"

	OD_fired = 0

	# Read all the PMTs coordinates (x,y,z)
	all_PMTs_x, all_PMTs_y, all_PMTs_z, all_PMTs_gain_corr = read_all_PMTs_coordinates("OSIRIS_cable_map_N.conf")

	with open(sys.argv[2], 'w'):
		pass

	print("### Welcome to WaveformAnalyzer ###")		

	if(synth_mode == 'false'):
	
		block_size = 100
		file = uproot.open(sys.argv[1])
		tree = file['EventTree']
		Entries = int(tree.num_entries) - 1 
	
		extracted_string, extracted_date = extract_date_and_formatted_date(sys.argv[1])
		print("Run: " , extracted_string, "\nDate:", extracted_date)		
		with open(sys.argv[2], 'a') as file:
			file.write(header)
	
			for start in tqdm(range(0, Entries, block_size), desc="Processing", leave=True):
				stop = min(start + block_size, Entries)
				events = tree.arrays(["eventId", "trgSec", "trgNsec", "IDdata.samples", "IDdata.GCUID", "IDdata.channelID"], entry_start=start, entry_stop=stop)
				trgTime = []

				for j in range(stop - start):
					#print(events["trgSec"][j],events["trgNsec"][j], events["IDdata.channelID"][j])
					#ak.to_numpy(events["trgSec"][j]) + 1e-9*ak.to_numpy(events["trgNsec"][j])
					trgTime.append(events["trgSec"][j] + 1e-9*events["trgNsec"][j])
					IDdata_GCUID = ak.to_numpy(events["IDdata.GCUID"][j])
					IDdata_channelID = ak.to_numpy(events["IDdata.channelID"][j])
					IDdata_samples_vector = ak.to_numpy(events["IDdata.samples"][j])

					for i in range(IDdata_channelID.shape[0]):
						
						if IDdata_channelID[i] % 2 == 0:  # Only high-gains
							continue
										
						waveform = Waveform(IDdata_samples_vector[i, :], threshold_method='std_dev', baseline_entries=50)
						fired_PMTs, integrated_charge, rise_time = waveform.analyze_waveform()

						if fired_PMTs > 0:
							coordinate, gain_corr = find_PMT_Coordinates(IDdata_GCUID[i], IDdata_channelID[i])
							integrated_charge = integrated_charge / float(gain_corr) * 100.

							if (IDdata_GCUID[i], IDdata_channelID[i]) in OD_ids:
								OD_fired = 1
								#print("GCU: ", IDdata_GCUID[i], "; channelID: ", IDdata_channelID[i])
							else:
								OD_fired = 0

							file.write(f"{extracted_string}\t{extracted_date}\t{start+j}\t{integrated_charge}\t{i}\t{rise_time}\t{trgTime[j]}\t{IDdata_channelID[i]}\t{IDdata_GCUID[i]}\t{coordinate[0]}\t{coordinate[1]}\t{coordinate[2]}\t{int(IDdata_channelID.shape[0]/2)}\t{gain_corr}\t{OD_fired}\t{IDdata_samples_vector.shape[0]}\n")
						

	if(synth_mode == 'true'):
	
		print("Synth mode on: generating synthetic waveforms")	
		
		extracted_string, extracted_date = '2024_X', '2024_MM_DD'
		
		
		with open(sys.argv[2], 'a') as file:
			file.write(header_synth)	
			Entries = int(sys.argv[3])	
			for start in tqdm(range(0, Entries), desc="Processing", leave=True):
				
				if(start % 2 ==0):
					synth_params = {'flat_height': 11000, 'gaussian_amplitude': -200,'gaussian_center': 250,'gaussian_width': 10}	
				else:
					synth_params = {'flat_height': 11000, 'gaussian_amplitude': 0,'gaussian_center': 250,'gaussian_width': 10}
					
				samples = generate_synthetic_waveform(synth_params)					
				waveform = Waveform(samples, threshold_method='std_dev', baseline_entries=50)
				fired_PMTs, integrated_charge, rise_time = waveform.analyze_waveform()
				#print(start,fired_PMTs, integrated_charge, rise_time)
				
				file.write(f"{extracted_string}\t{extracted_date}\t{start}\t{integrated_charge}\t{rise_time}\n")


if __name__ == "__main__":
	time_start = time.time()
	main()
	time_end = time.time()
	print(f"Completed in: {time_end-time_start:.2f} s")
