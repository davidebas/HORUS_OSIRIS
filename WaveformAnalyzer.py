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

synthetic_params = {
    'flat_height': 1.0,            # Altezza della componente piatta
    'gaussian_amplitude': 5.0,     # Amplitude della gaussiana
    'gaussian_center': 50,         # Centro della gaussiana
    'gaussian_width': 10           # Larghezza della gaussiana
}

def main():

	if len(sys.argv) != 4:
		print("Usage: python WF_Analyzer.py [Input_RootFile] [Output_TxtFile]")
		sys.exit(1)

	file = uproot.open(sys.argv[1])
	tree = file['EventTree']
	Entries = int(tree.num_entries) - 1 
	
	Entries = 2000

	# Read all the PMTs coordinates (x,y,z)
	all_PMTs_x, all_PMTs_y, all_PMTs_z, all_PMTs_gain_corr = read_all_PMTs_coordinates("OSIRIS_cable_map_N.conf")

	with open(sys.argv[2], 'w'):
		pass

	extracted_string, extracted_date = extract_date_and_formatted_date(sys.argv[1])

	print("### Welcome to WaveformAnalyzer ###")	
	print("Run: " , extracted_string, "\nDate:", extracted_date)
	
	synth_mode = 0
	synth_params = {}	
	if synth_mode == 1:
		synth_params = {
		'flat_height': 11000,            
		'gaussian_amplitude': -100,     
		'gaussian_center': 250,         
		'gaussian_width': 10           
		}


	with open(sys.argv[2], 'a') as file:
		for j in tqdm(range(0, Entries), total=Entries, desc="Processing", leave=True,   ):
			
			if synth_mode == 1:
				samples = generate_synthetic_waveform(synth_params)
				
				waveform = Waveform(samples, threshold_method='std_dev', baseline_entries=50)
				
				fired_PMTs, integrated_charge, idx = waveform.analyze_waveform()
				
				print(j,fired_PMTs, integrated_charge, idx)
			
			if synth_mode == 0:
				events = tree.arrays(["eventId","trgNsec","IDdata.samples","IDdata.GCUID","IDdata.channelID"], entry_start=j, entry_stop=j+1)
				trgNsec = ak.to_numpy(events["trgNsec"])
				IDdata_GCUID = ak.to_numpy(events["IDdata.GCUID"])
				IDdata_channelID = ak.to_numpy(events["IDdata.channelID"])			
				IDdata_samples_vector = ak.to_numpy(events["IDdata.samples"])

				for i in range(0, IDdata_channelID.shape[1]-1):
				

					if IDdata_channelID[0,i] % 2 == 0:  # Only high-gains
						continue

					waveform = Waveform(IDdata_samples_vector[0, i, :], threshold_method='std_dev', baseline_entries=50)
					
					#waveform = Waveform(samples, synthetic_waveform_params=synthetic_params)
	
					fired_PMTs, integrated_charge, idx = waveform.analyze_waveform()


					if fired_PMTs > 0:
						coordinate, gain_corr = find_PMT_Coordinates(IDdata_GCUID[0,i], IDdata_channelID[0,i])
						integrated_charge = integrated_charge / float(gain_corr) * 100.
						#print(IDdata_GCUID[0,i], IDdata_channelID[0,i], coordinate, gain_corr)
						#print(gain_corr)
						with open(sys.argv[2], 'a') as file:
							file.write(f"{extracted_string}\t{extracted_date}\t{j}\t{integrated_charge}\t{i}\t{idx}\t{trgNsec[0]}\t{IDdata_channelID[0, i]}	\t{IDdata_GCUID[0, i]}\t{coordinate[0]}\t{coordinate[1]}\t{coordinate[2]}\t{int(IDdata_channelID.shape[1]/2)}\t{gain_corr}\n")


def generate_synthetic_waveform(params):
    flat_height = params['flat_height']
    gaussian_amplitude = params['gaussian_amplitude']
    gaussian_center = params['gaussian_center']
    gaussian_width = params['gaussian_width']

    # Generazione della waveform sintetica
    x = np.arange(600)
    # Creazione del rumore di fondo con fluttuazioni casuali
    background_noise = np.random.normal(0, 5, size=600)  # Deviazione standard di 0.1 per esempio
    flat_component = np.ones_like(x) * flat_height + background_noise
    
    # Generazione del contributo gaussiano con fluttuazioni casuali
    gaussian_noise = np.random.normal(0, 0.05, size=600)  # Deviazione standard di 0.05 per esempio
    gaussian_component = np.exp(-0.5 * ((x - gaussian_center) / gaussian_width)**2) * (gaussian_amplitude + gaussian_noise)
    
    synthetic_samples = flat_component + gaussian_component
    
    #print(synthetic_samples)

    return synthetic_samples

if __name__ == "__main__":
	time_start = time.time()
	main()
	time_end = time.time()
	print(f"Completed in: {time_end-time_start:.2f} s")
