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

class Waveform:
	def __init__(self, samples, threshold_method, baseline_entries):
		self.samples = samples
		self.min_value = np.min(samples)
		self.baseline = np.mean(samples[:baseline_entries])  # Average of earliest baseline_entries of WF
		self.std_dev_baseline = np.std(samples[:baseline_entries])
		self.threshold_method = threshold_method

	def check_threshold(self, method):
		if method == 'baseline':
			threshold_absolute = 20
			return self.baseline - self.min_value > threshold_absolute
		elif method == 'std_dev':
			threshold_std_dev = 5
			return (self.baseline - self.min_value) > threshold_std_dev * self.std_dev_baseline
		else:
			raise ValueError("Threshold finder method not valid. Please use 'baseline' or 'std_dev'.")

	def analyze_waveform(self):
		fired_PMTs = 0
		integrated_charge = 0
		threshold_crossing_index = None
		
		#threshold_method = 'baseline'
		
		if self.check_threshold(self.threshold_method):
			diff_baseline_max = self.baseline - self.min_value

			# Find the index WF for the rise time
			for idx, val in enumerate(self.samples):
				if abs(self.baseline - val) > abs(diff_baseline_max) / 5:  # Rise time threshold
					fired_PMTs += 1
					integrated_charge += np.sum(self.samples[idx:] - self.baseline)
					threshold_crossing_index = idx  # Memorize the index associated to rise time
					break 

		return fired_PMTs, integrated_charge, threshold_crossing_index

def extract_date_and_formatted_date(input_string):
	match = re.search(r'(\d{8})_(\d{6})', input_string) # Isolate the date 
	if match:
		extracted_string = match.group(0)  # Full string
		extracted_date_string = match.group(1)  # Only date
	try:
		extracted_date = datetime.strptime(extracted_date_string, "%Y%m%d").date()
		formatted_date = extracted_date.strftime("%Y-%m-%d")
		return extracted_string, formatted_date
	except ValueError:
		return None, None

def read_all_PMTs_coordinates(nome_file):
	x = np.array([])
	y = np.array([])
	z = np.array([])
	theta = np.array([])
	with open(nome_file, 'r') as file:
		for line in file:
			columns = line.strip().split('\t')
			try:
				x = np.append(x, float(columns[-3]))
				y = np.append(y, float(columns[-2]))
				z = np.append(z, float(columns[-1]))
				theta = np.append(theta, np.sign(float(columns[-2])) * np.arccos(float(columns[-3]) / np.sqrt(float(columns[-3])**2 + float(columns[-2])**2)))
			except ValueError:
				continue
	return x, y, z, theta

def load_PMT_coordinates(file_path):
    PMT_coordinates = {}
    try:
        with open(file_path, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                if len(columns) == 14:
                    GCUID = int(columns[1])
                    GCUChannel_true = int(columns[3]) // 2
                    if GCUID not in PMT_coordinates:
                        PMT_coordinates[GCUID] = {}
                    PMT_coordinates[GCUID][GCUChannel_true] = (columns[11], columns[12], columns[13])
    except FileNotFoundError:
        print("PMT coordinates file not found.")
    return PMT_coordinates

PMT_coordinates_data = load_PMT_coordinates("OSIRIS_cable_map.conf")

def find_PMT_Coordinates(GCUID_user, GCUChannel_user):
    GCUID_user = int(GCUID_user)
    GCUChannel_user_true = int(GCUChannel_user) // 2
    if GCUID_user in PMT_coordinates_data and GCUChannel_user_true in PMT_coordinates_data[GCUID_user]:
        return PMT_coordinates_data[GCUID_user][GCUChannel_user_true]
    else:
        return None

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
	
	fired_PMTs_Cut = 0

	print(extracted_string, extracted_date)	

	with open(sys.argv[2], 'a') as file:
		for j in range(0, 1000):
			events = tree.arrays(["eventId","trgNsec","IDdata.samples","IDdata.GCUID","IDdata.channelID"], entry_start=j, entry_stop=j+1)
			trgNsec = ak.to_numpy(events["trgNsec"])
			IDdata_GCUID = ak.to_numpy(events["IDdata.GCUID"])
			IDdata_channelID = ak.to_numpy(events["IDdata.channelID"])			
			IDdata_samples_vector = ak.to_numpy(events["IDdata.samples"])

			print("Analyzing event " + str(j))

			for i in range(0, IDdata_channelID.shape[1]-1):

				if IDdata_channelID[0,i] % 2 == 0:  # Only high-gains
					continue

				waveform = Waveform(IDdata_samples_vector[0, i, :], threshold_method='std_dev', baseline_entries=20)

				fired_PMTs, integrated_charge, idx = waveform.analyze_waveform()
				time_2 = time.time()

				if fired_PMTs > fired_PMTs_Cut:
					coordinate = find_PMT_Coordinates(IDdata_GCUID[0,i], IDdata_channelID[0,i])
					time_3 = time.time()
					with open(sys.argv[2], 'a') as file:
						file.write(f"{extracted_string}\t{extracted_date}\t{j}\t{integrated_charge}\t{i}\t{idx}\t{trgNsec[0]}\t{IDdata_channelID[0, i]}	\t{IDdata_GCUID[0, i]}\t{coordinate[0]}\t{coordinate[1]}\t{coordinate[2]}\t{int(IDdata_channelID.shape[1]/2)}\n")



if __name__ == "__main__":
	time_start = time.time()
	main()
	time_end = time.time()
	print("Run time: " + str(time_end-time_start))
