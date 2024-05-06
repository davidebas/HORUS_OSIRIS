# utils.py
import re
import numpy as np
from datetime import datetime
from numba import njit,jit

def print_progress(j, total_events, next_percentage):
    percentage_done = (j + 1) / total_events * 100
    if percentage_done >= next_percentage:
        print(f"Analyzing event {j} - {percentage_done:.0f}% done")
        next_percentage += 5
        return next_percentage
    return next_percentage
    
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
	gain_corr = np.array([])
	with open(nome_file, 'r') as file:
		for line in file:
			columns = line.strip().split('\t')
			try:
				x = np.append(x, float(columns[11]))
				y = np.append(y, float(columns[12]))
				z = np.append(z, float(columns[13]))
				gain_corr = np.append(gain_corr,float(columns[14]))
				#theta = np.append(theta, np.sign(float(columns[-2])) * np.arccos(float(columns[-3]) / np.sqrt(float(columns[-3])**2 + float(columns[-2])**2)))
			except ValueError:
				continue
	return x, y, z, gain_corr

def load_PMT_coordinates(file_path):
    PMT_coordinates = {}
    gain_corr = {}
    try:
        with open(file_path, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                if len(columns) > 10:
                    GCUID = int(columns[1])
                    GCUChannel = int(columns[3])
                    if GCUID not in PMT_coordinates:
                        PMT_coordinates[GCUID] = {}
                        gain_corr[GCUID] = {}
                    PMT_coordinates[GCUID][GCUChannel] = (columns[11], columns[12], columns[13])
                    gain_corr[GCUID][GCUChannel] = float(columns[14])
    except FileNotFoundError:
        print("PMT coordinates file not found.")
    return PMT_coordinates, gain_corr

PMT_coordinates_data, gain_corr_data = load_PMT_coordinates("OSIRIS_cable_map_N.conf")

def find_PMT_Coordinates(GCUID_data, GCUChannel_data):
    GCUID_data = int(GCUID_data)
    GCUChannel_data_halved = int(GCUChannel_data) // 2
    #print(PMT_coordinates_data[GCUID_data])
    #print("############")
    #print(GCUID_data,GCUChannel_data_halved)
    #print(PMT_coordinates_data[GCUID_data])
    if GCUID_data in PMT_coordinates_data and GCUChannel_data_halved in PMT_coordinates_data[GCUID_data]:
        #print("OK", GCUID_data, GCUChannel_data, GCUChannel_data_halved, PMT_coordinates_data[GCUID_data][GCUChannel_data_halved], gain_corr_data[GCUID_data][GCUChannel_data_halved])
        #print("############")
        return PMT_coordinates_data[GCUID_data][GCUChannel_data_halved], gain_corr_data[GCUID_data][GCUChannel_data_halved]
    else:
        print("OOOH")
        return None
        

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

	return synthetic_samples

