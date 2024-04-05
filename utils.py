# utils.py
import re
import numpy as np
from datetime import datetime

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
	#theta = np.array([])
	with open(nome_file, 'r') as file:
		for line in file:
			columns = line.strip().split('\t')
			try:
				x = np.append(x, float(columns[-3]))
				y = np.append(y, float(columns[-2]))
				z = np.append(z, float(columns[-1]))
				#theta = np.append(theta, np.sign(float(columns[-2])) * np.arccos(float(columns[-3]) / np.sqrt(float(columns[-3])**2 + float(columns[-2])**2)))
			except ValueError:
				continue
	return x, y, z, 0

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

