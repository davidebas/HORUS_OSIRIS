# -*- coding: utf-8 -*-
# OSIRIS WAVEFORM BASIC ANALYZER
# Author: Davide Basilico davide.basilico@mi.infn.it  2024 March 27

import sys
import uproot
import re
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
plt.rcParams.update({'font.size': 10})

def extract_date_and_formatted_date(input_string):
    match = re.search(r'(\d{8})_(\d{6})', input_string)
    if match:
        extracted_string = match.group(0)  # Entire string
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
                theta = np.append(theta,np.sign(float(columns[-2]))*np.arccos(float(columns[-3]) / np.sqrt(float(columns[-3])**2+float(columns[-2])**2)))
                
            except ValueError:
                continue
                
    return x, y, z, theta

def find_PMT_Coordinates(file_path, GCUID_user, GCUChannel_user):
# il GCUID del rootfile va associato al BECPort_ID, cioe colonna 1, del file excel
    try:
        with open(file_path, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')                
                GCUChannel_user_true = 777
                
                if(GCUChannel_user==0 or GCUChannel_user==1):
                	GCUChannel_user_true = 0
                if(GCUChannel_user==2 or GCUChannel_user==3):
                	GCUChannel_user_true = 1
                if(GCUChannel_user==4 or GCUChannel_user==5):
                	GCUChannel_user_true = 2              	
              
                if (len(columns) == 14) and (int(columns[1]) == int(GCUID_user)) and (int(columns[3]) == int(GCUChannel_user_true)):
                    return (columns[11], columns[12], columns[13])
    
    except FileNotFoundError:
        print("PMT coordinates file not found.")
    
    return None

file = uproot.open(sys.argv[1])
tree = file['EventTree']
Entries = int(tree.num_entries) -1 

# Read all the PMTs coordinates (x,y,z)
all_PMTs_x, all_PMTs_y, all_PMTs_z, all_PMTs_theta = read_all_PMTs_coordinates("OSIRIS_cable_map.conf")

with open(sys.argv[2], 'w'):
    pass

#fig_WF, ax_WF = plt.subplots(1,1, figsize=(8, 5))

extracted_string, extracted_date = extract_date_and_formatted_date(sys.argv[1])

print(extracted_string, extracted_date)

#with open(sys.argv[2], 'a') as file:
#	file.write("Name\tDate\tEvent\tCharge\t{i}\t{idx}\t{trgNsec[0]}\t{IDdata_channelID[0, i]}\t{IDdata_GCUID[0, i]}\t{coordinate[0]}\t{coordinate[1]}\t{coordinate[2]}\t{int(IDdata_channelID.shape[1]/2)}\n")

for j in range(0,Entries):

	events = tree.arrays(["eventId","trgNsec","IDdata.samples","IDdata.GCUID","IDdata.channelID"], library="pd",entry_start=j, entry_stop=j+1)
	trgNsec = np.array(events["trgNsec"])
	IDdata_samples_vector = np.array(events["IDdata.samples"])
	IDdata_GCUID = np.array(events["IDdata.GCUID"])
	IDdata_channelID = np.array(events["IDdata.channelID"])
	eventId = np.array(events["eventId"])
	
	print("Analyzing event " + str(j) + "/" + str(Entries))
	
	coordinate_x = np.array([])
	coordinate_y = np.array([])
	coordinate_z = np.array([])
	somme_differenze = np.array([])
	
	fired_PMTs = 0
	
	#[1,156,600]
	
	for i in range(0,IDdata_channelID.shape[1]-1):
		
#	for i in Entries(1, 79):
		#plt.figure(figsize=(10, 5))
		#plt.plot(IDdata_samples_vector[0, i, :])  # Fissa la prima dimensione (1) e scorri lungo i rimanenti due assi
		baseline = np.mean(IDdata_samples_vector[0, i, :100])  # Calcola la media delle prime 100 entrate della WF
		min_value = np.min(IDdata_samples_vector[0, i, :])
		#plt.axhline(baseline, color='r', linestyle='--', label='Media delle prime N entrate della WF') 
		#plt.axhline(min_value, color='g', linestyle='--', label='Valore minimo')  # Linea del minimo della WF
		#plt.title(f'PMT {i+1}')
		#plt.xlabel('DAQ window time [ns]')
		#plt.ylabel('Waveform')
		#plt.grid(True)
		integrated_charge = 0
		
		

		if(baseline-min_value>20):	# this should be improved with a dedicated "cluster" algorithm!
			diff_baseline_max = baseline - min_value

			# Find the index of array[0, i, :] such that the difference between baseline and array[0, i, :] is (baseline-max)/N
			for idx, val in enumerate(IDdata_samples_vector[0, i, :]):
				if abs(baseline - val) > abs(diff_baseline_max) / 5:	# rise time determined via "CFD"
					#print(idx)
					
					if(IDdata_channelID[0,i] % 2 ==0):	# only high-gain channels should be considered
						break
						
					#if(IDdata_GCUID[0,i] == 17):
					#	plt.plot(IDdata_samples_vector[0, i, :], color='blue', alpha=0.2) 
					
					#plt.scatter(idx, val, color='red', label='trigger time')
					
					fired_PMTs += 1
					
					integrated_charge += np.sum(IDdata_samples_vector[0, i, idx:]-baseline)
					
					coordinate = find_PMT_Coordinates("OSIRIS_cable_map.conf", IDdata_GCUID[0,i], IDdata_channelID[0,i])
					
					#print(IDdata_GCUID[0,i], IDdata_channelID[0,i],coordinate[0],coordinate[1],coordinate[2])
					
					with open(sys.argv[2], 'a') as file:
						file.write(f"{extracted_string}\t{extracted_date}\t{j}\t{integrated_charge}\t{i}\t{idx}\t{trgNsec[0]}\t{IDdata_channelID[0, i]}\t{IDdata_GCUID[0, i]}\t{coordinate[0]}\t{coordinate[1]}\t{coordinate[2]}\t{int(IDdata_channelID.shape[1]/2)}\n")
					
					coordinate_x = np.append(coordinate_x, float(coordinate[0]))
					coordinate_y = np.append(coordinate_y, float(coordinate[1]))
					coordinate_z = np.append(coordinate_z, float(coordinate[2]))					
					somme_differenze = np.append(somme_differenze, integrated_charge)

					#print(IDdata_channelID[0, i],IDdata_GCUID[0, i],coordinate[0],coordinate[1],coordinate[2])
					
					#print(coordinate[0], coordinate[1], integrated_charge)
					
					
					break
	
	#ax_2.scatter(coordinate_x, coordinate_y, c=somme_differenze, cmap='viridis')
	
	#print("Event " + str(j) + " , fired PMTs : " + str(fired_PMTs))
	
	#if(len(coordinate_z) > 50):
	
	#	coordinate_theta = np.sign(coordinate_y)*np.arccos(coordinate_x / np.sqrt(coordinate_x**2+coordinate_y**2))
		
	#	fig, ax = plt.subplots(1,2, figsize=(12, 5))
	#	plt.colorbar(ax[0].scatter(coordinate_x, coordinate_y, c=-somme_differenze, s = 60, cmap='Reds'))		
	#	ax[0].scatter(all_PMTs_x, all_PMTs_y, marker='o', edgecolor='black', alpha = 0.3, s = 60, facecolor='none')
	#	ax[0].set_xlim(-4200,4200)
	#	ax[0].set_ylim(-4200,4200)	
	#	ax[0].set_xlabel('PMT $x$ [mm]')
	#	ax[0].set_ylabel('PMT $y$ [mm]')
	#	ax[0].set_title('')

	#	plt.colorbar(ax[1].scatter(coordinate_theta, coordinate_z, c=-somme_differenze, s = 60, cmap='Reds'))			
	#	ax[1].scatter(all_PMTs_theta, all_PMTs_z, marker='o', edgecolor='black', alpha = 0.3, s = 60, facecolor='none')		
	#	ax[1].set_xlabel('PMT $\\theta$ [rad]')
	#	ax[1].set_ylabel('PMT $z$ [mm]')
	#	ax[1].set_title('')
	#	plt.tight_layout()
	#	plt.show()
	
				
		#plt.close()
		
			#else:
			#	continue
		
			#break  # Esci dal ciclo esterno se l'indice è stato trovato
		#plt.show()

#plt.show()

