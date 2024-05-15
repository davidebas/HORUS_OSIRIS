# OSIRIS RECO + COINCINDENCE analysis  v.0.1.0
# Authors: Davide Basilico davide.basilico@mi.infn.it, Marco Beretta marco.beretta@mi.infn.it

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import uproot
import os
from tqdm import tqdm
from tqdm.auto import tqdm
tqdm.pandas()

base_file = sys.argv[1]
if(base_file[-9] =="_"):
    base_number = int(base_file[-5])
    base_filename = base_file[:-5] 
    file_list = [f"{base_filename}{i}.txt" for i in range(0, base_number + 1)]
    all_results = pd.DataFrame()

else:
    file_list = [base_file]

### Standard event reco

for file in file_list:
    print("open file: " + file)
    data_unclean = pd.read_csv(file, delimiter='\t')
    data = data_unclean.replace([np.inf, -np.inf], np.nan).dropna(subset=['charge'])

    # Group by index and calculation of the metrics
    grouped = data.groupby('index')
    total_charge = -grouped['charge'].sum()
    average_live_pmts = grouped['LivePMTs'].mean()
    trgTime = grouped['trgTime'].mean()
    trgTime_diff = trgTime.diff().shift(-1)
    normalized_charge = total_charge / average_live_pmts
    Fired_PMTs = grouped.size()

    # Center of charge calculation
    weights = data['charge'].abs()
    x_CM = (data['x_PMT'] * weights).groupby(data['index']).sum() / weights.groupby(data['index']).sum()
    y_CM = (data['y_PMT'] * weights).groupby(data['index']).sum() / weights.groupby(data['index']).sum()
    z_CM = (data['z_PMT'] * weights).groupby(data['index']).sum() / weights.groupby(data['index']).sum()

    data = data.merge(x_CM.rename('x_CM'), on='index')
    data = data.merge(y_CM.rename('y_CM'), on='index')
    data = data.merge(z_CM.rename('z_CM'), on='index')

    # Time of flight calculation (TOF)
    data['TOF'] = np.sqrt(
        (data['x_CM'] - data['x_PMT'])**2 +
        (data['y_CM'] - data['y_PMT'])**2 +
        (data['z_CM'] - data['z_PMT'])**2
    ) / 1000 / 1.55 / 3e8 * 1e9  # TOF in nanoseconds

    # Difference between rise time and TOF
    data['WF_RiseTime_diff'] = data['WF_RiseTime'] - data['TOF']

    # Dataframe creation
    results = pd.DataFrame({
        'Charge': total_charge,
        'Fired_PMTs': Fired_PMTs,
        'Charge_Norm': normalized_charge,
        'Energy': normalized_charge / 3650,
        'x_CM': x_CM,
        'y_CM': y_CM,
        'z_CM': z_CM,
        'trgTime': trgTime,
        'trgTime_diff': trgTime.diff().shift(-1),
        'TOF': data.groupby('index')['TOF'].apply(list),
        'WF_RiseTime': data.groupby('index')['WF_RiseTime'].apply(list),
        'WF_RiseTime_diff': data.groupby('index')['WF_RiseTime_diff'].apply(list)
    })

    if(base_file[-9] =="_"):
        all_results = pd.concat([all_results, results], ignore_index=True)
    else:
        all_results = results

### Coincidence anlysis

results['trgTime_diff'] = results['trgTime'].diff()

bismuto = []
polonio = []

Po_tau = 237E-6 # mean life time
hmTau = 5 #how many Tau
EB_min = 0.0 #Bismuth cuts
EB_max =3.5
EP_min = 0.6 #Polonium cuts
EP_max = 1.3
r_cut = 5000  #no cut for now
#accidental_offset = np.random.uniform(low=0.01, high=1.0, size=1000) #10e-3 modify this value to evalute the accidental backgrounds
offset = 0

results_filtered = all_results[all_results['Energy'] < EB_max]
coincidenze = []
num_coi = []

#for offset in tqdm(accidental_offset, desc="Progress"):

for j in range(len(all_results)): #tqdm( , desc="Processando eventi principali")
    for i in range(j + 1, len(results_filtered)):
        trg_diff = (results_filtered['trgTime'].iloc[i])- results_filtered['trgTime'].iloc[j]
        if trg_diff >= Po_tau * hmTau + offset or trg_diff < offset:
            break        
        pos_diff = np.sqrt(
            (results_filtered['x_CM'].iloc[i] - results_filtered['x_CM'].iloc[j])**2 +
            (results_filtered['y_CM'].iloc[i] - results_filtered['x_CM'].iloc[j])**2 +
            (results_filtered['z_CM'].iloc[i] - results_filtered['z_CM'].iloc[j])**2
        )
        
        if pos_diff < r_cut:
            if EB_min < results_filtered['Energy'].iloc[j] < EB_max:
                if EP_min <= results_filtered['Energy'].iloc[i] <= EP_max:
                    coincidenze.append((results_filtered['Energy'].iloc[j], results_filtered['Energy'].iloc[i]))

### Plotting part

bismuto = [b for b, p in coincidenze]  
polonio = [p for b, p in coincidenze]  

num_coincidenze = len(coincidenze)
#num_coi.append(len(coincidenze))

plt.figure(figsize=(10, 6)) 

plt.hist(bismuto, bins=100, alpha=0.5, label='Bismuth', color='blue')
plt.hist(polonio, bins=20, alpha=0.5, label='Polonium', color='red')

if(base_file[-9] =="_"):
    plt.text(0.4, 0.9, f'Run: {base_filename}*.txt\nCoincidenze: {num_coincidenze}', transform=plt.gca().transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.7, edgecolor='black'))
else:
    plt.text(0.4, 0.9, f'Run: {base_file}\nCoincidenze: {num_coincidenze}', transform=plt.gca().transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.7, edgecolor='black'))


plt.xlabel('Energy (MeV)')  
plt.ylabel('Counts')  
plt.title('Bi - Po spectra')  
plt.legend(loc='upper left')  
plt.yscale('log')

### Saving part

save_name = base_file[:-4] 
percorso_file = os.path.join("Bi-Po_plot", f"{save_name}.png")  
plt.savefig(percorso_file)

plt.show()
#plt.hist(num_coi, bins=100, alpha=0.5, label='Accidentals', color='blue')
#plt.show()
