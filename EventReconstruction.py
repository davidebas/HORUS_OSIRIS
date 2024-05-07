import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import uproot
from tqdm import tqdm
from tqdm.auto import tqdm
tqdm.pandas()

c = 3.0e8
n_LS = 1.55

print("### Event Reconstruction ###")		

def weighted_center(group, coord):
    return (group[coord] * group['charge'].abs()).sum() / group['charge'].abs().sum()

plt.rcParams.update({'font.size': 12})

data_unclean = pd.read_csv(sys.argv[1], delimiter='\t')

data = data_unclean.replace([np.inf, -np.inf], np.nan).dropna(subset=['charge'])

grouped = data.groupby('index')
total_charge = -grouped['charge'].sum()
average_live_pmts = grouped['LivePMTs'].mean()
trgTime = grouped['trgTime'].mean()
trgTime_diff = trgTime.diff().shift(-1)
normalized_charge = total_charge / average_live_pmts
Fired_PMTs = grouped.size()

weights = data['charge'].abs()

print("-- Position reco based on charge barycenter")
x_CM = (data['x_PMT'] * weights).groupby(data['index']).sum() / weights.groupby(data['index']).sum()
y_CM = (data['y_PMT'] * weights).groupby(data['index']).sum() / weights.groupby(data['index']).sum()
z_CM = (data['z_PMT'] * weights).groupby(data['index']).sum() / weights.groupby(data['index']).sum()

data = data.merge(x_CM.rename('x_CM'), on='index')
data = data.merge(y_CM.rename('y_CM'), on='index')
data = data.merge(z_CM.rename('z_CM'), on='index')

print("-- TOF calculation")
data['TOF'] = np.sqrt((data['x_CM'] - data['x_PMT'])**2 +
                      (data['y_CM'] - data['y_PMT'])**2 +
                      (data['z_CM'] - data['z_PMT'])**2) / 1000 / n_LS / c * 1.0e9

data['WF_RiseTime_diff'] = data['WF_RiseTime'] - data['TOF']                      

TOF = data.groupby('index')['TOF'].apply(list)
WF_RiseTime_diff_per_event = data.groupby('index')['WF_RiseTime_diff'].apply(list)
WF_RiseTime_per_event = data.groupby('index')['WF_RiseTime'].apply(list)

print("-- Saving results into output rootfile")
results = pd.DataFrame({
    'Charge': total_charge,
    'Fired_PMTs': Fired_PMTs,        
    'Charge_Norm': normalized_charge,
    'x_CM': x_CM,
    'y_CM': y_CM,
    'z_CM': z_CM,
    'trgTime' : trgTime,
    'trgTime_diff' : trgTime_diff,
    'TOF' : TOF,
    'WF_RiseTime' : WF_RiseTime_per_event,
    'WF_RiseTime_diff' : WF_RiseTime_diff_per_event
})

#for idx, row in results.iterrows():
#    print(f"Event index: {idx}, integrated charge: {row['Charge']}")

output_filename = sys.argv[2]

with uproot.recreate(output_filename) as f:
    f["RecEvents"] = {
        "Index": results.index.to_numpy(),
        "x": results['x_CM'].to_numpy() if 'x_CM' in results else None,
        "y": results['y_CM'].to_numpy() if 'y_CM' in results else None,
        "z": results['z_CM'].to_numpy() if 'z_CM' in results else None,
        "Fired_PMTs": results['Fired_PMTs'].to_numpy(),
        "Charge": results['Charge'],
        "Charge_Norm": results['Charge_Norm'],
        'trgTime' : trgTime,
        'trgTime_diff' : trgTime_diff,
        "TOF": results['TOF'],
        'WF_RiseTime' : results['WF_RiseTime'],
        'WF_RiseTime_diff' : results['WF_RiseTime_diff']
    }

print(f"Output ROOT file created: {output_filename}")




#plt.figure(figsize=(10, 5))
#plt.hist(results['WF_RiseTime_diff'], bins=100, color='blue', edgecolor='blue', linewidth=1.2)
#plt.xlabel('Charge_Norm [A.U.]')
#plt.ylabel('Events')
#plt.yscale('log')
#plt.title('Charge_Norm Histogram')
#plt.show()


















