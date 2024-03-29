# -*- coding: utf-8 -*-
# OSIRIS EVENT RECONSTRUCTION
# Author: Davide Basilico davide.basilico@mi.infn.it  2024 March 27

import numpy as np
import uproot
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
import matplotlib.colors as mcolors
plt.rcParams.update({'font.size': 12})

def y_fit(x, a, b,c,d):			# function I used to fit the firedPMT-charge scatter plot -- can be changed, not really needed
    return (a*x*x+b*x+c+np.exp(d*x)) 
    
occurrences = {}
sums = {}
occurrences_HG = {}
sums_HG = {}
sums_HG_CN = {}
x_CM_per_numero = {}
y_CM_per_numero = {}
z_CM_per_numero = {}
Charge_Norm = {}

#### READING THE WF ANALYZER OUTPUT
data = np.genfromtxt(sys.argv[1], skip_header=0)  #in case, skip the first line which should be an header and change skip_header=0
NameFile = data[:, 0]
Date = data[:, 1]
index = data[:, 2]
charge = data[:, 3]
WF_RiseTime = data[:, 5]
trgNsec = data[:, 6]
ID_channel = data[:, 7]
xx = data[:, 9]
yy = data[:, 10]
zz = data[:, 11]
LivePMTs = data[:, 12]

for numero in np.unique(index):
    indici = index == numero
    x_CM = np.sum(charge[indici] * xx[indici]) / np.sum(charge[indici])
    x_CM_per_numero[numero] = x_CM
    y_CM = np.sum(charge[indici] * yy[indici]) / np.sum(charge[indici])
    y_CM_per_numero[numero] = y_CM
    z_CM = np.sum(charge[indici] * zz[indici]) / np.sum(charge[indici])
    z_CM_per_numero[numero] = z_CM
    print(numero, np.sum(charge[indici]/LivePMTs[indici]), x_CM,y_CM,z_CM)

WF_RiseTime_diff_aligned = np.zeros_like(WF_RiseTime)
WF_RiseTime_diff = np.zeros_like(WF_RiseTime)
TOF = np.zeros_like(WF_RiseTime)

Multiplicity_Cut = 2
unique_indices, counts = np.unique(index, return_counts=True)
groups_to_process = unique_indices[counts > Multiplicity_Cut]

TOF = np.zeros_like(WF_RiseTime)
WF_RiseTime_diff = np.zeros_like(WF_RiseTime)
WF_RiseTime_diff_aligned = np.zeros_like(WF_RiseTime)
min_WF_RiseTime_per_evento = np.zeros_like(WF_RiseTime)

#### Calculation of rise time
for numero in groups_to_process:
    indici = index == numero
    if np.sum(indici) > Multiplicity_Cut:   
        min_WF_RiseTime_per_evento[indici] = np.min(WF_RiseTime[indici])

#### Calculation of TOF and rise time differences (only for selected events)
for i in range(len(WF_RiseTime)):
    if index[i] in groups_to_process:
        TOF[i] = np.sqrt((x_CM_per_numero[index[i]] - xx[i])**2 + (y_CM_per_numero[index[i]] - yy[i])**2 + (z_CM_per_numero[index[i]] - zz[i])**2) / 1000 / (1.55) / 3.0e8 * 1.0e9

        WF_RiseTime_diff[i] = WF_RiseTime[i] - TOF[i]
        WF_RiseTime_diff_aligned[i] = WF_RiseTime_diff[i] - min_WF_RiseTime_per_evento[i]
  
#### Calculation of number of fired PMTs and integrated charge for each event
for i in range(len(index)):
    idx = int(index[i])
    occurrences[idx] = occurrences.get(idx, 0) + 1
    sums[idx] = sums.get(idx, 0) + charge[i]
    
    if ID_channel[i] % 2 != 0:
        occurrences_HG[idx] = occurrences_HG.get(idx, 0) + 1
        sums_HG[idx] = sums_HG.get(idx, 0) + charge[i]
        sums_HG_CN[idx] = sums_HG_CN.get(idx, 0) + charge[i]/LivePMTs[i]


for idx in sums:
    sums[idx] /= occurrences[idx]

occurrences_counts = list(occurrences.values())
sums_values = list(sums.values())
sums_values = [-x for x in sums_values]

occurrences_counts_HG = list(occurrences_HG.values())
sums_values_HG = list(sums_HG.values())
sums_values_HG = [-x for x in sums_values_HG]

sums_values_HG_CN = list(sums_HG_CN.values())
sums_values_HG_CN = [-x for x in sums_values_HG_CN]


x_CM_array = np.array(list(x_CM_per_numero.values())) / 1000
y_CM_array = np.array(list(y_CM_per_numero.values())) / 1000
z_CM_array = np.array(list(z_CM_per_numero.values())) / 1000


#### PLOTTING

# Reconstructed energy spectra + scatter plot firedPMTs:charge
plt.figure(figsize=(10, 5))

plt.subplot(1, 3, 1)
plt.hist(occurrences_counts_HG, bins=max(occurrences_counts_HG)+1, color='blue', edgecolor='blue', linewidth=1.2)
plt.xlabel('"Fired" Waveforms')
plt.ylabel('events')
plt.yscale('log')
plt.title('"Fired" PMTs histogram')

plt.subplot(1, 3, 2)
plt.hist(sums_values_HG, bins=1000, range=(0,1e6), color='blue', edgecolor='blue', linewidth=1.2)
plt.xlabel('Charge [A.U.]')
plt.ylabel('Events')
plt.yscale('log')
plt.title('Charge histogram, high gain channels')

print(len(occurrences_counts_HG))
print(len(sums_values_HG))

plt.subplot(1, 3, 3)
range_xy = [[0,80], [0,1e6]]
plt.hist2d(occurrences_counts_HG, sums_values_HG, bins = 80, range=range_xy, cmap='jet')
plt.xlabel('Fired PMTs')
plt.ylabel('Charge [A.U.]')


# Reconstructed coordinate: scatter plots
plt.figure(figsize=(10, 5))

x_bins = np.linspace(x_CM_array.min(), x_CM_array.max(), 200)
y_bins = np.linspace(y_CM_array.min(), y_CM_array.max(), 200)
z_bins = np.linspace(y_CM_array.min(), y_CM_array.max(), 200)

x_range = [np.percentile(x_CM_array, 4), np.percentile(x_CM_array, 96)]
y_range = [np.percentile(y_CM_array, 4), np.percentile(y_CM_array, 96)]
z_range = [np.percentile(z_CM_array, 4), np.percentile(z_CM_array, 96)]

r_CM_array = np.sqrt(x_CM_array**2 + y_CM_array**2 + z_CM_array**2)

range_xy = [[-4.3, 4.3], [-4.3, 4.3]]

plt.subplot(1, 3, 1)
plt.hist2d(x_CM_array, y_CM_array, bins = 100, range=range_xy, cmap='jet')
plt.xlabel('x [m]')
plt.ylabel('y [m]')

plt.subplot(1, 3, 2)
plt.xlabel('y [m]')
plt.ylabel('z [m]')
plt.hist2d(y_CM_array, z_CM_array, bins = 100, range=range_xy, cmap='jet')

plt.subplot(1, 3, 3)
plt.xlabel('x [m]')
plt.ylabel('z [m]')
plt.hist2d(x_CM_array, z_CM_array, bins = 100, range=range_xy, cmap='jet')

plt.tight_layout()

# Reconstructed coordinate: 1D Histos
plt.figure(figsize=(10,5))


plt.subplot(2, 2, 1)
plt.hist(x_CM_array, bins = 200, range=(-3,3) )
plt.xlabel('x [m]')
plt.subplot(2, 2, 2)
plt.xlabel('y [m]')
plt.hist(y_CM_array, bins = 200, range=(-3,3))
plt.subplot(2, 2, 3)
plt.hist(z_CM_array, bins = 200, range=(-3,3))
plt.xlabel('z [m]')
plt.subplot(2, 2, 4)
plt.hist(r_CM_array, bins = 200, range=(0,5e0))
plt.xlabel('r [m]')
#plt.show()

plt.figure(figsize=(13,7))


plt.subplot(2,3,1)
plt.hist(WF_RiseTime_diff_aligned, bins = 200, range = (-100,700) )
plt.xlabel('t scint aligned [ns]')
plt.ylabel('entries')
plt.subplot(2,3,4)
plt.hist(WF_RiseTime_diff_aligned, bins = 200, range = (-100,700) )
plt.yscale('log')
plt.xlabel('t scint aligned [ns]')
plt.ylabel('entries [log scale]')
plt.subplot(2,3,2)
plt.hist(WF_RiseTime_diff, bins = 200, range = (100,600) , color = "purple")
plt.xlabel('t scint [ns]')
plt.ylabel('entries')
plt.subplot(2,3,5)
plt.hist(WF_RiseTime_diff, bins = 200, range = (100,600) , color = "purple")
plt.yscale('log')
plt.xlabel('t scint [ns]')
plt.ylabel('entries [log scale]')

plt.subplot(2,3,3)
plt.hist(TOF, bins = 200, range = (-5,40) , color = "darkorange")
plt.xlabel('TOF [ns]')
plt.yscale('log')
plt.ylabel('entries')

plt.subplot(2,3,6)
plt.hist(min_WF_RiseTime_per_evento, bins = 200, range = (100,600) , color = "seagreen" )
plt.xlabel('minimum time rise (for each event)')
plt.ylabel('entries')
plt.tight_layout()
#plt.show()

Charge_Norm = np.array(Charge_Norm)

print(Charge_Norm)
print(Charge_Norm.size)
print(Charge_Norm.shape)

output_filename = "output.root"

# Apre il file ROOT in modalità scrittura
with uproot.recreate(output_filename) as f:
    # Crea un nuovo file
    f["RecEvents"] = {
    	"index" : np.unique(index),
        "x": x_CM_array,
        "y": y_CM_array,
        "z": z_CM_array,
        "r": r_CM_array,
        "FiredPMTs": occurrences_counts_HG,
        "Charge": np.array(sums_values_HG)/1000,
        "Charge_Norm" : np.array(sums_values_HG_CN)/1000
    }

print("ROOT file:", output_filename)
