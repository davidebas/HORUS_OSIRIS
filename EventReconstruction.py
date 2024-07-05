import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import uproot
from tqdm import tqdm
from tqdm.auto import tqdm
tqdm.pandas()
import os
import argparse
from argparse import RawTextHelpFormatter

c = 3.0e8
n_LS = 1.55

def parse_arguments():
    prs = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    prs.add_argument("-i", "--InputFile", default="null", help="Input file name")
    prs.add_argument("-o", "--OutputFile", default="null", help="Output file name; default: rec-input_file.root")
    prs.add_argument("-m", "--muon", default="true", help="Apply muon veto; default: true")
    prs.add_argument("-a", "--All", default="false", help="Analyse all the data file with a specific pattern; default: false")
    prs.add_argument("-d", "--Dir", default="false", help="Directory containing files to be analyzed; default: false")
    prs.add_argument("-Muon_Veto_Threshold", "--Threshold_OD_Fired", default=5, help="Muon veto threshold for OD multiplicity; default: 5")    
    return prs.parse_args()

def load_data(input_file, all_files, folder_path):
    df_list = []

    if input_file != "null" and all_files == "false":
        data_unclean = pd.read_csv(input_file, delimiter='\t')
        print("Open only: ", input_file)
        return data_unclean

    if all_files != "false":
    	folder = os.path.dirname(input_file)
    	prefix = os.path.basename(input_file).rsplit('_', 1)[0]  # Ottieni il prefisso dal nome del file di input

    	files = [f for f in os.listdir(folder) if f.startswith(prefix)]
    	files.sort()
    	for i, filename in enumerate(files):
    	    file_path = os.path.join(folder, filename)
    	    df = pd.read_csv(file_path, delimiter='\t')
    	    if i != 0:
    	        df["index"] = df["index"] + df_list[-1]["index"].max() + 1
    	    df_list.append(df)
    	    print("Loading: ", filename)
    	    print("Events from ", df["index"].min(), " to ", df["index"].max())
    	data_unclean = pd.concat(df_list, ignore_index=True)
    	return data_unclean


    if folder_path != "false":
        files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]
        files.sort()
        for i, filename in enumerate(files):
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path, delimiter='\t')
            if i != 0:
                df["index"] = df["index"] + df_list[-1]["index"].max() + 1
            df_list.append(df)
            print(f"Open: {filename}")
            print(df["index"].min(), df["index"].max())
        data_unclean = pd.concat(df_list, ignore_index=True)
        print(data_unclean["index"].max())
        return data_unclean

    return None

def process_data(data):
    data = data.replace([np.inf, -np.inf], np.nan).dropna(subset=['charge'])
    grouped = data.groupby('index')
    total_charge = -grouped['charge'].sum()
    average_live_pmts = grouped['LivePMTs'].mean()
    trgTime = grouped['trgTime'].mean()
    trgTime_diff = trgTime.diff().shift(-1)
    normalized_charge = total_charge / average_live_pmts
    Fired_PMTs = grouped.size()
    OD_fired = grouped['OD'].sum()
    Shape_Ch = grouped['Shape_Ch'].max()

    data_OD = data[data['OD'] == 1]
    grouped_OD = data_OD.groupby('index')
    total_charge_OD = -grouped_OD['charge'].sum()
    average_live_pmts_OD = grouped_OD['LivePMTs'].mean()
    normalized_charge_OD = total_charge_OD / average_live_pmts_OD

    data_ID = data[data['OD'] == 0]
    grouped_ID = data_ID.groupby('index')
    total_charge_ID = -grouped_ID['charge'].sum()
    average_live_pmts_ID = grouped_ID['LivePMTs'].mean()
    normalized_charge_ID = total_charge_ID / average_live_pmts_ID

    weights = data['charge'].abs()

    print("-- Position reco based on charge barycenter")
    x_CM = (data['x_PMT'] * weights).groupby(data['index']).sum() / weights.groupby(data['index']).sum()
    y_CM = (data['y_PMT'] * weights).groupby(data['index']).sum() / weights.groupby(data['index']).sum()
    z_CM = (data['z_PMT'] * weights).groupby(data['index']).sum() / weights.groupby(data['index']).sum()

    data = data.merge(x_CM.rename('x_CM'), on='index')
    data = data.merge(y_CM.rename('y_CM'), on='index')
    data = data.merge(z_CM.rename('z_CM'), on='index')

    print("-- TOF calculation")
    data['TOF'] = np.sqrt((data['x_CM'] - data['x_PMT']) ** 2 +
                          (data['y_CM'] - data['y_PMT']) ** 2 +
                          (data['z_CM'] - data['z_PMT']) ** 2) / 1000 / n_LS / c * 1.0e9

    TOF = data.groupby('index')['TOF'].apply(list)
    data['WF_RiseTime_diff'] = data['WF_RiseTime'] - data['TOF']
    WF_RiseTime_diff_per_event = data.groupby('index')['WF_RiseTime_diff'].apply(list)
    WF_RiseTime_per_event = data.groupby('index')['WF_RiseTime'].apply(list)
    mean_WF_RiseTime_per_event = data.groupby('index')['WF_RiseTime'].mean()

    trgTime_aligned = grouped['trgTime'].mean() - (mean_WF_RiseTime_per_event - 250) * 1e-9
    trgTime_diff_aligned = trgTime_aligned.diff().shift(-1)

    results = pd.DataFrame({
        'Charge': total_charge,
        'Fired_PMTs': Fired_PMTs,
        'Charge_Norm': normalized_charge,
        'Charge_Norm_OD': normalized_charge_OD,
        'Charge_Norm_ID': normalized_charge_ID,
        'x_CM': x_CM,
        'y_CM': y_CM,
        'z_CM': z_CM,
        'trgTime': trgTime,
        'trgTime_diff': trgTime_diff,
        'trgTime_aligned': trgTime_aligned,
        'trgTime_diff_aligned': trgTime_diff_aligned,
        'TOF': TOF,
        'WF_RiseTime': WF_RiseTime_per_event,
        'WF_RiseTime_diff': WF_RiseTime_diff_per_event,
        'OD_fired': OD_fired,
        'Shape_Ch': Shape_Ch
    })
    return results

def apply_muon_veto(results,Threshold_OD_Fired):
    print("-- Applying the muon veto")
    od_mult = results['OD_fired'] >= Threshold_OD_Fired
    events_to_remove = results[od_mult]

    time_threshold = 20E-6
    for i, (index, row) in enumerate(tqdm(events_to_remove.iterrows(), desc="\tApplying time cut: ")):
        if i == 0:
            print("\n\n\tIt could take time. This is the total number of iterations: ", len(events_to_remove.index.to_numpy()), ". Be patient...\n")
        time_diff_condition = (results['trgTime'] - row['trgTime']).abs() < time_threshold
        results = results.loc[~time_diff_condition]

    results = results.loc[~od_mult]
    return results

def calculate_energies(results):
    print("-- Calculating Energy_Prompt and Energy_Delayed")
    results['Energy_Prompt'] = results['Charge_Norm'] / 3800.
    results['Energy_Delayed'] = results['Charge_Norm'].shift(-1) / 3800.
    results = results.dropna(subset=['Energy_Prompt'])
    return results

def save_results(results, output_file, input_file, all_files, folder_path):
    if output_file == "null":
        if all_files == "true" or folder_path != "false":
            output_filename = "rec-" + input_file[:-9] + ".root"
        else:
            output_filename = "rec-" + input_file[:-4] + ".root"
    else:
        output_filename = output_file

    with uproot.recreate(output_filename) as f:
        f["RecEvents"] = {
            "Index": results.index.to_numpy(),
            "x": results['x_CM'].to_numpy() if 'x_CM' in results else None,
            "y": results['y_CM'].to_numpy() if 'y_CM' in results else None,
            "z": results['z_CM'].to_numpy() if 'z_CM' in results else None,
            "Fired_PMTs": results['Fired_PMTs'].to_numpy(),
            "Charge": results['Charge'],
            "Charge_Norm": results['Charge_Norm'],
            "Charge_Norm_OD": results['Charge_Norm_OD'],
            "Charge_Norm_ID": results['Charge_Norm_ID'],
            'trgTime': results['trgTime'],
            'trgTime_diff': results['trgTime_diff'],
            'trgTime_aligned': results['trgTime_aligned'],
            'trgTime_diff_aligned': results['trgTime_diff_aligned'],
            "TOF": results['TOF'],
            'WF_RiseTime': results['WF_RiseTime'],
            'WF_RiseTime_diff': results['WF_RiseTime_diff'],
            'Energy_Prompt': results['Energy_Prompt'].to_numpy(),
            'Energy_Delayed': results['Energy_Delayed'].to_numpy(),
            'OD_fired': results['OD_fired'],
            'Shape_Ch': results['Shape_Ch']
        }
    print(f"Output ROOT file created: {output_filename}")

def main():
    args = parse_arguments()

    data_unclean = load_data(args.InputFile, args.All, args.Dir)
    if data_unclean is None:
        print("No data loaded.")
        return

    print("### Event Reconstruction ###")
    results = process_data(data_unclean)

    if args.muon == "true":
        results = apply_muon_veto(results,args.Threshold_OD_Fired)

    results = calculate_energies(results)
    save_results(results, args.OutputFile, args.InputFile, args.All, args.Dir)

if __name__ == "__main__":
    main()

