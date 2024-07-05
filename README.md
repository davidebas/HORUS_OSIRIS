# HORUS: Hybrid OSIRIS data Reconstruction and Understanding Software
- v0.2
- 2024 July 05
- Authors: Davide Basilico (davide.basilico@mi.infn.it) and Marco Beretta (marco.beretta@mi.infn.it)
- Basic event analyzer for OSIRIS (waveform analysis + event reconstruction) with additional tools for the determination of background levels.

1) Conversion from binary format to ROOTfile: eb2root parser (credits to Kai Loo)

./eb2root -i datFile.dat -o RootfileEBParser.root

2) Waveform Analysis:

python3 WaveformAnalyzer.py RootfileEBParser.root Output_WaveformAnalyzer.txt

3) Event reconstruction:

- Position reco based on charge barycenter
- TOF calculation
- Applying the muon veto

usage: EventReconstruction.py [-h] [-i INPUTFILE] [-o OUTPUTFILE] [-m MUON] [-a ALL] [-d DIR]

options:
  -h, --help            show this help message and exit
  -i INPUTFILE, --InputFile INPUTFILE
                        Input file name
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        Output file name; default: rec-input_file.root
  -m MUON, --muon MUON  Apply muon veto; default: true
  -a ALL, --All ALL     Analyse all the data file with a specific pattern; default: false
  -d DIR, --Dir DIR     Directory containing files to be analyzed; default: false
  
The output ROOTfile from the EventReconstruction step includes:
- reconstructed position (coordinates x,y,z) for each event
- number of overall fired PMTs (Fired_PMTs) and OD fired PMTs (OD_fired) for each event
- reconstructed energy: Charge, Charge normalized to the number of Live PMTs, Charge normalized to the number of Live PMTs for ID PMTs (Charge, Charge_Norm, Charge_Norm_ID) for each event
- reconstructed energy for each couple of subsequent events, after the calibration (Charge/3800): Energy_Prompt and Energy_Delayed
- Time of flight for each fired PMT: TOF
- Rise Time for each fired PMT (WF_RiseTime) and the subtracted Rise Time - TOF (WF_RiseTime_diff)
- number of channels for each event (Shape_Ch)
