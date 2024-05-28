# HORUS: Hybrid OSIRIS data Reconstruction and Understanding Software
Basic event analyzer for OSIRIS (waveform analysis + event reconstruction) with additional tools for the determination of background levels.

-- Waveform analysis:
python3 WaveformAnalyzer.py InputRootfileFromEBParser.root OutputWaveformAnalyzer.txt

-- Event reconstruction:
python3 EventReconstruction.py OutputWaveformAnalyzer.txt
