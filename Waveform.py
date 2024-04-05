# waveform_analyzer.py
import numpy as np

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
