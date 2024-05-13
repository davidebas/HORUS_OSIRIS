import numpy as np
import matplotlib.pyplot as plt

class Waveform:
	def __init__(self, samples, threshold_method, baseline_entries, cleaning):
		self.samples = samples
		self.min_value = np.min(samples)
		self.baseline = np.mean(samples[:baseline_entries])  # Average of earliest baseline_entries of WF
		self.std_dev_baseline = np.std(samples[:baseline_entries])
		self.threshold_method = threshold_method
		self.cleaning = cleaning

	def clean_data(self, samples):
		freq_domain = np.fft.fft(self.samples)
		freqs = np.fft.fftfreq(len(self.samples))

		# Plotting to visualize the process
		#fig, axs = plt.subplots(3, 1, figsize=(10, 12))
		#axs[0].plot(self.samples)

		filtered_freq_domain = np.where(abs(freqs) < 0.10, freq_domain, 0)
		self.samples = np.sqrt( (np.fft.ifft(filtered_freq_domain).real)**2+ (np.fft.ifft(filtered_freq_domain).imag)**2)

		#axs[0].set_title('Original Waveform')
		#axs[1].plot(freqs, np.abs(freq_domain))
		#axs[1].set_yscale('log')
		#axs[1].set_xlim(0,0.2)
		#axs[1].set_title('Frequency Domain')
		#axs[2].plot(self.samples)
		#axs[2].set_title('Deconvoluted Waveform')
		#plt.tight_layout()
		#plt.show()
		
		return self.samples

	def check_threshold(self, method):
		
		if(self.cleaning == 'true'):
			self.samples = self.clean_data(self.samples)
		
		if method == 'baseline':
			threshold_absolute = 20
			return self.baseline - self.min_value > threshold_absolute
		elif method == 'std_dev':
			threshold_std_dev = 5
			return (self.baseline - self.min_value) > threshold_std_dev * self.std_dev_baseline
		else:
			raise ValueError("Threshold finder method not valid. Please use 'baseline', 'std_dev', or 'deconvolution'.")


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
					end_idx = min(idx + 100, len(self.samples))
					integrated_charge += np.sum(self.samples[idx:end_idx] - self.baseline)
					threshold_crossing_index = idx  # Memorize the index associated to rise time
					break 

			#plt.figure()
			#plt.plot(self.samples)
			#plt.xlabel('Time')
			#plt.ylabel('Amplitude')
			#plt.title('Waveform')
			#plt.show()

		return fired_PMTs, integrated_charge, threshold_crossing_index
