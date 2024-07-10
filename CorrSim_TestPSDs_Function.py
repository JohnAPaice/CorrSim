## Creates model Power Spectra

import CorrSim_Functions as corrfunc
import importlib
import math
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas

importlib.reload(corrfunc)

## ========== Fourier Parameters ==========

## ----- Power Spectrum Parameters -----

'''
	-Power Spectral Parameters-
Here, one defines the shape of the power spectra, and the coherence of the two bands.
There are currently two types of power spectral models implemented: broken_powerlaw, and lorentzian

	- broken_powerlaw:	Define a single broken power law for both the power spectra and the coherence.
						This is the simpler of the two models.

						The broken power law for both power spectra is the same shape, and are defined by a power
							and a break frequency (their differing normalisations will come later, based on F_rms).
						For the coherence, you must define a power, break frequency and a constant.

	- lorentzians:		Define lorentzians for both data sets.
						This is the more complex model, but the fine coherence control makes it more realistic.

						Lorentzians are defined using 3 parameters (Normalisation, Midpoint and Width)
						Band 2 Lorentzians take a 4th parameter, the coherence fraction (how coherent they are with Band 1)
						In both cases, this is done via a single 1-D list which is later parsed automatically.

Note that the overall normalisation of the power spectra is based on the F_rms values for each band.
	Therefore, the white noise parameter, as well as the 'Norm' parameter in the lorentzians, will only
	adjust the relative contributions of each component, and not the overall normalisation.

'''


def TestPSDs(output_dir, fileprefix, obs_length, time_res, power_model_type,
	ps_power, break_freq, coh_constant, coh_power, coh_break_freq, Lorentz_params_A, Lorentz_params_B):

	##					Make Output Directory
	##	--------------------------------------------------------

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	## ===== Calculations =====

	print("TESTING POWER SPECTRAL MODELS\n")

	print("Creating model Power Spectra...")

	## --- Define frequencies

	freq_min	= 1 / obs_length						## Minimum Frequency
	freq_max	= 1 / (2*time_res)						## Maximum Frequency
	freq_res	= 1 / obs_length						## Frequency Resolution

	frequencies	= np.arange(0, freq_max+(freq_res*0.1), freq_res)


	##		Power Spectrum Parameters
	##	----------------------------------

	## Type 1: Broken Powerlaw
	if power_model_type == "broken_powerlaw":

		print("Using: Broken Powerlaw")

		power_model = corrfunc.power_model_broken_powerlaw(frequencies, ps_power, break_freq, coh_constant, coh_power, coh_break_freq)

	## Type 2: Lorentzians
	elif power_model_type == "lorentzians":

		print("Using: Lorentzians")

		power_model = corrfunc.power_model_lorentzians(frequencies, Lorentz_params_A, Lorentz_params_B)

		Lorentzians_A				= power_model[4]
		Lorentzians_B_coh			= power_model[5]
		Lorentzians_B_inc			= power_model[6]
		num_lorentz_A				= power_model[7]
		num_lorentz_B				= power_model[8]


	power_spectrum_A		= power_model[0]
	power_spectrum_B_coh	= power_model[1]
	power_spectrum_B_inc	= power_model[2]

	power_spectrum_B		= power_spectrum_B_coh + power_spectrum_B_inc





	###
	### ===== Plotting =====
	###

	## Rebin the Data - for faster plotting
	freqs_rebinned			= corrfunc.log_rebin(frequencies, 400)
	frequencies_binned		= freqs_rebinned[0]
	cuts					= freqs_rebinned[1]
	power_spectrum_A_binned = corrfunc.log_rebin_average(power_spectrum_A, cuts, frequencies_binned)
	power_spectrum_B_binned = corrfunc.log_rebin_average(power_spectrum_B, cuts, frequencies_binned)


	### ==== Plot Model Power Spectra ===
	##  Plots the overall model power spectra:
	##  - In F_rms^2 power
	##  - F_rms^2 power mutliplied by frequency.

	print("Plotting overall models...")


	fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(7,7))

	## --- 1: Power
	ax = axs[0]
	ax.set_xscale('log')

	ax.plot(frequencies_binned,  power_spectrum_A_binned, linestyle="-", label="Band A", color="blue")
	ax.plot(frequencies_binned,  power_spectrum_B_binned, linestyle="-", label="Band B", color="red")

	ax.legend(loc="lower center")

	y_min = min([min((power_spectrum_B)[1:len(frequencies)]), min((power_spectrum_A)[1:len(frequencies)])])
	y_max = max([max((power_spectrum_B)[1:len(frequencies)]), max((power_spectrum_A)[1:len(frequencies)])])
	ax.set_xlim(frequencies[1]/2, max(frequencies)*2)
	ax.set_ylim(y_min/10, y_max*2)

	ax.set_xscale("log")
	ax.set_yscale("log")

	ax.set_xlabel('Frequency (Hz)')
	ax.set_ylabel('Frac. RMS Power')

	ax.grid(linestyle=":", zorder=5)


	## --- 2: Frequency * Power
	ax = axs[1]
	ax.set_xscale('log')

	ax.plot(frequencies_binned,  power_spectrum_A_binned*frequencies_binned, linestyle="-", label="Band A", color="blue")
	ax.plot(frequencies_binned,  power_spectrum_B_binned*frequencies_binned, linestyle="-", label="Band B", color="red")

	ax.legend(loc="lower center")

	y_min = min([min(frequencies_binned*power_spectrum_B_binned), min(frequencies_binned*power_spectrum_A_binned)])
	y_max = max([max(frequencies_binned*power_spectrum_B_binned), max(frequencies_binned*power_spectrum_A_binned)])
	ax.set_xlim(frequencies_binned[1]/2, max(frequencies_binned)*2)
	ax.set_ylim(y_min/10, y_max*2)

	ax.set_xscale("log")
	ax.set_yscale("log")

	ax.set_xlabel('Frequency (Hz)')
	ax.set_ylabel('Frequency * Frac. RMS Power')

	ax.grid(linestyle=":", zorder=5)

	plt.subplots_adjust(bottom=0.08, left=0.1, top=0.98, right=0.98, hspace=0)

	fileOut = "./" + output_dir + "/" + fileprefix + "Model Power Spectra.png"
	plt.savefig(fileOut)

	plt.close('all')



	###
	### ===== Plot Power Spectra Components =====
	##  Plots the overall model power spectra:
	##  - A (F_rms^2 power)					B (F_rms^2 power)
	##  - A (F_rms^2 power * frequency)		B (F_rms^2 power * frequency)

	print("Plotting Components...")

	fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=False, figsize=(10,7))

	## --- 1a: Power A
	ax = axs[0,0]

	ax.set_xlabel('Frequency (Hz)')
	ax.set_ylabel('Frac. RMS Power')

	if power_model_type == "broken_powerlaw":

		ax.plot(frequencies_binned,  power_spectrum_A_binned, linestyle="-", color="blue", alpha=0.5)

	elif power_model_type == "lorentzians":

		for i in range(num_lorentz_A):
			Lorentzians_A_binned = corrfunc.log_rebin_average(Lorentzians_A[:,i], cuts, frequencies_binned)
			ax.plot(frequencies_binned, Lorentzians_A_binned, linestyle = '--', color="blue", alpha=0.5)

	ax.plot(frequencies_binned, power_spectrum_A_binned, linewidth = 3, linestyle = '-', color="blue")

	line1 = mlines.Line2D([], [], color='blue',		linestyle='-',	label='A (Overall)')
	line2 = mlines.Line2D([], [], color='blue',		linestyle='--',	label='A (Components)')

	if power_model_type == "broken_powerlaw":
		ax.legend(handles=[line1], loc="lower center")
	elif power_model_type == "lorentzians":
		ax.legend(handles=[line1, line2], loc="lower center")

	y_min = min([min((power_spectrum_B)[1:len(frequencies)]), min((power_spectrum_A)[1:len(frequencies)])])
	y_max = max([max((power_spectrum_B)[1:len(frequencies)]), max((power_spectrum_A)[1:len(frequencies)])])
	ax.set_xlim(frequencies[1]/2, max(frequencies)*2)
	ax.set_ylim(y_min/10, y_max*2)

	print(y_min)
	print(y_max)

	ax.set_xscale("log")
	ax.set_yscale("log")

	ax.grid(linestyle=":", zorder=5, linewidth=0.5)

	print("Plotted Series A Components (1)...")


	## --- 1b: Power B
	ax = axs[0,1]

	ax.set_xlabel('Frequency (Hz)')

	if power_model_type == "broken_powerlaw":

		power_spectrum_B_coh_binned = corrfunc.log_rebin_average(power_spectrum_B_coh, cuts, frequencies_binned)
		power_spectrum_B_inc_binned = corrfunc.log_rebin_average(power_spectrum_B_inc, cuts, frequencies_binned)

		ax.plot(frequencies_binned, power_spectrum_B_coh_binned, linestyle = '-.', alpha=0.5, color="purple")
		ax.plot(frequencies_binned, power_spectrum_B_inc_binned, linestyle = ':',  alpha=0.5, color="red")

	elif power_model_type == "lorentzians":
		for i in range(num_lorentz_B):

			Lorentzians_B_coh_binned = corrfunc.log_rebin_average(Lorentzians_B_coh[:,i], cuts, frequencies_binned)
			Lorentzians_B_inc_binned = corrfunc.log_rebin_average(Lorentzians_B_inc[:,i], cuts, frequencies_binned)

			ax.plot(frequencies_binned, Lorentzians_B_coh_binned,linestyle = '-.', alpha=0.5, color="purple")
			ax.plot(frequencies_binned, Lorentzians_B_inc_binned,linestyle = ':',  alpha=0.5, color="red")

	ax.plot(frequencies_binned, power_spectrum_B_binned, linewidth = 3, linestyle = '-', color="red")

	line4 = mlines.Line2D([], [], color='red',		linestyle='-',	label='B (Overall)')
	line5 = mlines.Line2D([], [], color='purple',	linestyle='-.',	label='B (Coherent)')
	line6 = mlines.Line2D([], [], color='red',		linestyle=':',	label='B (Incoherent)')

	ax.legend(handles=[line4, line5, line6], loc="lower center")

	ax.set_xlim(frequencies[1]/2, max(frequencies)*2)
	ax.set_ylim(y_min/10, y_max*2)

	print(y_min)
	print(y_max)

	ax.set_xscale("log")
	ax.set_yscale("log")

	ax.grid(linestyle=":", zorder=5, linewidth=0.5)

	print("Plotted Series B Components (1)...")


	## --- 2a: Frequency * Power A
	ax = axs[1,0]

	ax.set_xlabel('Frequency (Hz)')
	ax.set_ylabel('Frequency * Frac. RMS Power')

	if power_model_type == "broken_powerlaw":

		ax.plot(frequencies_binned,  frequencies_binned*power_spectrum_A_binned, linestyle="-", color="blue", alpha=0.5)

	elif power_model_type == "lorentzians":

		for i in range(num_lorentz_A):
			Lorentzians_A_binned = corrfunc.log_rebin_average(Lorentzians_A[:,i], cuts, frequencies_binned)
			ax.plot(frequencies_binned, frequencies_binned*Lorentzians_A_binned, linestyle = '--', color="blue", alpha=0.5)

	ax.plot(frequencies_binned, frequencies_binned*power_spectrum_A_binned, linewidth = 3, linestyle = '-', color="blue")

	line1 = mlines.Line2D([], [], color='blue',		linestyle='-',	label='A (Overall)')
	line2 = mlines.Line2D([], [], color='blue',		linestyle='--',	label='A (Components)')

	if power_model_type == "broken_powerlaw":
		ax.legend(handles=[line1], loc="lower center")
	elif power_model_type == "lorentzians":
		ax.legend(handles=[line1, line2], loc="lower center")

	y_min = min([min(frequencies_binned*power_spectrum_B_binned), min(frequencies_binned*power_spectrum_A_binned)])
	y_max = max([max(frequencies_binned*power_spectrum_B_binned), max(frequencies_binned*power_spectrum_A_binned)])
	ax.set_xlim(frequencies_binned[1]/2, max(frequencies_binned)*2)
	ax.set_ylim(y_min/10, y_max*2)

	print(y_min)
	print(y_max)

	ax.set_xscale("log")
	ax.set_yscale("log")

	ax.grid(linestyle=":", zorder=5, linewidth=0.5)

	print("Plotted Series A Components (2)...")


	## --- 2b: Frequency * Power B
	ax = axs[1,1]

	ax.set_xlabel('Frequency (Hz)')

	if power_model_type == "broken_powerlaw":

		ax.plot(frequencies_binned, frequencies_binned*power_spectrum_B_coh_binned, linestyle = '-.', alpha=0.5, color="purple")
		ax.plot(frequencies_binned, frequencies_binned*power_spectrum_B_inc_binned, linestyle = ':',  alpha=0.5, color="red")

	elif power_model_type == "lorentzians":

		for i in range(num_lorentz_B):
			Lorentzians_B_coh_binned = corrfunc.log_rebin_average(Lorentzians_B_coh[:,i], cuts, frequencies_binned)
			Lorentzians_B_inc_binned = corrfunc.log_rebin_average(Lorentzians_B_inc[:,i], cuts, frequencies_binned)

			ax.plot(frequencies_binned, frequencies_binned*Lorentzians_B_coh_binned,linestyle = '-.', alpha=0.5, color="purple")
			ax.plot(frequencies_binned, frequencies_binned*Lorentzians_B_inc_binned,linestyle = ':',  alpha=0.5, color="red")

	ax.plot(frequencies_binned, frequencies_binned*power_spectrum_B_binned, linewidth = 3, linestyle = '-', color="red")

	line4 = mlines.Line2D([], [], color='red',		linestyle='-',	label='B (Overall)')
	line5 = mlines.Line2D([], [], color='purple',	linestyle='-.',	label='B (Coherent)')
	line6 = mlines.Line2D([], [], color='red',		linestyle=':',	label='B (Incoherent)')

	ax.legend(handles=[line4, line5, line6], loc="lower center")

	ax.set_xlim(frequencies_binned[1]/2, max(frequencies_binned)*2)
	ax.set_ylim(y_min/10, y_max*2)

	print(y_min)
	print(y_max)

	ax.set_xscale("log")
	ax.set_yscale("log")

	ax.grid(linestyle=":", zorder=5, linewidth=0.5)

	print("Plotted Series B Components (2)...")


	plt.subplots_adjust(bottom=0.08, left=0.08, top=0.98, right=0.98, hspace=0, wspace=0)

	fileOut = "./" + output_dir + "/" + fileprefix + "Power Spectrum Components.png"
	plt.savefig(fileOut)#, bbox_inches='tight')

	plt.close('all')

	print("Plotted model Power Spectra.")
