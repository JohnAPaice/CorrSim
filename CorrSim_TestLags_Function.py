## Creates model Lags

import CorrSim_Functions as corrfunc
import importlib
import math
import matplotlib.pyplot as plt
import numpy as np
import os

importlib.reload(corrfunc)

'''
	-Test Model Lags-
In this code, you can quickly test what shape your lags create.

--- Step 1: ---
Choose whether to define your lags using time or phase lags; time_or_phase = "time"/"phase"

--- Step 2: ---
Define the overall lag. This defines what the lag should be in any undefined sections.

--- Step 3: ---
Set the lags dependant upon frequency, using either phase or time space.
This is done by 'drawing' the 'shape' of the lags using a variety of distributions.
Phase lags are described in semi-log space. Time lags are described in log-log space.

The following functions are used:
corrfunc.lag_section_phase(distribution, start_freq, close_freq, start_lag, close_lag=0, \
	extra_freq=0, extra_lag=0)
corrfunc.lag_section_time(distribution, start_freq, close_freq, start_lag, close_lag=0, \
	extra_freq=0, extra_lag=0)

	-Inputs:-

To define a lag, call corrfunc.lag_section_time or corrfunc.lag_section_time,
and define the distribution and the required parameters.

NOTE: Make sure to define your sections IN ORDER OF INCREASING FREQUENCY.
If you don't, it will create incorrect lags!


The following lag distributions are available:
	- Constant			= Constant lag, defined by start_lag.
							Requires: start_freq, close_freq, start_lag

	- Linear			= When time_or_phase="Phase": Linear in semi-log space.
							When time_or_phase="Time": Linear in log-log space.
							Requires: start_freq, close_freq, start_lag, close_lag

	- Power				= Exponentially increasing lag (Linear in normal space).
							Requires: start_freq, close_freq, start_lag, close_lag

	- Polynomial		= When time_or_phase="Phase": Polynomial in semi-log space.
							When time_or_phase="Time": Polynomial in log-log space.
							Works by solving a second-order polynomial for the three coordinates given.
							Requires: start_freq, close_freq, start_lag, close_lag, extra_freq, extra_lag

	- Constant_Time		= ! Note: Only for corrfunc.lag_section_phase !
							Creates phase lags that give a constant time (set by Start Lag)
							Requires: start_freq, close_freq, start_lag

	- Constant_Phase	= ! Note: Only for corrfunc.lag_section_phase_time !
							Creates time lags that give a constant phase (set by Start Lag)
							Requires: start_freq, close_freq, start_lag


The output is seven numbers for each section. These are:
	- start		= Starting frequency of the section
	- close		= Ending frequency of the section
	- A			= Constant (see below)
	- B			= Constant (see below)
	- C			= Constant (see below)
	- P			= Constant (see below)
	- log		= Boolean; defines whether or not distributions have to be calculated using logs (1=Yes, 0=No)

	A, B, C, P are constants in the expression (Ax^2 + Bx^P + C)

'''


def TestLags(output_dir, fileprefix, obs_length, time_res, time_or_phase, overall_lag, model_lag_array):


	##					Make Output Directory
	##	--------------------------------------------------------

	if not os.path.exists(output_dir):
	    os.makedirs(output_dir)


	## ===== Calculations =====

	print("TESTING LAG MODELS\n")

	print("Creating model lags...")

	## --- Define frequencies

	freq_min	= 1 / obs_length						## Minimum Frequency
	freq_max	= 1 / (2*time_res)						## Maximum Frequency
	freq_res	= 1 / obs_length						## Frequency Resolution

	frequencies	= np.arange(0, freq_max+(freq_res*0.1), freq_res)


	## --- Calculate model lags

	model_lags = np.zeros(len(frequencies)) + overall_lag			## Define the basic model lags

	section_pointer = 0												## Create a variable that points to a section
	for i in range(len(model_lags)):								## For each frequency bin...
		if frequencies[i] < model_lag_array[section_pointer, 0]:	## If we're not yet in the current section,
			next													## Continue until we reach one.
		elif frequencies[i] > model_lag_array[section_pointer, 1]:	## If we're past the current section.
			section_pointer = section_pointer + 1					## Increment to the next section

		if section_pointer == len(model_lag_array[:,1]):			## If we're out of defined sections.
			break													## Break.
		elif frequencies[i] >= model_lag_array[section_pointer, 0] and frequencies[i] <= model_lag_array[section_pointer, 1]:	## If we're in a section,
			if (model_lag_array[section_pointer, 6] == 1):													## If we need to use logs,
				model_lags[i] = model_lag_array[section_pointer, 2]*math.log(frequencies[i],10)**2 +	\
					model_lag_array[section_pointer, 3]*math.log(frequencies[i],10)**model_lag_array[section_pointer, 5] + \
					model_lag_array[section_pointer, 4]														## Calculate the lag based on the section's parameters.
			else:																					## If we don't need to use logs,
				model_lags[i] = model_lag_array[section_pointer, 2]*frequencies[i]**2 + \
					model_lag_array[section_pointer, 3]*frequencies[i]**model_lag_array[section_pointer, 5] + \
					model_lag_array[section_pointer, 4]														## Calculate the lag based on the section's parameters.

	if time_or_phase == 'time':
		print("Using: Time Lags")

		model_time_lags		= model_lags
		model_phase_lags	= np.array([])
		for j in range(len(model_lags)):
			model_phase_lags = np.append(model_phase_lags, model_lags[j] * 2 * np.pi * frequencies[j])

	elif time_or_phase == 'phase':
		print("Using: Phase Lags")

		model_time_lags		= np.array([])
		model_phase_lags	= model_lags

		model_time_lags = np.append(model_time_lags, np.nan)
		for j in range(1, len(model_lags)):
			model_time_lags = np.append(model_time_lags, model_lags[j] / (2 * np.pi * frequencies[j]))

	# model_time_lags_inv = [-x for x in model_time_lags]



	## Rebin the Data - for faster plotting
	freqs_rebinned			= corrfunc.log_rebin(frequencies, 400)
	frequencies_binned		= freqs_rebinned[0]
	cuts					= freqs_rebinned[1]
	model_phase_lags_binned	= corrfunc.log_rebin_average(model_phase_lags, cuts, frequencies_binned)
	model_time_lags_binned	= corrfunc.log_rebin_average(model_time_lags, cuts, frequencies_binned)

	frequencies_binned		= np.asarray(frequencies_binned)
	model_phase_lags_binned	= np.asarray(model_phase_lags_binned)
	model_time_lags_binned	= np.asarray(model_time_lags_binned)



	### ===== Plot Model Lags =====

	print("Plotting...")

	fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6,6))

	## --- Phase Lags ---
	ax = axs[0]
	ax.set_xscale('log')

	ax.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
	ax.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
	ax.yaxis.set_major_formatter(plt.FuncFormatter(corrfunc.format_func))
	ax.set_ylim([-3.8, 3.8])

	ax.plot(frequencies_binned, model_phase_lags_binned-(4*np.pi),	lw=2, zorder=3, linestyle=':',  color="red", label = "-4π")
	ax.plot(frequencies_binned, model_phase_lags_binned-(2*np.pi),	lw=2, zorder=3, linestyle='--', color="firebrick", label = "-2π")
	ax.plot(frequencies_binned, model_phase_lags_binned, 			lw=2, zorder=3, linestyle='-',  color="black", label = "±0")
	ax.plot(frequencies_binned, model_phase_lags_binned+(2*np.pi),	lw=2, zorder=3, linestyle='--', color="mediumblue", label = "+2π")
	ax.plot(frequencies_binned, model_phase_lags_binned+(4*np.pi),	lw=2, zorder=3, linestyle=':',  color="dodgerblue", label = "+4π")

	ax.axhline(+np.pi, linestyle="-", linewidth=1.5, color="darkgrey", zorder=2)
	ax.axhline(-np.pi, linestyle="-", linewidth=1.5, color="darkgrey", zorder=2)

	ax.legend()

	ax.set_ylabel("Phase Lag (Radians)")
	ax.grid(linestyle='--')

	if time_or_phase == "phase":
		for spine in ['top', 'right', 'bottom', 'left']:
			ax.spines[spine].set_linewidth(2)
			ax.spines[spine].set_color("red")
			ax.spines[spine].set_zorder(10)
	else:
		ax.spines["bottom"].set_linewidth(2)
		ax.spines["bottom"].set_color("red")
		ax.spines["bottom"].set_zorder(10)


	## --- Time Lags ---
	ax = axs[1]
	ax.set_xscale('log')
	ax.set_yscale('log', nonpositive="clip")

	ax.plot(frequencies_binned,  model_time_lags_binned, linestyle="-",  lw=2, color="firebrick", label = "+ve Lags")
	ax.plot(frequencies_binned, -model_time_lags_binned, linestyle="--", lw=2, color="blue", label = "-ve Lags")

	ax.legend()

	ax.set_ylabel("Time Lag (s)")
	ax.set_xlabel("Frequency (Hz)")
	ax.grid(linestyle='--')

	plt.subplots_adjust(bottom=0.1, left=0.15, top=0.98, right=0.98, hspace=0)

	if time_or_phase == "time":
		for spine in ['top', 'right', 'bottom', 'left']:
			ax.spines[spine].set_linewidth(2)
			ax.spines[spine].set_color("red")
			ax.spines[spine].set_zorder(10)
	else:
		ax.spines["top"].set_linewidth(2)
		ax.spines["top"].set_color("red")
		ax.spines["top"].set_zorder(10)

	## --- Export ---

	fileOut = "./" + output_dir + "/" + fileprefix + "Test Lags.png"
	plt.savefig(fileOut)

	plt.close('all')

	print("Plotted model lags.")
