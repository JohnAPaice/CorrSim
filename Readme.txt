###  CorrSim - A Correlated Timing Observation Simulator  ###

>                      By John A Paice                      <
> With invaluable help from Ranjeev Misra and Poshak Gandhi <

Thank you for downloading CorrSim! I hope you find it useful in your projects.

The main aim of this code is to:
1: Help in planning multiwavelength timing observations, and
2: Aid in the interpretation of these observations.

If this project was indeed useful to you in your research, please include the following acknowledgement:
“This work made use of the \texttt{CorrSim} software.”
With a footnote linking to the appropriate gitlab page:
https://gitlab.com/astro_johnapaice/CorrSim/

To cite the work, please cite Paice et al. (in Prep)


## CorrSim’s function:

CorrSim takes a variety of input parameters, and with them, creates two lightcurves - Series A and Series B. These lightcurves are related, as defined by given coherence and time lag parameters. This is meant to simulate a multiwavelength observation.

CorrSim’s main role compared to other simulation software is that it gives you fine control over how these lightcurves are related to each other, in a Fourier and Timing sense.

Firstly, it allows you to define the power spectra of the lightcurves, either through a broken power law or a series of Lorentzians. Through this, you can also set how coherent the lightcurves are with each other - i.e. how much they are correlated with each other as a function of Fourier frequency.

Secondly, it allows you to define the phase lags or time lags between the lightcurves, with respect to Fourier frequency. These are defined by essentially ‘plotting’ the lags on a graph.


## Dependencies:

CorrSim requires the following python modules:
copy
importlib
math
matplotlib
numpy
os
pandas
scipy
stingray
time
tkinter
warnings


## How to use CorrSim:

There are two ways to use CorrSim. One is with the GUI, while the other is by calling the function directly.

Directions for both methods are found later on in this readme - Ctrl-F 'GUI Instructions' for the GUI, and 'Function Instructions' for the function.


## Tips on speeding up/slowing down CorrSim:

(Note that this is only tested on my machine - so take this with a grain of salt, and try some testing yourself if you're interested.
This is also all just for the calculations - time to plot the graphs afterwards is not considered)

The time CorrSim takes to run is linearly correlated with the number of bins (So 1e5 bins will take ~10 times longer than 1e4 bins)

The following is the fastest way to run the code (aside from changing the number of bins):
- No noise sources
- Define Power Spectra using a Broken Power Law
- No phase/time lags defined (Or only define low-frequencies)

Compared to code run with the above, here is how each change affects the runtime:
- Red Noise (Band A);		x2
- Red Noise (Band B);		x3
- Defining 1 Lorentzian;		x2.5
- Defining 6 Lorentzians;		x3
- Defining high-frequency lags;	x2

Combining the above has the effect of adding the effects of the multipliers; e.g. Defining 6 Lorentzians and also the high-frequency lags has a combined multiplier of x4.

CorrSim_GUI can also be run from within ipython; but be warned that it can slow down if it's closed and re-opened multiple times. Restarting ipython can solve this.


## Warnings about use:

- The code automatically creates text files containing the Fourier models and lightcurves. For large numbers of bins (>1e6), these can be very large (>100 MB). Be aware of the size of the files you are creating!

- Particularly long or high-resolution observations can be taxing on your system, and may make CorrSim take a long time to run. Press 'Calculate Bins' to make sure you don't go over the soft 1e6 warning. If your number of bins has a very large factor, it can also drastically slow down the calculations; you can optimise the number of bins by checking 'Optimise Bins', which adjusts your observation length so that your number of bins is the nearest ‘7-smooth’ number, i.e. a number with no factor higher than 7.





### How to use CorrSim with the GUI (GUI Instructions):

Run the CorrSim GUI using:

> python3 CorrSim_GUI.py
(Or your Python 3 version of choice. CorrSim is not compatible with previous versions of Python.)

Input parameters are set within this GUI. Parameter names can be hovered over with the mouse to bring up an explanatory tooltip, and are also described at the end of this file in ‘Glossary: Input Parameters’. You can use File -> Load/Save to load/save parameters. An example file (Example_Input_Parameters.txt) is provided.

To help you create Power Spectra and Phase/Time lags, you can test the models using the 'Test Power Spectra' and 'Test Lags' buttons. These will just plot your inputs, so you can check what you’re telling the program to create.

By running the program, it will automatically create files containing:
- A log file containing your parameters
- The simulated power spectra (F_rms Normalised) as a text file
- The simulated lightcurves as a text file
- The simulated lightcurves as a .png

Optionally, the code can also give you the following:
- The Model CCF as both a text file and a png
- The CCF calculated from the simulated lightcurves as both a text file and a plot.
- The Model Fourier data (Power Spectra, Coherence, and Phase/Time lags) as both a text file and a plot.
- The Fourier data calculated from the simulated lightcurves as both a text file and a plot.

It will also give you extra plots if multiple boxes are checked:
- A plot of the CCF calculated from the simulated lightcurves, with the model overlaid on top (If both Model and Calculate CCF boxes are checked).
- A plot of the Fourier data calculated from the simulated lightcurves, with the models overlaid on top (If both Model and Calculate Fourier boxes are checked).
- A large compilation plot of the inputs and outputs (If all four plotting boxes are checked).


## Glossary of Parameters:


—-- Observation & Source Parameters —--

— Observation Parameters: —

Observation Length (s) (Number, >0)
Length of simulated observation (seconds). Particularly relevant for the number of bins, which is the biggest factor in the size of the output files and the speed of CorrSim.

Time Resolution (s) (Number, >0)
Time resolution of simulated observation (seconds). Particularly relevant for the number of bins, which is the biggest factor in the size of the output files and the speed of CorrSim.

Optimise Bins (y/n)
Optional; The observation length is changed so that the number of bins a 7-smooth number (i.e. a number with no factors greater than 7). This makes the code run relatively fast.


— Source Parameters: —

Mean Counts/s (Number, >0)
Mean count rate in Series A and B (Counts/s)

F_rms (Number, >0)
Fractional RMS in Series A and B. Defines how variable the lightcurves are.
For example: for the X-ray lightcurve of X-ray Binaries, F_rms = ~0.3 in the hard state, and ~0.04 in the Soft State.
For the Optical lightcurve of X-ray Binaries, F_rms = ~0.1.


— Noise Parameters: —

Apply Red Noise (y/n)
Applies ‘red noise’ to either Series A, B, both, or neither.
Red noise is also known as ‘Brownian noise’, and is proportional to frequency^-2.
CorrSim applies red noise to the power spectra before they’re converted to lightcurves, using the methods of Timmer & Koenig 1995 (1995A&A...300..707T):
For each Fourier Frequency, two Gaussian distributed random numbers are drawn.
These are multiplied by sqrt(0.5*(power_spectrum)).
The result is the real and imaginary part of the Fourier transformed data.
Red noise is dependent on the Fractional RMS and Slope.

Red Noise F_rms (Number, >0)
Fractional RMS of the red noise; the larger the number, the greater the noise.

Red Noise Slope (Number)
Defines the dependence of the noise on frequency (Noise is proportional to f^(slope)).


Apply Poisson Noise (y/n)
Applies ‘poisson noise’ to either Series A, B, both, or neither.
Poisson noise, also known as ‘shot noise’, comes from the random nature of emitted photons.
CorrSim applies poisson noise to the lightcurves by drawing a random number from a poisson distribution for each bin, with the mean centred on the flux in that bin.


Apply Readout Noise (y/n)
Applies ‘readout noise’ to either Series A, B, both, or neither.
Readout noise comes from fluctuations in reading out charge from a CCD.
This is dependant on ‘read_noise_A’ and ‘read_noise_B’, defined below.

Readout Noise (Number, >0)
Readout noise from a CCD in electrons. Only applicable if Apply Readout Noise is ticked in either A or B.


Apply Scintillation Noise (y/n)
Applies ‘scintillation noise’ to either Series A, B, both, or neither.
This noise comes from atmospheric effects, and is thus only applicable to simulations from ground-based observatories.
Scintillation noise is dependant on a telescopes diameter, altitude, and exposure time, a target’s altitude, the height of atmospheric turbulence, and some empirical value. Default values are for the NTT at La Silla, Chile, with values taken from Osborn et al. 2015 (DOI: 10.1093/mnras/stv1400).

Telescope Diameter (Number, >0)
Diameter of the observing telescope (metres).
An increase in this value decreases the scintillation noise.

Telescope Altitude (Number, >0)
Altitude of the observing telescope (metres)
An increase in this value decreases the scintillation noise.

Exposure Time (Number, < time_resolution)
Exposure time of the observing telescope. This is different from time resolution, as it is does NOT include the time to read out information from the CCD. 
An increase in this value decreases the scintillation noise.

Turbulence Height (Number, >0)
Height of turbulence in the atmosphere (metres). A typical value here is around 8000.
An increase in this value increases the scintillation noise.

Target Altitude (Number, 0-90)
Altitude of the source (degrees).
An increase in this value decreases the scintillation noise.

Empirical Value (Number, >0)
An empirical coefficient that varies depending on the site of the telescope. This is defined as C_Y in Osborn et al. 2015 (DOI: 10.1093/mnras/stv1400), where several sites are listed. The mean value is 1.5.
An increase in this value increases the scintillation noise.


--— Define Power Spectra --—

The power spectra can be defined either through a broken powerlaw model or a Lorentzian model.
CorrSim_TestPowerSpectra.py can be used to experiment and test inputs until it produces what you like!

Broken Powerlaw OR Lorentzians
Sets the type of power spectra to be used.
‘Broken Powerlaw’ uses a simple model of constant power that breaks at some frequency, and then becomes a power law component that decreases with increasing frequency. The coherence is described in the same way. 
‘Lorentzians’ allows you to define a series of Lorentzians that sum to the power spectrum for each series. Series B includes an extra parameter for each Lorentzian to define how coherent it is with Series A.


Parameters for the Broken Powerlaw model:


Power Spectra:

Index (Number)
Index of the Power Law component of the power spectra.

Break Freq. (Number, >0)
Break Frequency of the power law (the frequency at which it transitions from a constant to a power law).


Coherence:

Constant (Number, >0)
The fraction that Series B is coherent with Series A during the constant component of the coherence.

Index (Number)
Index of the Power Law component of the coherence.

Break Freq. (Number, >0)
Break Frequency of the coherence (the frequency at which it transitions from a constant to a power law)


Parameters for the Lorentzian model:

This array of boxes produces the Lorentzians. Up to six can be used in each series. Each Lorentzian is defined by:

L(f) = N/π (0.5*W)/((f-M)^2+(0.5*W)^2)

Where f is frequency, and N, W, and M are the Normalisation (Norm), Width, and Midpoint (Mid) respectively.
The Series B Lorentzians have an extra parameters, Coh, which defines how coherent each Lorentzian is with Series A.


—-- Define Lags —--

Time OR Phase
Lags can be defined either as time lags or phase lags.
Time lags are defined in log-log space, while phase lags are defined in semi-log space.

Overall Lag (Number)
The time/phase lag to use for all frequencies not otherwise defined.


Parameters for the lag model:

Lags are defined by defining several sections of different distributions. Think of it as 'drawing' the lag plot; define the distribution you want (e.g. Linear, Polynomial), the starting frequency, the ending frequency, and then any required lags:

(Unused):		This line isn't used.						
Const. Time:		Constant lag (L_1) in time between F_1 and F_2.			(Requires: Freq_1, Lag_1, Freq_2)
Const. Phase:		Constant lag (L_1) in phase between F_1 and F_2.			(Requires: Freq_1, Lag_1, Freq_2)
Linear:			Linear distribution between (F_1, L_1) and (F_2, L_2).		(Requires: Freq_1, Lag_1, Freq_2, Lag_2)
Power:			Exponentially increasing lag between (F_1, L_1) and (F_2, L_2).	(Requires: Freq_1, Lag_1, Freq_2, Lag_2)
Polynomial:		Solves a second-order polynomial for the three coordinates given.	(Requires: Freq_1, Lag_1, Freq_2, Lag_2, Freq_3, Lag_3)

Note that the program needs 0 values for all unused boxes, and will often fail if not provided them.


—-- Plotting Parameters —--

Plot Red Noise (y/n)
Plots the model power spectra, the simulated red noise, and a summation of the two for each band.


Plot Model CCF (y/n)
Plots the model Cross-Correlation Function (A time-domain representation of how the two lightcurves correlate as a function of time lag. 1 = perfect correlation, 0 = no correlation, -1 = perfect anti-correlation).

Calculate CCF (y/n)
Calculates and plots the Cross-Correlation Function (CCF) of the simulated lightcurves. This CCF is made by splitting up the lightcurves into segments (defined by Segment Size (s)), running CCF analysis on each of them, and then averaging the result.

Segment Size (s) (Number)
Define the length (in seconds) the lightcurves are cut up into before the Cross-Correlation Function is calculated. Must be less than Observation Length.

Binning (Bins) (Integer)
Optional; if >1, average the lightcurves over this many bins before the Cross-Correlation Function is calculated. This can be done to speed up calculations.


Plot Fourier Models (y/n)	
Plot the input models of the Power Spectra, Coherence, Phase Lags, and Time Lags, all dependant on Fourier frequency.
Power Spectra are representations of how much signals vary.
Coherence is how much the two signals correlate (0 < c(f) < 1)
Phase Lags are a representation of how much each frequency is offset in phase.
Time Lags are a translation of phase lags to the time domain.

Calculate Fourier (y/n)
Calculates and plots the power spectra, coherence, and phase/time lags between the two series. It accomplishes this by splitting up the lightcurves in a number of segments (size given by Segment Size (bins)), carrying out the analysis on each segment, and then averaging the result. These will be different from the inputted values, due to the effects of noise and finite sampling.

Segment Size (bins) (Integer)
Define the length (in bins) the lightcurves are cut up into before the Fourier analysis is carried out. Ideally, this should be a power of two.

Rebinning (Number, >0)
Define the amount of logarithmic rebinning to do when calculating and plotting the Fourier analysis. This goes straight to the Stingray software, which uses the convention of Rebinning-1; i.e. a value of 0.1 will rebin logarithmically with a factor of 1.1.

Reference Freq. (Hz) (Number)
This is used when calculating the phase lags. Phase lags are inherently constrained between +/-pi. If the actual phase lags are outside of this range, e.g. between pi and 2pi, then analysis will show that they are between -pi and 0 radians due to the cyclical nature of sine waves. A ‘reference frequency’ is thus defined to be the frequency at which the code calibrates the rest of the phase lags; this should be the frequency at which the measured phase lag is between +/-pi.

Remove White Noise (y/n)
This is used for plotting the recovered power spectra and coherence.
Note that this is a VERY simple approximation - it assumes that the last point in the power spectrum is white-noise dominated and removes 99% of that value from the power spectra.
This is ONLY a good approximation when the spectrum is white-noise dominated at the highest Fourier frequency.


—-- Final Parameters —--

Output Directory (String)
Defines the output directory; e.g. ‘CorrSim_Outputs’ will save all outputs to the folder ./CorrSim_Outputs. This folder will be created if it doesn't already exist.

File prefix (String)
Defines a prefix for the filenames of all outputted files.





### How to use the CorrSim function directly (Function Instructions):

You can also use the CorrSim function directly by calling it in python. In a directory that has all the files, import it as:

> import CorrSim

And then call the function:

> CorrSim.CorrSim()

The file CorrSim.py has a brief description of the arguments, some of which differ slightly from the parameters in the GUI. All are fairly straight-forward, except for scin_noise_A, scin_noise_B, Lorentz_params_A, Lorentz_params_B, and model_lag_array. They are listed here in detail:

## Glossary of Arguments

    output_dir              (String)                            
Defines the output directory.              
(Default: "CorrSim_Outputs")

    fileprefix              (String)                            
Defines a prefix for all output files       
(Default: "CorrSim - ")

    obs_length              (Number, >0)                        
Length of observation in seconds            
(Default: 1000)

    time_res                (Number, >0)                        
Time resolution of observation in seconds   
(Default: 0.1)

    mean_counts_per_sec_A   (Number, >0)                        
Mean count rate in Series A                 
(Default: 1000)

    mean_counts_per_sec_B   (Number, >0)                        
Mean count rate in Series B                 
(Default: 5000)

    F_rms_A                 (Number, >0)                        
Fractional RMS in Series A                  
(Default: 0.3)

    F_rms_B                 (Number, >0)                        
Fractional RMS in Series B                  
(Default: 0.1)

    apply_red               ("A"/"B"/"AB"/"n")                  
Apply Red Noise to A/B/AB/None              
(Default: "n")

    F_rms_rednoise_A        (Number, >0)
Fractional RMS of the red noise; the larger the number, the greater the noise.
(Default: 0.2)

    F_rms_rednoise_B        (Number, >0)
Fractional RMS of the red noise; the larger the number, the greater the noise.
(Default: 0.03)

    rednoise_slope_A        (Number, >0)
Defines the dependence of the noise on frequency (Noise is proportional to f^(slope)).
(Default: -2)

    rednoise_slope_B        (Number, >0)
Defines the dependence of the noise on frequency (Noise is proportional to f^(slope)).
(Default: -2)

    apply_poisson           ("A"/"B"/"AB"/"n")                  
Apply Poisson Noise to A/B/AB/None          
(Default: "n")

    apply_scintillation     ("A"/"B"/"AB"/"n")                  
Apply Scintillation Noise to A/B/AB/None    
(Default: "n")

    apply_readout           ("A"/"B"/"AB"/"n")                  
Apply Readout Noise to A/B/AB/None          
(Default: "n")

    remove_whitenoise       (1/0)                               
Automatically remove white noise?           
(Default: 0)

    scin_noise_A            (Number, >0)                        
Scintillation Noise parameter for B. This is calculated by running CorrSim_MakeScintillation.py, which will prompt you for various parameters (explained in the GUI documentation above).
(Default: 0)

    scin_noise_B            (Number, >0)                        
Scintillation Noise parameter for B. This is calculated by running CorrSim_MakeScintillation.py, which will prompt you for various parameters (explained in the GUI documentation above).
(Default: 0)

    read_noise_A            (Number, >0)                        
Read Noise parameter for A, in units of electrons.                 
(Default: 0)

    read_noise_B            (Number, >0)                        
Read Noise parameter for B, in units of electrons.                 
(Default: 0)

    plot_model_ccf          (1/0)                               
Plot model CCF?                             
(Default: 0)

    calculate_ccf           (1/0)                               
Calculate CCF?                              
(Default: 0)

    ccf_segment             (Integer)                           
CCF Segment Size in seconds                 
(Default: 30)

    ccf_binning             (Number, >0)                        
CCF binning (bins)                          
(Default: 0)

    plot_model_fourier      (1/0)                               
Plot model Fourier products?                
(Default: 0)

    calculate_fourier       (1/0)                               
Calculate Fourier products?                 
(Default: 0)

    fourier_bins            (Number, >0)                        
Fourier Segment Size (bins)                 
(Default: 512)

    fourier_rebin           (Number, >0)                        
Fourier Rebinning                           
(Default: 0.3)

    reference_freq          (Number, >0)                        
Reference Frequency (Hz)                    
(Default: 1)

    power_model_type        ("broken_powerlaw"/"lorentzians")   
Which Power Model type to use               
(Default: "broken_powerlaw")

    ps_power                (Number)                            
Index of the Power Spectra                  
(Default: -1.5)

    break_freq              (Number, >0)                        
Break Frequency of the Power Spectra        
(Default: 1)

    coh_constant            (Number, >0)                        
Constant value of the Coherence             
(Default: 0.1)

    coh_power               (Number)                            
Index of the Coherence                      
(Default: -0.5)

    coh_break_freq          (Number, >0)                        
Break Frequency of the Coherence            
(Default: 0.5)

    Lorentz_params_A        (List)                              
Lorentzian Parameters for Series A. This must be a list divisible by 3. Each set of three is a Lorentzian, with the first number being the Normalisation, the second being the Width, and the third being the Midpoint.
(Default: [])

    Lorentz_params_B        (List)                              
Lorentzian Parameters for Series B. This must be a list divisible by 4. Each set of four is a Lorentzian, with the first number being the Normalisation, the second being the Width, the third being the Midpoint, and the fourth being the Lorentzian's coherence (between 0 and 1) with Series A.
(Default: [])

    time_or_phase           ("time"/"phase")                    
Define lags in Time or Phase                
(Default: "phase")

    overall_lag             (Number, >0)                        
Define overall lag                          
(Default: 0)

    model_lag_array         (ndarray)                       
Array for the lag parameters. This requires running CorrSim_MakeLags.py, which will prompt you for several values. When it has finished running, it will print an numpy array to the terminal - use this (or the variable it gives - model_lag_array) as the argument for model_lag_array.
(Default: np.zeros((6, 7)))

    write_data              (1/0)                               
Write the data to text files?               
(Default: "y")





## Acknowledgements

CorrSim uses the Stingray software package for some of its Fourier analysis, which can be found at https://github.com/StingraySoftware/stingray and is described in Huppenkothen et al. 2019 (DOI: 10.3847/1538-4357/ab258d). John A Paice thanks A Stevens and D Huppenkothen for help with the software. John A Paice also thanks D. Ashton for spectral timing help.
Some observational parameters given as examples were taken by GTC/HiPERCAM and NICER; we thank the teams of both of those telescopes, particularly our HiPERCAM observers Vik Dhillon and Stuart Littlefair. Full acknowledgements and analysis are described in Paice et al. 2019 (DOI: 10.1093/mnrasl/slz148