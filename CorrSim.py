## Main code for CorrSim. Use alongside CorrSim_Functions.

import CorrSim_Functions as corrfunc
import importlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
import pandas as pd
import scipy
import time
import warnings

importlib.reload(corrfunc)

## Ignore specific Stingray warning
warnings.filterwarnings("ignore", message="SIMON says: Errorbars on cross spectra are not thoroughly tested. Please report any inconsistencies.")

def CorrSim(output_dir="CorrSim_Outputs", fileprefix="CorrSim - ", \
            obs_length=1000, time_res=0.1, mean_counts_per_sec_A=1000, \
            mean_counts_per_sec_B=5000, F_rms_A=0.3, F_rms_B=0.1, \
            apply_red="n", F_rms_rednoise_A=0.2, F_rms_rednoise_B=0.03, \
            rednoise_slope_A=-2, rednoise_slope_B=-2, apply_poisson="n", \
            apply_readout="n", read_noise_A=0, read_noise_B=0, \
            apply_scintillation="n", scin_noise_A=0, scin_noise_B=0, \
            plot_rednoise=0, \
            plot_model_ccf=0, calculate_ccf=0, ccf_segment=30, ccf_binning=0, \
            plot_model_fourier=0, calculate_fourier=0, fourier_bins=512, \
            fourier_rebin=0.3, reference_freq=1, remove_whitenoise=0, \
            power_model_type="broken_powerlaw", ps_power=-1.5, break_freq=1, \
            coh_constant=0.1, coh_power=-0.5, coh_break_freq=0.5, \
            Lorentz_params_A=[], Lorentz_params_B=[], time_or_phase="phase", \
            overall_lag=0, model_lag_array=np.zeros((6, 7)), write_data=1):
    """
    Simulate two lightcurves based on a series of fourier parameters, noise
    sources, coherence models, and lags.

    Hello, and welcome to CorrSim!

    This code is primarily intended to be used through the CorrSim GUI;
    you can run that in this folder using:
    > python3 CorrSim_xGUI.py
    (Or your python3 version of choice)

    If you wish to call the function directly, then in a python session
    in this folder:
    > import CorrSim as CorrMain
    > CorrMain.CorrSim()
    The arguments are described below.
    For full descriptions, see the glossary in the Readme file.

    Parameters
    ----------
        output_dir : str
            Defines the output directory                (Default: "CorrSim_Outputs")
        fileprefix : str
            Defines a prefix for all output files       (Default: "CorrSim - ")
        obs_length : float (>0)
            Length of observation in seconds            (Default: 1000)
        time_res : float (>0)
            Time resolution of observation in seconds   (Default: 0.1)
        mean_counts_per_sec_A : float (>0)
            Mean count rate in Series A                 (Default: 1000)
        mean_counts_per_sec_B : float (>0)
            Mean count rate in Series B                 (Default: 5000)
        F_rms_A : float (>0)
            Fractional RMS in Series A                  (Default: 0.3)
        F_rms_B : float (>0)
            Fractional RMS in Series B                  (Default: 0.1)
        apply_red : "A"/"B"/"AB"/"n"
            Apply Red Noise to A/B/AB/None              (Default: "n")
        F_rms_rednoise_A : float (>0)
            Fractional RMS of the Red Noise in Band A   (Default: 0.2)
        F_rms_rednoise_B : float (>0)
            Fractional RMS of the Red Noise in Band B   (Default: 0.03)
        rednoise_slope_A : float (>0)
            Slope of the red noise in Band A            (Default: 2)
        rednoise_slope_B : float (>0)
            Slope of the red noise in Band B            (Default: 2)
        apply_poisson : "A"/"B"/"AB"/"n"
            Apply Poisson Noise to A/B/AB/None          (Default: "n")
        apply_readout : "A"/"B"/"AB"/"n"
            Apply Readout Noise to A/B/AB/None          (Default: "n")
        read_noise_A : float (>0)
            Read Noise parameter for A                  (Default: 0)
        read_noise_B : float (>0)
            Read Noise parameter for B                  (Default: 0)
        apply_scintillation : "A"/"B"/"AB"/"n"
            Apply Scintillation Noise to A/B/AB/None    (Default: "n")
        scin_noise_A : float (>0)
            Scintillation Noise parameter for A         (Default: 0)
        scin_noise_B : float (>0)
            Scintillation Noise parameter for B         (Default: 0)
        plot_rednoise : 1/0
            Plot the input rednoise?                    (Default: 0)
        plot_model_ccf : 1/0
            Plot model CCF?                             (Default: 0)
        calculate_ccf : 1/0
            Calculate CCF?                              (Default: 0)
        ccf_segment : float (>0)
            CCF Segment Size in seconds                 (Default: 30)
        ccf_binning : float (>0)
            CCF binning (bins)                          (Default: 0)
        plot_model_fourier : 1/0
            Plot model Fourier products?                (Default: 0)
        calculate_fourier : 1/0
            Calculate Fourier products?                 (Default: 0)
        fourier_bins : float (>0)
            Fourier Segment Size (bins)                 (Default: 512)
        fourier_rebin : float (>0)
            Fourier Rebinning                           (Default: 0.3)
        reference_freq : float (>0)
            Reference Frequency (Hz)                    (Default: 1)
        remove_whitenoise : 1/0
            Automatically remove white noise?           (Default: 0)
        power_model_type : "broken_powerlaw"/"lorentzians"
            Which Power Model type to use               (Default: "broken_powerlaw")
        ps_power : float
            Index of the Power Spectra                  (Default: -1.5)
        break_freq : float (>0)
            Break Frequency of the Power Spectra        (Default: 1)
        coh_constant : float (>0)
            Constant value of the Coherence             (Default: 0.1)
        coh_power : float
            Index of the Coherence                      (Default: -0.5)
        coh_break_freq : float (>0)
            Break Frequency of the Coherence            (Default: 0.5)
        Lorentz_params_A : list
            Lorentzian Parameters for Series A          (Default: [])
        Lorentz_params_B : list
            Lorentzian Parameters for Series B          (Default: [])
        time_or_phase : "time"/"phase"
            Define lags in Time or Phase                (Default: "phase")
        overall_lag : float (>0)
            Define overall lag                          (Default: 0)
        model_lag_array : ndarray
            Array for the lag parameters                (Default: np.zeros((6, 7)))
        write_data : 1/0
            Write the data to text files?               (Default: 1)
    """


    print("\n--------------- Starting CorrSim ---------------")

    ##                    Make Output Directory
    ##    --------------------------------------------------------

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    ##                    Set-up Calculations
    ##    --------------------------------------------------------

    ## Start of calculations
    start = time.time()

    mean_counts_A    = mean_counts_per_sec_A*time_res
    mean_counts_B    = mean_counts_per_sec_B*time_res

    time_arr    = np.arange(0, obs_length, time_res)    ## Create a set of time values
    num_bins    = len(time_arr)


    ## --- Define frequencies

    freq_min    = 1 / obs_length                        ## Minimum Frequency
    freq_max    = 1 / (2*time_res)                        ## Maximum Frequency
    freq_res    = 1 / obs_length                        ## Frequency Resolution


    frequencies    = np.arange(0, freq_max+(freq_res*0.1), freq_res)

    if len(frequencies) == len(time_arr)/2:
        frequencies = np.append(frequencies, frequencies[-1]+freq_res)

    ## In case of odd numbers of bins, make even instead.
    if len(time_arr) % 2 == 1:
        time_arr        = time_arr[0:-1]
        frequencies    = np.arange(0, freq_max, freq_res)

    ##        Power Spectrum Parameters
    ##    ----------------------------------

    print("Setting up model power spectra...")

    if power_model_type == "broken_powerlaw":

        power_model = corrfunc.power_model_broken_powerlaw(frequencies, ps_power, break_freq, coh_constant, coh_power, coh_break_freq)


    elif power_model_type == "lorentzians":

        power_model = corrfunc.power_model_lorentzians(frequencies, Lorentz_params_A, Lorentz_params_B)

        Lorentzians_A                = power_model[4]
        Lorentzians_B_coh            = power_model[5]
        Lorentzians_B_inc            = power_model[6]
        num_lorentz_A                = power_model[7]
        num_lorentz_B                = power_model[8]


    unnorm_power_spectrum_A        = power_model[0]
    unnorm_power_spectrum_B_coh    = power_model[1]
    unnorm_power_spectrum_B_inc    = power_model[2]
    coherence_model                = power_model[3]


    ##        Normalise Power Spectra
    ##    ----------------------------------

    print("Normalising power spectra...")

    F_rms_A_orig = F_rms_A
    F_rms_B_orig = F_rms_B

    ## Is Red Noise active? If so, be sure to adjust the normalisation so that the F_rms matches up:
    if apply_red == "A" or apply_red == "AB":
        if F_rms_A < F_rms_rednoise_A:
            raise Exception("F_rms_A must be greater than F_rms_rednoise_A.")
        F_rms_A = np.sqrt(F_rms_A**2 - F_rms_rednoise_A**2)
    if apply_red == "B" or apply_red == "AB":
        if F_rms_B < F_rms_rednoise_B:
            raise Exception("F_rms_B must be greater than F_rms_rednoise_B.")
        F_rms_B = np.sqrt(F_rms_B**2 - F_rms_rednoise_B**2)

    normalised_power_spectra    = corrfunc.normalise_power_spectra(unnorm_power_spectrum_A, F_rms_A, unnorm_power_spectrum_B_coh, \
        unnorm_power_spectrum_B_inc, F_rms_B, obs_length, time_res, num_bins, mean_counts_A, mean_counts_B)

    power_spectrum_A        = normalised_power_spectra[0]
    power_spectrum_B        = normalised_power_spectra[1]
    power_spectrum_B_coh    = normalised_power_spectra[2]
    power_spectrum_B_inc    = normalised_power_spectra[3]
    normalisation_A         = normalised_power_spectra[4]
    frac_rms_norm_A         = normalised_power_spectra[5]
    normalisation_B         = normalised_power_spectra[6]
    frac_rms_norm_B         = normalised_power_spectra[7]

    power_spectrum_model_A  = power_spectrum_A
    power_spectrum_model_B  = power_spectrum_B


    ##            Make Model Lags
    ##    ------------------------------

    model_lags = corrfunc.make_model_lags(frequencies, model_lag_array, overall_lag)

    if time_or_phase == 'time':
        model_time_lags        = model_lags
        model_phase_lags    = []
        for i in range(len(model_lags)):
            model_phase_lags.append(model_lags[i] * 2 * np.pi * frequencies[i])
        model_phase_lags = np.array(model_phase_lags)

    elif time_or_phase == 'phase':
        model_time_lags        = []
        model_phase_lags    = model_lags

        model_time_lags.append(np.nan)
        for i in range(1, len(model_lags)):
            model_time_lags.append(model_lags[i] / (2 * np.pi * frequencies[i]))
        model_time_lags = np.array(model_time_lags)



    ##          Calculate Model CCF
    ##  ----------------------------------

    if plot_model_ccf == 1:

        print("Calculating Model CCF...")

        model_ccf_full = corrfunc.calc_model_CCF(coherence_model, power_spectrum_model_A, frac_rms_norm_A, \
            power_spectrum_model_B, frac_rms_norm_B, model_lags, time_res, F_rms_A, F_rms_B)

        model_ccf_lags    = model_ccf_full[0]
        model_ccf         = model_ccf_full[1]



    ##          Convert Power Spectra to Lightcurve
    ##  ---------------------------------------------------

    print("Converting Power Spectra to Lightcurves:")

    ##  --- Use the Timmer-Koenig method to create Real and Imaginary parts of all three power spectra ---
    # print("Creating Re/Im parts of the Fourier Transformed Lightcurves (FTLs)...")
    FTL_Re_A,     FTL_Im_A      = corrfunc.create_ps_real_imag(frequencies, power_spectrum_A)
    FTL_Re_B_inc, FTL_Im_B_inc  = corrfunc.create_ps_real_imag(frequencies, power_spectrum_B_inc)

    ## --- If there are an odd number of frequencies, make the last imaginary bin equal to 0 ---
    if len(frequencies) % 2 == 1:
        FTL_Im_A[       len(frequencies)-1]    = 0
        FTL_Im_B_inc[   len(frequencies)-1]    = 0

    ##  --- Make complex numbers from FTL_B ---
    # print("Combining into complex numbers...")
    FTL_A_temp      = FTL_Re_A      + FTL_Im_A      *1j
    FTL_B_inc       = FTL_Re_B_inc  + FTL_Im_B_inc  *1j

    ##  --- Create a lag-shifted version of FTL_A ---
    print("Shifting the phase lags of Band A...")
    FTL_Norm        = np.append(0, np.sqrt(power_spectrum_B_coh[1:]/power_spectrum_A[1:]))
    FTL_B_coh       = FTL_Norm * FTL_A_temp*np.exp(-1j*model_phase_lags)
    FTL_A           = FTL_A_temp

    ##  --- Redefine the first bin to contain the mean count rate ---
    # print("Fixing Mean Count Rates...")
    FTL_A[0]        = mean_counts_A*num_bins+0j
    FTL_B_coh[0]    = mean_counts_B*num_bins+0j

    ##  --- At long last, combine Band B together! ---
    # print("Combining Band B FTLs...")
    FTL_B_total   = FTL_B_coh + FTL_B_inc

    ##  --- Reverse and Append all the Fourier-Transformed Lightcurves ---
    FTL_A       = corrfunc.complex_reverse_and_append(FTL_A)
    FTL_B_total = corrfunc.complex_reverse_and_append(FTL_B_total)

    # ##            Create Lightcurves
    # ##    ----------------------------------

    print("Creating Lightcurves...")
    flux_A  = corrfunc.ftl_to_lightcurve(FTL_A,         correct_negative=True)
    flux_B  = corrfunc.ftl_to_lightcurve(FTL_B_total,   correct_negative=True)

    ## --- Calculate F_rms ---
    print("Given F_rms_A = ", F_rms_A_orig)
    print("Given F_rms_B = ", F_rms_B_orig)
    print("[Pre-noise] F_rms_A = std(flux_A) / mean(flux_A) = {:0.3f}".format(np.std(flux_A)/np.mean(flux_A)))
    print("[Pre-noise] F_rms_B = std(flux_B) / mean(flux_B) = {:0.3f}".format(np.std(flux_B)/np.mean(flux_B)))


    ##         Add Noise
    ##  -----------------------

    ## --- Red Noise

    if apply_red == "A" or apply_red == "AB":
        print("Applying Red Noise for Band A...")

        # First, create real and imaginary parts of some simple red noise, and make them into a power spectrum.
        rednoise_1, rednoise_2  = corrfunc.create_red_noise_nums(power_spectrum_A)
        FTL_real, FTL_imag      = corrfunc.create_red_noise_real_imag(rednoise_1, rednoise_2, frequencies, mean_counts_A,  num_bins, 1, rednoise_slope_A)
        red_noise_A             = corrfunc.real_imag_to_ps(frequencies, FTL_real, FTL_imag)

        ## Next, we find the RMS of the power spectrum, and normalise it to what we actually want.
        red_noise_A_normalised, red_noise_A_norm, frac_rms_rd_norm_A = corrfunc.normalise_single_power_spectra(red_noise_A, F_rms_rednoise_A, obs_length, time_res, num_bins, mean_counts_A)

        ## Once we know this, we can go back and make some new red noise from scratch, this time with the correct normalisation!
        ## We do this so that the real and imaginary parts are correctly normalised, not just the power spectra.
        FTL_rednoise_real_A, FTL_rednoise_imag_A    = corrfunc.create_red_noise_real_imag(rednoise_1, rednoise_2, frequencies, mean_counts_A,  num_bins, (red_noise_A_norm/frac_rms_rd_norm_A)**0.5, rednoise_slope_A)
        red_noise_A_normalised_2                    = corrfunc.real_imag_to_ps(frequencies, FTL_rednoise_real_A, FTL_rednoise_imag_A)

        FTL_A = FTL_A + (FTL_rednoise_real_A + FTL_rednoise_imag_A * 1j)

        ## --- Re-make Lightcurves ---
        flux_A  = corrfunc.ftl_to_lightcurve(FTL_A,         correct_negative=True)

    if apply_red == "B" or apply_red == "AB":
        print("Applying Red Noise for Band B...")

        ## First, create real and imaginary parts of some simple red noise, and make them into a power spectrum.
        rednoise_1, rednoise_2  = corrfunc.create_red_noise_nums(power_spectrum_B)

        FTL_real, FTL_imag      = corrfunc.create_red_noise_real_imag(rednoise_1, rednoise_2, frequencies, mean_counts_B,  num_bins, 1, rednoise_slope_B)
        red_noise_B             = corrfunc.real_imag_to_ps(frequencies, FTL_real, FTL_imag)

        ## Next, we find the RMS of the power spectrum, and normalise it to what we actually want.
        red_noise_B_normalised, red_noise_B_norm, frac_rms_rd_norm_B = corrfunc.normalise_single_power_spectra(red_noise_B, F_rms_rednoise_B, obs_length, time_res, num_bins, mean_counts_B)

        ## Once we know this, we can go back and make some new red noise from scratch, this time with the correct normalisation!
        ## We do this so that the real and imaginary parts are correctly normalised, not just the power spectra.
        FTL_rednoise_real_B, FTL_rednoise_imag_B    = corrfunc.create_red_noise_real_imag(rednoise_1, rednoise_2, frequencies, mean_counts_B,  num_bins, (red_noise_B_norm/frac_rms_rd_norm_B)**0.5, rednoise_slope_B)
        red_noise_B_normalised_2                    = corrfunc.real_imag_to_ps(frequencies, FTL_rednoise_real_B, FTL_rednoise_imag_B)

        FTL_B_total = FTL_B_total + (FTL_rednoise_real_B + FTL_rednoise_imag_B * 1j)

        ## --- Re-make Lightcurves ---
        flux_B  = corrfunc.ftl_to_lightcurve(FTL_B_total,   correct_negative=True)

    ## --- Calculate F_rms ---
    if apply_red == "A" or apply_red == "AB":
        print("[Post-red noise] F_rms_A = std(flux_A) / mean(flux_A) = {:0.3f}".format(np.std(flux_A)/np.mean(flux_A)))
    if apply_red == "B" or apply_red == "AB":
        print("[Post-red noise] F_rms_B = std(flux_B) / mean(flux_B) = {:0.3f}".format(np.std(flux_B)/np.mean(flux_B)))


    ## --- Poisson Noise

    poisson_noise_multiplier_A = 0
    poisson_noise_multiplier_B = 0

    if apply_poisson == "A" or apply_poisson == "AB":
        print("Applying Poisson Noise to Band A...")
        flux_A = corrfunc.apply_poisson_noise(flux_A)
        poisson_noise_multiplier_A = 1

    if apply_poisson == "B" or apply_poisson == "AB":
        print("Applying Poisson Noise to Band B...")
        flux_B = corrfunc.apply_poisson_noise(flux_B)
        poisson_noise_multiplier_B = 1


    ## --- Scintillation Noise
    ## Scintillation noise is calculated through the GUI.
    ## Otherwise, calculate it manually, and define it in scin_noise_A/B
    scintillation_noise_multiplier_A = 0
    scintillation_noise_multiplier_B = 0

    # scin_noise_A = 0
    if apply_scintillation == "A" or apply_scintillation == "AB":
        print("Applying Scintillation Noise to Band A...")
        flux_A = corrfunc.apply_scintillation_noise(scin_noise_A, flux_A)
        scintillation_noise_multiplier_A = 1

    # scin_noise_B = 0
    if apply_scintillation == "B" or apply_scintillation == "AB":
        print("Applying Scintillation Noise to Band B...")
        flux_B = corrfunc.apply_scintillation_noise(scin_noise_B, flux_B)
        scintillation_noise_multiplier_B = 1


    ## --- Readout Noise
    readout_noise_multiplier_A = 0
    readout_noise_multiplier_B = 0

    if apply_readout == "A" or apply_readout == "AB":
        print("Applying Readout Noise to Band A...")
        flux_A            = corrfunc.apply_readout_noise(read_noise_A, flux_A)
        readout_noise_multiplier_A = 1

    if apply_readout == "B" or apply_readout == "AB":
        print("Applying Readout Noise to Band B...")
        flux_B            = corrfunc.apply_readout_noise(read_noise_B, flux_B)
        readout_noise_multiplier_B = 1


    ## --- Correct any negative counts ---
    # flux_A[flux_A<0] = 0
    # flux_B[flux_B<0] = 0


    ## --- What's the new F_rms? ---
    if apply_poisson == "A" or apply_poisson == "AB" or \
       apply_scintillation == "A" or apply_scintillation == "AB" or \
       apply_readout == "A" or apply_readout == "AB":
       print("[Post-all noise] F_rms_A = std(flux_A) / mean(flux_A) = {:0.3f}".format(np.std(flux_A)/np.mean(flux_A)))
    if apply_poisson == "B" or apply_poisson == "AB" or \
       apply_scintillation == "B" or apply_scintillation == "AB" or \
       apply_readout == "B" or apply_readout == "AB":
       print("[Post-all noise] F_rms_B = std(flux_B) / mean(flux_B) = {:0.3f}".format(np.std(flux_B)/np.mean(flux_B)))


    ## --- Read/Scintillation Factors ---
    poisson_noise_factor_A = 2/mean_counts_per_sec_A
    poisson_noise_factor_B = 2/mean_counts_per_sec_B
    # read_noise_factor_A = (read_noise_A/np.mean(flux_A))**2/500
    # read_noise_factor_B = (read_noise_B/np.mean(flux_B))**2/500
    # scin_noise_factor_A = ((scin_noise_A*(time_res)) / (np.mean(flux_A)/(mean_counts_per_sec_A*time_res))**2 )*2
    # scin_noise_factor_B = ((scin_noise_B*(time_res)) / (np.mean(flux_B)/(mean_counts_per_sec_B*time_res))**2 )*2
    read_noise_factor_A = ((read_noise_A**2 * obs_length * 2) / (num_bins * np.mean(flux_A)**2))
    read_noise_factor_B = ((read_noise_B**2 * obs_length * 2) / (num_bins * np.mean(flux_B)**2))
    scin_noise_factor_A = (((np.sqrt(scin_noise_A)*np.mean(flux_A))**2 * obs_length * 2) / (num_bins * np.mean(flux_A)**2))
    scin_noise_factor_B = (((np.sqrt(scin_noise_B)*np.mean(flux_B))**2 * obs_length * 2) / (num_bins * np.mean(flux_B)**2))

    ## Optional: Plot the red noise, to make sure it's okay.
    if plot_rednoise == 1:
        print("Plotting Red Noise...")

        plt.figure(figsize=(10,5))

        if apply_red == "A" or apply_red == "AB":
            power_spectrum_A_w_rn                    = corrfunc.real_imag_to_ps(frequencies, FTL_A.real, FTL_A.imag)
            plt.step(frequencies, red_noise_A_normalised_2*frac_rms_rd_norm_A, where="mid", label="A (Red Noise)", color="darkturquoise")
            plt.plot(frequencies,  power_spectrum_A*frac_rms_norm_A, linestyle="-", label="A (Model)", zorder=15, color="blue")
            plt.step(frequencies,  power_spectrum_A_w_rn*frac_rms_norm_A, where="mid", label="A (Model + Red Noise)", color="cornflowerblue")

        if apply_red == "B" or apply_red == "AB":
            power_spectrum_B_w_rn                    = corrfunc.real_imag_to_ps(frequencies, FTL_B_total.real, FTL_B_total.imag)
            plt.step(frequencies, red_noise_B_normalised_2*frac_rms_rd_norm_B, where="mid", label="B (Red Noise)", color="sandybrown")
            plt.plot(frequencies,  power_spectrum_B*frac_rms_norm_B, linestyle="-", label="B (Model)", zorder=15, color="red")
            plt.step(frequencies,  power_spectrum_B_w_rn*frac_rms_norm_B, where="mid", label="B (Model + Red Noise)", color="lightcoral")

        if apply_red == "AB":
            handles, labels = plt.gca().get_legend_handles_labels()
            order = [0,2,4,1,3,5]
            plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
        else:
            plt.legend()

        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Frac. RMS Power')
        plt.grid(linestyle=":", zorder=5)

        fileOut = "./" + output_dir + "/" + fileprefix + "Power Spectrum Red Noise.png"
        plt.savefig(fileOut, bbox_inches='tight')

        plt.close('all')

    ## --- End of calculations ---
    print("")
    print("=======")
    print("With " + str(num_bins) + " bins...", sep="")
    print(f'Time: {time.time() - start}s')
    print("=======")
    print("")


    if write_data == 1:

        ##            Write the Data
        ##    ----------------------------------


        print("Writing Data...")

        ## Model Fourier Data
        if plot_model_fourier == 1:
            fileOut = "./" + output_dir + "/" + fileprefix + "Model Fourier Data.txt"
            output = pd.DataFrame({    'Frequencies' : frequencies, 'Band_A_Power' : power_spectrum_model_A*frac_rms_norm_A,
                                    'Band_B_Power' : power_spectrum_model_B*frac_rms_norm_B, 'Coherence' : coherence_model,
                                    'Phase Lags' : model_phase_lags,
                                    'Time Lags' : model_time_lags})
            output.to_csv(fileOut, index=False)

        ## Model CCF
        if plot_model_ccf == 1:
            fileOut = "./" + output_dir + "/" + fileprefix + "Model CCF Data.txt"
            output = pd.DataFrame({'model_ccf_lags' : model_ccf_lags, 'model_ccf' : model_ccf})
            output.to_csv(fileOut, index=False)


        ## Power Spectra
        fileOut = "./" + output_dir + "/" + fileprefix + "Power Spectrum Data.txt"
        output = pd.DataFrame({'Frequencies' : frequencies, 'Band_A_Power' : power_spectrum_A, 'Band_B_Power' : power_spectrum_B})
        output.to_csv(fileOut, index=False)

        ## Lightcurves
        fileOut = "./" + output_dir + "/" + fileprefix + "Lightcurve Data.txt"
        output = pd.DataFrame({'A_Time' : time_arr, 'Band_A_Flux' : flux_A, 'Band_B_Flux' : flux_B})
        output.to_csv(fileOut, index=False)


        ##            Write a short log file of the Details
        ##    -----------------------------------------------------

        fileOut = output_dir + "/" + fileprefix + "Log File.txt"

        with open(fileOut, 'w') as the_file:
            the_file.write("File Prefix:\t\t" + fileprefix + '\n')
            the_file.write("Observation Length:\t" + str(obs_length) + "s\n")
            the_file.write("Time Resolution:\t" + str(time_res) + "s\n")
            the_file.write("Number of Bins:\t\t" + str(num_bins) + "s\n")
            the_file.write("Mean Counts/s (A):\t" + str(mean_counts_per_sec_A) + "\n")
            the_file.write("Mean Counts/s (B):\t" + str(mean_counts_per_sec_B) + "\n")
            the_file.write("F_rms_A:\t\t" + str(F_rms_A_orig) + "\n")
            the_file.write("F_rms_B:\t\t" + str(F_rms_B_orig) + "\n")
            the_file.write("Apply Red Noise:\t" + apply_red + "\n")
            the_file.write("Apply Poisson Noise:\t" + apply_poisson + "\n")
            the_file.write("Apply Readout Noise:\t" + apply_readout + "\n")
            the_file.write("Readout Noise (A):\t" + str(read_noise_A) + "\n")
            the_file.write("Readout Noise (B):\t" + str(read_noise_B) + "\n")
            the_file.write("Apply Scin. Noise:\t" + apply_scintillation + "\n")
            the_file.write("\n")
            the_file.write("Min Freq:\t\t" + str(freq_min) + "\n")
            the_file.write("Max Freq:\t\t" + str(freq_max) + "\n")
            the_file.write("Freq Res:\t\t" + str(freq_res) + "\n")
            the_file.write("\n")
            the_file.write("Frac_RMS_Norm_1:\t" + str(frac_rms_norm_A) + "\n")
            the_file.write("Frac_RMS_Norm_2:\t" + str(frac_rms_norm_B) + "\n")
            the_file.write("Returned F_rms_1:\t" + str(np.std(flux_A)/np.mean(flux_A)) + "\n")
            the_file.write("Returned F_rms_2:\t" + str(np.std(flux_B)/np.mean(flux_B)) + "\n")
            the_file.write("\n")
            the_file.write("Max Flux_1:\t\t" + str(max(flux_A)) + "\n")
            the_file.write("Min Flux_1:\t\t" + str(min(flux_A)) + "\n")
            the_file.write("Mean Flux_1:\t\t" + str(np.mean(flux_A)) + "\n")
            the_file.write("\n")
            the_file.write("Max Flux_2:\t\t" + str(max(flux_B)) + "\n")
            the_file.write("Min Flux_2:\t\t" + str(min(flux_B)) + "\n")
            the_file.write("Mean Flux_2:\t\t" + str(np.mean(flux_B)) + "\n")
            the_file.write("\n")
            the_file.write("Power Spectra:\t\t" + power_model_type + "\n")
            the_file.write("- Broken Powerlaw values:\n")
            the_file.write("Power Spectral Index\t:" + str(ps_power) + "\n")
            the_file.write("Pow. Spec. Break Freq.:\t" + str(break_freq) + "\n")
            the_file.write("Coherence Constant:\t" + str(coh_constant) + "\n")
            the_file.write("Coherence Power:\t" + str(coh_power) + "\n")
            the_file.write("Coherence Break Freq.:\t" + str(coh_break_freq) + "\n")
            the_file.write("- Lorentzian values:\n")
            the_file.write("Series A:\t\t" + str(Lorentz_params_A) + "\n")
            the_file.write("Series B:\t\t" + str(Lorentz_params_B) + "\n")
            the_file.write("\n")
            the_file.write("Lag Method:\t\t" + time_or_phase + "\n")
            the_file.write("Overall Lag:\t\t" + str(overall_lag) + "\n")
            the_file.write("Model Lag Array 1:\t" + str(model_lag_array[0,:]) + "\n")
            the_file.write("Model Lag Array 2:\t" + str(model_lag_array[1,:]) + "\n")
            the_file.write("Model Lag Array 3:\t" + str(model_lag_array[2,:]) + "\n")
            the_file.write("Model Lag Array 4:\t" + str(model_lag_array[3,:]) + "\n")
            the_file.write("Model Lag Array 5:\t" + str(model_lag_array[4,:]) + "\n")
            the_file.write("Model Lag Array 6:\t" + str(model_lag_array[5,:]) + "\n")


        ###
        ### ===== Plot Lightcurves =====
        ###

        insetA_1 = 0
        insetA_2 = 9.5

        insetB_1 = 2
        insetB_2 = 3

        fig = plt.figure(figsize=(10,5))

        plot_margin = 0.1




        ## --- Main (Band A) ---
        left = 0.0
        bottom = 0.755
        width = 1.0
        height = 0.245
        main_A_ax = fig.add_axes([left, bottom, width, height], [])

        main_A_ax.plot(time_arr, flux_A, linestyle='-', lw=1, color="blue")

        main_A_ax.set_xlabel('Time (s)')
        main_A_ax.set_ylabel('Counts/Bin')

        main_A_ax.xaxis.tick_top()
        main_A_ax.xaxis.set_label_position("top")

        main_A_ax.grid(linestyle=":", zorder=5)

        main_A_ax.axvline(insetA_1, linestyle="--", linewidth=1, color="blue", zorder=9)
        main_A_ax.axvline(insetA_2, linestyle="--", linewidth=1, color="blue", zorder=9)

        ## --- Main (Band B) ---
        left = 0.0
        bottom = 0.51
        width = 1.0
        height = 0.245
        main_B_ax = fig.add_axes([left, bottom, width, height], [])

        main_B_ax.plot(time_arr, flux_B, linestyle='-', lw=1, color="red")

        main_B_ax.set_ylabel('Counts/Bin')
        main_B_ax.xaxis.set_ticklabels([])

        main_B_ax.yaxis.tick_right()
        main_B_ax.yaxis.set_label_position("right")

        main_B_ax.grid(linestyle=":", zorder=5)

        main_B_ax.axvline(insetA_1, linestyle="--", linewidth=1, color="blue", zorder=9)
        main_B_ax.axvline(insetA_2, linestyle="--", linewidth=1, color="blue", zorder=9)


        ## --- Zoom 1 (Band A) ---
        left = 0.0
        bottom = 0.245
        width = 0.69
        height = 0.245
        zoomA_A_ax = fig.add_axes([left, bottom, width, height], [])

        zoomA_A_ax.plot(time_arr, flux_A, linestyle='-', lw=1, color="blue")

        zoomA_A_ax.set_xlim([insetA_1, insetA_2])
        zoomA_A_ax.xaxis.tick_top()
        zoomA_A_ax.xaxis.set_ticklabels([])

        zoomA_A_ax.set_ylabel('Counts/Bin')

        zoomA_A_ax.grid(linestyle=":", zorder=5)

        zoomA_A_ax.axvline(insetB_1, linestyle="--", linewidth=1, color="blue", zorder=9)
        zoomA_A_ax.axvline(insetB_2, linestyle="--", linewidth=1, color="blue", zorder=9)

        ## --- Zoom 1 (Band B) ---
        left = 0.0
        bottom = 0.0
        width = 0.69
        height = 0.245
        zoomA_B_ax = fig.add_axes([left, bottom, width, height], [])

        zoomA_B_ax.plot(time_arr, flux_B, linestyle='-', lw=1, color="red")

        zoomA_B_ax.set_xlim([insetA_1, insetA_2])

        zoomA_B_ax.set_xlabel('Time (s)')
        zoomA_B_ax.set_ylabel('Counts/Bin')

        zoomA_B_ax.grid(linestyle=":", zorder=5)

        zoomA_B_ax.axvline(insetB_1, linestyle="--", linewidth=1, color="blue", zorder=9)
        zoomA_B_ax.axvline(insetB_2, linestyle="--", linewidth=1, color="blue", zorder=9)


        ## --- Zoom 2 (Band A) ---
        left = 0.71
        bottom = 0.245
        width = 0.29
        height = 0.245
        zoomB_A_ax = fig.add_axes([left, bottom, width, height], [])

        zoomB_A_ax.plot(time_arr, flux_A, linestyle='-', lw=1, color="blue")

        zoomB_A_ax.set_xlim([insetB_1, insetB_2])
        zoomB_A_ax.xaxis.tick_top()
        zoomB_A_ax.xaxis.set_ticklabels([])

        zoomB_A_ax.yaxis.tick_right()
        zoomB_A_ax.yaxis.set_label_position("right")
        zoomB_A_ax.set_ylabel('Counts/Bin')

        zoomB_A_ax.grid(linestyle=":", zorder=5)


        ## --- Zoom 2 (Band B) ---
        left = 0.71
        bottom = 0.0
        width = 0.29
        height = 0.245
        zoomB_B_ax = fig.add_axes([left, bottom, width, height], [])

        zoomB_B_ax.plot(time_arr, flux_B, linestyle='-', lw=1, color="red")

        zoomB_B_ax.set_xlim([insetB_1, insetB_2])

        zoomB_B_ax.yaxis.tick_right()
        zoomB_B_ax.yaxis.set_label_position("right")
        zoomB_B_ax.set_xlabel('Time (s)')
        zoomB_B_ax.set_ylabel('Counts/Bin')

        zoomB_B_ax.grid(linestyle=":", zorder=5)



        fileOut = "./" + output_dir + "/" + fileprefix + "Lightcurves.png"
        plt.savefig(fileOut, bbox_inches='tight')

        plt.close('all')

    ###                Plotting
    ###    ----------------------------------



    ###
    ### ===== Plot Power Spectra Components =====
    ###


    ###
    ### ===== Plot Model CCF =====
    ###

    if plot_model_ccf == 1:

        plt.figure(figsize=(8,5))

        plt.plot(model_ccf_lags, model_ccf, color="k")

        plt.xlabel('Lag (s)')
        plt.ylabel('CCF Coefficient')

        plt.xlim(-10,10)

        plt.grid(linestyle=":", zorder=5)

        plt.axhline(0, linestyle="--", linewidth=1, zorder=9)
        plt.axvline(0, linestyle="--", linewidth=1, zorder=9)

        fig = plt.gcf()

        savename = "./" + output_dir + "/" + fileprefix + "Model CCF (10s lag).png"
        fig.savefig(savename, bbox_inches='tight')

        plt.xlim(-0.5,0.5)

        savename = "./" + output_dir + "/" + fileprefix + "Model CCF (2s lag).png"
        plt.savefig(savename, bbox_inches='tight')


    ###
    ### ===== Calculate CCF =====
    ###

    bin_kwarg = ""
    if ccf_binning > 1:
        bin_kwarg = " (Binned " + str(ccf_binning) + ")"

    if calculate_ccf == 1:

        print("Calculating CCF...")

        ccf_av = corrfunc.simulated_ccf(flux_A, flux_B, time_res, ccf_segment, ccf_binning)

        lower = int(len(ccf_av)/3)
        upper = int(2*len(ccf_av)/3)
        mean_ccf_error = np.mean(ccf_av[lower:upper,2])
        stdv_ccf_error = np.std(ccf_av[lower:upper,2])

        print(" ")
        print("------------------------------")
        print("Mean CCF Error: ", '{:g}'.format(float('{:.{p}g}'.format(mean_ccf_error, p=3))), " +/- ", '{:g}'.format(float('{:.{p}g}'.format(stdv_ccf_error, p=3))))
        print("------------------------------")
        print(" ")

        ## Write the data
        fileOut = "./" + output_dir + "/" + fileprefix + "CCF Data (" + str(ccf_segment) + "s" + bin_kwarg + ").txt"# + bin_kwarg
        output = pd.DataFrame({'lag' : ccf_av[:,0], 'ccf' : ccf_av[:,1], 'ccf_err' : ccf_av[:,2]})
        output.to_csv(fileOut, index=False)


        ###
        ### ===== Plot Simulated CCF =====
        ###

        plt.figure(figsize=(10,5))

        plt.step(ccf_av[:,0], ccf_av[:,1], where='mid', lw=1, color="black")

        plt.xlim(-ccf_segment/3, ccf_segment/3)

        plt.xlabel('Lag (s)')
        plt.ylabel('CCF Coefficient')

        plt.grid(linestyle=":", zorder=5)

        fileOut = "./" + output_dir + "/" + fileprefix + "Averaged CCF (" + str(ccf_segment) + "s" + bin_kwarg + ").png"
        plt.savefig(fileOut, bbox_inches='tight')

        plt.close('all')


        ###
        ### ===== Plot Simulated CCF with Model (If requested) =====
        ###

        if plot_model_ccf == 1:
            plt.figure(figsize=(10,5))

            plt.step(ccf_av[:,0], ccf_av[:,1], where='mid', lw=1, color="red", label="Simulated CCF")

            plt.plot(model_ccf_lags, model_ccf, color="k", label="Model CCF")

            plt.xlim(-ccf_segment/3, ccf_segment/3)

            plt.xlabel('Lag (s)')
            plt.ylabel('CCF Coefficient')

            plt.grid(linestyle=":", zorder=5)

            plt.legend()

            fileOut = "./" + output_dir + "/" + fileprefix + "Averaged CCF + Model (" + str(ccf_segment) + "s" + bin_kwarg + ").png"

            plt.savefig(fileOut, bbox_inches='tight')

            plt.close('all')


    ###
    ### ===== Calculate Fourier =====
    ###

    if calculate_fourier == 1:

        print("Beginning Fourier Calculations...")

        segment_size = time_res * fourier_bins

        print(str(int(np.floor(num_bins/fourier_bins))) + " segments to be used, " + str(fourier_bins) + " bins (" + str(segment_size) + " seconds) each.")

        fourier_data = corrfunc.fourier(time_arr, flux_A, flux_B, time_res, \
            segment_size, fourier_rebin, reference_freq, apply_poisson, remove_whitenoise)

        fourier_output = np.transpose(fourier_data)

        fileOut = "./" + output_dir + "/" + fileprefix + "Fourier Data (bins " + str(fourier_bins) + " rebin " + str(fourier_rebin) + ").txt"
        output = pd.DataFrame(fourier_output)
        output_header = ["Frequency", "Band_A_Power", "Band_A_Error", "Band_B_Power", "Band_B_Error", \
            "Coherence", "Coherence_Error", "Phase_Lags", "Phase_Lags_Error", "Time_Lags", "Time_Lags_Inv", "Time_Lags_Error"]
        output.to_csv(fileOut, index=False, header=output_header)

        fourier_freqs            = fourier_data[0]
        powerspectra_A            = fourier_data[1]
        powerspectra_A_err        = fourier_data[2]
        powerspectra_B            = fourier_data[3]
        powerspectra_B_err        = fourier_data[4]
        coherence                = fourier_data[5]
        coherence_err            = fourier_data[6]
        phase_lags                = fourier_data[7]
        phase_lags_err            = fourier_data[8]
        time_lags                = fourier_data[9]
        time_lags_inv            = fourier_data[10]
        time_lags_err            = fourier_data[11]

        print(" ")
        print("-------------------------------")

        print("Median Absolute Fourier Errors (as percentage of measured values):")

        print("Powerspectrum A:\t",    round(np.nanmedian(abs((powerspectra_A_err/powerspectra_A)*100)),                 2),    "%", sep="")
        print("Powerspectrum B:\t",    round(np.nanmedian(abs((powerspectra_B_err/powerspectra_B)*100)),                 2),    "%", sep="")
        print("Coherence:\t\t",        round(np.nanmedian(abs((coherence_err/coherence)*100)),                         2),    "%", sep="")
        print("Phase Lags:\t\t",    round(np.nanmedian(abs((np.array(phase_lags_err)/np.array(phase_lags))*100)),    2),    "%", sep="")
        print("Time Lags:\t\t",        round(np.nanmedian(abs((np.array(time_lags_err) /np.array(time_lags)) *100)),    2),    "%", sep="")

        print("-------------------------------")
        print(" ")
        print("Plotting Fourier...")
        print(" ")

        fig, axs = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(5,8), gridspec_kw={'height_ratios': [10,7,7,10]})

        ## --- Power Spectra ---

        ax = axs[0]

        ax.errorbar(fourier_freqs, fourier_freqs * powerspectra_A, yerr=fourier_freqs * powerspectra_A_err, \
            ms=4, lw=0.5, capsize=3, fmt='bo', zorder=9, label="Band A")
        ax.errorbar(fourier_freqs, fourier_freqs * powerspectra_B, yerr=fourier_freqs * powerspectra_B_err, \
            ms=4, lw=0.5, capsize=3, fmt='ro', zorder=9, label="Band B")

        ax.set_ylabel("Freq * RMS Squared")
        ax.set_yscale('log', nonpositive='clip')
        ax.grid(linestyle='--')

        # ax.text(reference_freq, max(fourier_freqs * powerspectra_A), "Band A", horizontalalignment='center', \
        #     verticalalignment='center', style='italic', color='blue')
        # ax.text(reference_freq, max(fourier_freqs * powerspectra_B), "Band B", horizontalalignment='center', \
        #     verticalalignment='center', style='italic', color='red')

        ax.legend()

        ## --- Coherence ---
        ax = axs[1]
        ax.set_xscale('log')
        ax.set_yscale('log', nonpositive="clip")
        ax.set_ylim([min(coherence)/2, 2])

        ax.errorbar(fourier_freqs, coherence, yerr=coherence_err, lw=0.5, fmt='mo', zorder=9)
        ax.set_ylabel("Coherence")
        ax.grid(linestyle='--')


        ## --- Phase Lags ---
        ax = axs[2]
        ax.set_xscale('log')

        ax.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(corrfunc.format_func))
        ax.set_ylim([-3.8, 3.8])

        ax.errorbar(fourier_freqs, phase_lags, yerr=phase_lags_err, lw=0.5, fmt='mo')

        ax.set_ylabel("Phase Lag (Radians)")
        ax.grid(linestyle='--')


        ## --- Time Lags ---
        ax = axs[3]
        ax.set_xscale('log')
        ax.set_yscale('log', nonpositive="clip")


        ax.errorbar(fourier_freqs, time_lags,     yerr=time_lags_err, lw=0.5, fmt='mo', label = "+ve Lags")
        ax.errorbar(fourier_freqs, time_lags_inv, yerr=time_lags_err, lw=0.5, fmt='mo', fillstyle = "none", label = "-ve Lags")

        ax.legend()

        ax.set_ylabel("Time Lag (s)")
        ax.set_xlabel("Frequency (Hz)")
        ax.grid(linestyle='--')

        ## ----- PLOT END -----

        plt.subplots_adjust(bottom=0.1, left=0.15, top=0.98, right=0.98, hspace=0)

        savename = "./" + output_dir + "/" + fileprefix + "Fourier.pdf"

        plt.savefig(savename, bbox_inches='tight')

        plt.close('all')

    ## Find first and last for easy Red Noise plotting
    freq_first_last = np.array([frequencies[1], frequencies[-1]])


    if plot_model_fourier == 1:


        fig, axs = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(5,8), gridspec_kw={'height_ratios': [10,7,7,10]})

        ## --- Power Spectra ---


        model_frac_rms_norm_A = ((2 * time_res) / (num_bins * np.mean(flux_A)**2))                ## Find F_rms normalisation (Vaughan 2003, Pg. 12)
        model_frac_rms_norm_B = ((2 * time_res) / (num_bins * np.mean(flux_B)**2))                ## Find F_rms normalisation (Vaughan 2003, Pg. 12)


        ax = axs[0]

        ax.plot(frequencies,    frequencies*
            (power_spectrum_model_A*frac_rms_norm_A
            +poisson_noise_multiplier_A*poisson_noise_factor_A
            +readout_noise_multiplier_A*read_noise_factor_A
            +scintillation_noise_multiplier_A*scin_noise_factor_A),
            linestyle="-", zorder =14, color="cornflowerblue", label="Series A")
        ax.plot(frequencies,    frequencies*
            (power_spectrum_model_B*frac_rms_norm_B
            +poisson_noise_multiplier_B*poisson_noise_factor_B
            +readout_noise_multiplier_B*read_noise_factor_B
            +scintillation_noise_multiplier_B*scin_noise_factor_B),
            linestyle="-", zorder =14, color="salmon", label="Series B")

        ## Red Noise
        if apply_red == "A" or apply_red == "AB":
            ax.plot(freq_first_last, freq_first_last*(frac_rms_rd_norm_A*freq_first_last)**rednoise_slope_A, linestyle=":", zorder = 15, color="black", label="Red noise A")
            ax.plot(freq_first_last, freq_first_last*red_noise_A_norm*freq_first_last**rednoise_slope_A, linestyle="--", zorder = 13, color="cornflowerblue", label="Red noise A")
        if apply_red == "B" or apply_red == "AB":
            ax.plot(freq_first_last, freq_first_last*(frac_rms_rd_norm_B*freq_first_last)**rednoise_slope_B, linestyle=":", zorder = 13, color="salmon", label="Red noise B")

        ax.legend()

        ax.set_ylabel("Freq * RMS Squared")
        ax.set_yscale('log', nonpositive='clip')
        ax.grid(linestyle='--')

        ## --- Coherence ---
        ax = axs[1]
        ax.set_xscale('log')
        ax.set_yscale('log', nonpositive="clip")
        ax.set_ylim([min(coherence_model[1:-1])/2, 2])

        ax.plot(frequencies, coherence_model, linestyle="-", zorder =14, color="black")

        ax.set_ylabel("Coherence")
        ax.grid(linestyle='--')


        ## --- Phase Lags ---
        ax = axs[2]
        ax.set_xscale('log')

        ax.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(corrfunc.format_func))
        ax.set_ylim([-3.8, 3.8])

        ax.plot(frequencies, model_phase_lags-(4*np.pi), lw=1, zorder=8, linestyle='-', color="k")
        ax.plot(frequencies, model_phase_lags-(2*np.pi), lw=1, zorder=8, linestyle='-', color="k")
        ax.plot(frequencies, model_phase_lags, lw=1, zorder=8, linestyle='-', color="k")
        ax.plot(frequencies, model_phase_lags+(2*np.pi), lw=1, zorder=8, linestyle='-', color="k")
        ax.plot(frequencies, model_phase_lags+(4*np.pi), lw=1, zorder=8, linestyle='-', color="k")

        ax.set_ylabel("Phase Lag (Radians)")
        ax.grid(linestyle='--')


        ## --- Time Lags ---
        ax = axs[3]
        ax.set_xscale('log')
        ax.set_yscale('log', nonpositive="clip")

        ax.plot(frequencies,  np.array(model_time_lags), lw=1, zorder=8, linestyle='-', color="k")
        ax.plot(frequencies, -np.array(model_time_lags), lw=1, zorder=8, linestyle='--', color="k")

        ax.set_ylabel("Time Lag (s)")
        ax.set_xlabel("Frequency (Hz)")
        ax.grid(linestyle='--')

        ## ----- PLOT END -----

        plt.subplots_adjust(bottom=0.1, left=0.15, top=0.98, right=0.98, hspace=0)

        savename = "./" + output_dir + "/" + fileprefix + "Fourier_Models.pdf"

        plt.savefig(savename, bbox_inches='tight')

        plt.close('all')


    if plot_model_fourier == 1 and calculate_fourier == 1:

        fig, axs = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(5,8), gridspec_kw={'height_ratios': [10,7,7,10]})

        ## --- Power Spectra ---

        ax = axs[0]

        ax.errorbar(fourier_freqs, fourier_freqs * powerspectra_A, yerr=fourier_freqs * powerspectra_A_err, \
            ms=4, lw=0.5, capsize=3, fmt='bo', zorder=9, label="Series A")
        ax.errorbar(fourier_freqs, fourier_freqs * powerspectra_B, yerr=fourier_freqs * powerspectra_B_err, \
            ms=4, lw=0.5, capsize=3, fmt='ro', zorder=9, label="Series B")

        ax.plot(frequencies,    frequencies*
            (power_spectrum_model_A*frac_rms_norm_A
            +poisson_noise_multiplier_A*poisson_noise_factor_A
            +readout_noise_multiplier_A*read_noise_factor_A
            +scintillation_noise_multiplier_A*scin_noise_factor_A),
            linestyle="-", zorder =14, color="cornflowerblue")
        ax.plot(frequencies,    frequencies*
            (power_spectrum_model_B*frac_rms_norm_B
            +poisson_noise_multiplier_B*poisson_noise_factor_B
            +readout_noise_multiplier_B*read_noise_factor_B
            +scintillation_noise_multiplier_B*scin_noise_factor_B),
            linestyle="-", zorder =14, color="salmon")

        ## --- Red Noise ----


        if apply_red == "A" or apply_red == "AB":
            ## Create Red Noise Model
            unnorm_red_noise_model_A    = frequencies[1:]**(rednoise_slope_A)

            ## Normalise
            red_noise_model_A, red_noise_model_norm_A, red_noise_model_f_rms_norm_A = corrfunc.normalise_single_power_spectra( \
                unnorm_red_noise_model_A, F_rms_rednoise_A, obs_length, time_res, num_bins, mean_counts_A)

            ## Plot Model
            ax.plot(frequencies[1:], frequencies[1:]*red_noise_model_A*red_noise_model_f_rms_norm_A, linestyle="--", zorder = 13, color="cornflowerblue")#, label="Red Noise (A)")

        if apply_red == "B" or apply_red == "AB":
            ## Create Red Noise Model
            unnorm_red_noise_model_B    = frequencies[1:]**(rednoise_slope_B)

            ## Normalise
            red_noise_model_B, red_noise_model_norm_B, red_noise_model_f_rms_norm_B = corrfunc.normalise_single_power_spectra( \
                unnorm_red_noise_model_B, F_rms_rednoise_B, obs_length, time_res, num_bins, mean_counts_B)

            ## Plot Model
            ax.plot(frequencies[1:], frequencies[1:]*red_noise_model_B*red_noise_model_f_rms_norm_B, linestyle="--", zorder = 13, color="salmon")#, label="Red Noise (B)")


        ax.set_ylabel("Freq * RMS Squared")
        ax.set_yscale('log', nonpositive='clip')
        ax.grid(linestyle='--')

        ax.legend()

        ## --- Coherence ---
        ax = axs[1]
        ax.set_xscale('log')
        ax.set_yscale('log', nonpositive="clip")
        ax.set_ylim([min(coherence)/2, 2])

        ax.errorbar(fourier_freqs, coherence, yerr=coherence_err, lw=0.5, fmt='mo', zorder=9)
        ax.plot(frequencies, coherence_model, linestyle="-", zorder =14, color="black")

        ax.set_ylabel("Coherence")
        ax.grid(linestyle='--')


        ## --- Phase Lags ---
        ax = axs[2]
        ax.set_xscale('log')

        ax.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(corrfunc.format_func))
        ax.set_ylim([-3.8, 3.8])

        ax.errorbar(fourier_freqs, phase_lags, yerr=phase_lags_err, lw=0.5, fmt='mo')
        ax.plot(frequencies, model_phase_lags-(4*np.pi), lw=1, zorder=8, linestyle='-', color="k")
        ax.plot(frequencies, model_phase_lags-(2*np.pi), lw=1, zorder=8, linestyle='-', color="k")
        ax.plot(frequencies, model_phase_lags, lw=1, zorder=8, linestyle='-', color="k")
        ax.plot(frequencies, model_phase_lags+(2*np.pi), lw=1, zorder=8, linestyle='-', color="k")
        ax.plot(frequencies, model_phase_lags+(4*np.pi), lw=1, zorder=8, linestyle='-', color="k")

        ax.set_ylabel("Phase Lag (Radians)")
        ax.grid(linestyle='--')


        ## --- Time Lags ---
        ax = axs[3]
        ax.set_xscale('log')
        ax.set_yscale('log', nonpositive="clip")

        ax.errorbar(fourier_freqs, time_lags,     yerr=time_lags_err, lw=0.5, fmt='mo', label = "+ve Lags")
        ax.errorbar(fourier_freqs, time_lags_inv, yerr=time_lags_err, lw=0.5, fmt='mo', fillstyle="none", label = "-ve Lags")

        ax.plot(frequencies,  np.array(model_time_lags), lw=1, zorder=8, linestyle='-', color="k")
        ax.plot(frequencies, -np.array(model_time_lags), lw=1, zorder=8, linestyle='--', color="k")

        ax.legend()

        ax.set_ylabel("Time Lag (s)")
        ax.set_xlabel("Frequency (Hz)")
        ax.grid(linestyle='--')

        ## ----- PLOT END -----

        plt.subplots_adjust(bottom=0.1, left=0.15, top=0.98, right=0.98, hspace=0)

        savename = "./" + output_dir + "/" + fileprefix + "Fourier_Models_&_Data.pdf"

        plt.savefig(savename, bbox_inches='tight')

        plt.close('all')




    ## Compilation plot
    if plot_model_ccf == 1 and calculate_ccf == 1 and plot_model_fourier == 1 and calculate_fourier == 1:

        fig = plt.figure(figsize=(6, 13))
        plot_margin = 0.1

        ## ----- Lightcurves -----
        insetA_1 = 0
        insetA_2 = np.max(time_arr)

        insetB_1 = np.max(time_arr)/10
        insetB_2 = 2*np.max(time_arr)/10

        insetC_1 = np.max(time_arr)/10 + np.max(time_arr)/200
        insetC_2 = np.max(time_arr)/10 + 2*np.max(time_arr)/200

        ## --- Inset 1 ---
        ## Band A
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.0, 0.94,  1.0, 0.06], [])

        ax.plot(time_arr, flux_A, linestyle='-', lw=1, color="blue")
        ax.set_xlim([insetA_1, insetA_2])

        ax.minorticks_on()
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Counts')

        ax.grid(linestyle=":", zorder=5)

        ax.axvline(insetB_1, linestyle="--", linewidth=1, color="blue", zorder=9)
        ax.axvline(insetB_2, linestyle="--", linewidth=1, color="blue", zorder=9)

        ## Band B
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.0, 0.88,  1.0, 0.06], [])

        ax.plot(time_arr, flux_B, linestyle='-', lw=1, color="red")
        ax.set_xlim([insetA_1, insetA_2])

        ax.minorticks_on()
        ax.xaxis.set_ticklabels([])
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_ylabel('Counts')

        ax.grid(linestyle=":", zorder=5)

        ax.axvline(insetB_1, linestyle="--", linewidth=1, color="blue", zorder=9)
        ax.axvline(insetB_2, linestyle="--", linewidth=1, color="blue", zorder=9)


        ## --- Inset 2 ---
        ## Band A
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.0, 0.80, 0.49, 0.06], [])
        ax.plot(time_arr, flux_A, linestyle='-', lw=1, color="blue")
        ax.set_xlim([insetB_1, insetB_2])

        ax.minorticks_on()
        ax.xaxis.tick_top()
        ax.xaxis.set_ticklabels([])
        ax.set_ylabel('Counts')

        ax.grid(linestyle=":", zorder=5)

        ax.axvline(insetC_1, linestyle="--", linewidth=1, color="blue", zorder=9)
        ax.axvline(insetC_2, linestyle="--", linewidth=1, color="blue", zorder=9)

        ## Band B
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.0,  0.74,  0.49, 0.06], [])
        ax.plot(time_arr, flux_B, linestyle='-', lw=1, color="red")
        ax.set_xlim([insetB_1, insetB_2])

        ax.minorticks_on()
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Counts')

        ax.grid(linestyle=":", zorder=5)

        ax.axvline(insetC_1, linestyle="--", linewidth=1, color="blue", zorder=9)
        ax.axvline(insetC_2, linestyle="--", linewidth=1, color="blue", zorder=9)


        ## --- Inset 3 ---
        ## Band A
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.51, 0.80, 0.49, 0.06], [])
        ax.plot(time_arr, flux_A, linestyle='-', lw=1, color="blue")
        ax.set_xlim([insetC_1, insetC_2])

        ax.minorticks_on()
        ax.xaxis.tick_top()
        ax.xaxis.set_ticklabels([])
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_ylabel('Counts')

        ax.grid(linestyle=":", zorder=5)

        ## Band B
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.51,  0.74,  0.49, 0.06], [])
        ax.plot(time_arr, flux_B, linestyle='-', lw=1, color="red")
        ax.set_xlim([insetC_1, insetC_2])

        ax.minorticks_on()
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Counts')

        ax.grid(linestyle=":", zorder=5)




        ## ----- Compilation: Y limits and Boxes -----


        ax = fig.add_axes([ 0.0,  0.66,  0.48, 0.02], [])
        ax.axis('off')
        ax.text(0.5, 1.0, "--- Inputs ---", horizontalalignment='center', verticalalignment='center', fontsize=14)

        ax = fig.add_axes([0.52,  0.66,  0.48, 0.02], [])
        ax.axis('off')
        ax.text(0.5, 1.0, "--- Outputs ---", horizontalalignment='center', verticalalignment='center', fontsize=14)

        ## --- Power Spectra ---
        ps_ylimits    = [
            min([
                 min(frequencies[1:len(power_spectrum_model_A)]*(power_spectrum_model_A[1:len(power_spectrum_model_A)]*frac_rms_norm_A
                  +poisson_noise_multiplier_A*poisson_noise_factor_A
                  +readout_noise_multiplier_A*read_noise_factor_A
                  +scintillation_noise_multiplier_A*scin_noise_factor_A)),
                 min(frequencies[1:len(power_spectrum_model_B)]*(power_spectrum_model_B[1:len(power_spectrum_model_B)]*frac_rms_norm_B
                  +poisson_noise_multiplier_B*poisson_noise_factor_B
                  +readout_noise_multiplier_B*read_noise_factor_B
                  +scintillation_noise_multiplier_B*scin_noise_factor_B)),
                 min(fourier_freqs*powerspectra_A),
                 min(fourier_freqs*powerspectra_B)
            ])*0.5,
            max([
                 max(frequencies[1:len(power_spectrum_model_A)]*(power_spectrum_model_A[1:len(power_spectrum_model_A)]*frac_rms_norm_A
                  +poisson_noise_multiplier_A*poisson_noise_factor_A
                  +readout_noise_multiplier_A*read_noise_factor_A
                  +scintillation_noise_multiplier_A*scin_noise_factor_A)),
                 max(frequencies[1:len(power_spectrum_model_B)]*(power_spectrum_model_B[1:len(power_spectrum_model_B)]*frac_rms_norm_B
                  +poisson_noise_multiplier_B*poisson_noise_factor_B
                  +readout_noise_multiplier_B*read_noise_factor_B
                  +scintillation_noise_multiplier_B*scin_noise_factor_B)),
                 max(fourier_freqs*powerspectra_A),
                 max(fourier_freqs*powerspectra_B)
            ])*1.5
            ]

        ps_rect = patches.Rectangle(
            ( min(fourier_freqs)*0.7
            , ps_ylimits[0])
            , (max(fourier_freqs)-min(fourier_freqs))*1.3
            , ps_ylimits[1]-ps_ylimits[0]
            , linewidth=1, linestyle='--', edgecolor='k', facecolor='none')

        ## --- Coherence ---
        coh_ylimits    = [
            min([
                min(coherence_model[1:len(coherence_model)]),
                min(coherence)
            ])*0.7,
            max([
                max(coherence_model[1:len(coherence_model)]),
                max(coherence)
            ])*1.3]

        coh_rect = patches.Rectangle(
            ( min(fourier_freqs)*0.7
            , coh_ylimits[0])
            , (max(fourier_freqs)-min(fourier_freqs))*1.3
            , coh_ylimits[1]-coh_ylimits[0]
            , linewidth=1, linestyle='--', edgecolor='k', facecolor='none')

        ## --- Phase Lags ---
        phase_rect = patches.Rectangle(
            ( min(fourier_freqs)*0.7
            , -3.8)
            , (max(fourier_freqs)-min(fourier_freqs))*1.3
            , 7.6
            , linewidth=1, linestyle='--', edgecolor='k', facecolor='none')

        ## --- Time Lags ---
        if max(time_lags_inv) > 0:
            time_ylimits    = [
                np.nanmin([
                    min([i for i in time_lags if i > 0] if max(time_lags) > 0 else [np.nan]),
                    min(i for i in time_lags_inv if i > 0)
                ])*0.7,
                np.nanmax([
                    max([i for i in time_lags if i > 0] if max(time_lags) > 0 else [np.nan]),
                    max(i for i in time_lags_inv if i > 0)
                ])*1.3]
        else:
            time_ylimits    = [
                min([
                    min(i for i in time_lags if i > 0),
                ])*0.7,
                max([
                    max(i for i in time_lags if i > 0),
                ])*1.3]

        # time_rect = patches.Rectangle(
        #     ( min(fourier_freqs)*0.7
        #     , time_ylimits[0])
        #     , (max(fourier_freqs)-min(fourier_freqs))*1.3
        #     , time_ylimits[1]-time_ylimits[0]
        #     , linewidth=1, linestyle='--', edgecolor='k', facecolor='none')

        ## ----- Compilation: Input -----


        ## --- CCF ---
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.0, 0.57, 0.48, 0.10], [])

        ax.plot(model_ccf_lags, model_ccf, color="k")

        ax.set_xlabel('Lag (s)')
        ax.set_ylabel('CCF Coefficient')

        ax.set_xlim(-ccf_segment/3, ccf_segment/3)

        ax.grid(linestyle=":", zorder=5)

        ax.axhline(0, linestyle="--", linewidth=1, zorder=9)
        ax.axvline(0, linestyle="--", linewidth=1, zorder=9)



        ## --- Power Spectra ---
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.0,  0.42,  0.48, 0.08], [])
        ax.plot(frequencies,  frequencies*
            (power_spectrum_model_A*frac_rms_norm_A
            +poisson_noise_multiplier_A*poisson_noise_factor_A
            +readout_noise_multiplier_A*read_noise_factor_A
            +scintillation_noise_multiplier_A*scin_noise_factor_A),
            linewidth=2, linestyle="-", label="Band A", color="blue")
        ax.plot(frequencies,  frequencies*
            (power_spectrum_model_B*frac_rms_norm_B
            +poisson_noise_multiplier_B*poisson_noise_factor_B
            +readout_noise_multiplier_B*read_noise_factor_B
            +scintillation_noise_multiplier_B*scin_noise_factor_B),
            linewidth=2, linestyle="-", label="Band B", color="red")

        ax.hlines(0, frequencies[1]/2, max(frequencies)*2, linestyle="--", linewidth=1, zorder=9)

        ax.xaxis.set_ticklabels([])
        ax.xaxis.tick_top()

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlim(frequencies[1]/2, max(frequencies)*2)
        ax.set_ylim(ps_ylimits)

        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Freq * Power')

        ax.grid(linestyle=":", zorder=5)

        ax.add_patch(ps_rect)


        ## --- Coherence ---
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.0, 0.34,  0.48, 0.08], [])
        ax.plot(frequencies,  coherence_model, linestyle="-", linewidth=2, color="darkviolet")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel('Coherence')
        ax.grid(linestyle=":", zorder=5)
        ax.set_ylim(coh_ylimits)
        ax.set_xlim(frequencies[1]/2, max(frequencies)*2)

        ax.add_patch(coh_rect)


        ## --- Phase & Time Lags ---
        if time_or_phase == 'time':
            model_phase_lags = []
            for j in range(len(model_lags)):
                model_phase_lags.append(model_lags[j] * 2 * np.pi * frequencies[j])


            ## --- Phase Lags ---
            ##                   left, base, wide, high
            ax = fig.add_axes([ 0.0,  0.26,  0.48, 0.08], [])

            ax.plot(frequencies,  model_phase_lags, linestyle="-", linewidth=2, color="darkviolet")
            ax.plot(frequencies,  np.array(model_phase_lags)+(4*np.pi), linestyle=":", linewidth=2, color="darkviolet")
            ax.plot(frequencies,  np.array(model_phase_lags)+(2*np.pi), linestyle=":", linewidth=2, color="darkviolet")
            ax.plot(frequencies,  np.array(model_phase_lags)-(2*np.pi), linestyle="--", linewidth=2, color="darkviolet")
            ax.plot(frequencies,  np.array(model_phase_lags)-(4*np.pi), linestyle="--", linewidth=1, color="darkviolet")
            ax.set_xlim(frequencies[1]/2, max(frequencies)*2)

            ax.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
            ax.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
            ax.yaxis.set_major_formatter(plt.FuncFormatter(corrfunc.format_func))
            ax.set_ylim([-3.8, 3.8])

            ax.set_ylabel('Phase Lags (rad)')
            ax.axhline(0, linestyle=":", linewidth=1, zorder=9)
            ax.set_xscale("log")
            ax.set_xlabel("Frequency (Hz)")
            ax.grid(linestyle=":", zorder=5)
            ax.add_patch(phase_rect)


            ## --- Time Lags ---
            ##                   left, base, wide, high
            ax_time = fig.add_axes([ 0.0, 0.18,  0.48, 0.08], [])

            ax_time.plot(frequencies,  model_lags, linestyle="-",  linewidth=2, color="darkviolet")
            ax_time.plot(frequencies, -model_lags, linestyle="--", linewidth=2, color="darkviolet")
            ax_time.set_xlim(frequencies[1]/2, max(frequencies)*2)
            ax_time.set_ylim([min(abs(model_lags))/2, max(np.append(model_lags, -model_lags))*2])
            ax_time.set_ylabel('Time Lags (s)')
            ax_time.set_yscale("log")
            ax_time.set_xscale("log")
            ax_time.set_xlabel("Frequency (Hz)")
            ax_time.grid(linestyle=":", zorder=5)
            # ax_time.add_patch(time_rect)

        elif time_or_phase == 'phase':

            ## --- Phase Lags ---
            ##                   left, base, wide, high
            ax = fig.add_axes([ 0.0,  0.26,  0.48, 0.08], [])

            ax.plot(frequencies, model_lags, linestyle="-", linewidth=2, color="darkviolet")
            ax.plot(frequencies,  model_lags+(4*np.pi), linestyle=":", linewidth=1, color="darkviolet")
            ax.plot(frequencies,  model_lags+(2*np.pi), linestyle=":", linewidth=2, color="darkviolet")
            ax.plot(frequencies,  model_lags-(2*np.pi), linestyle="--", linewidth=2, color="darkviolet")
            ax.plot(frequencies,  model_lags-(4*np.pi), linestyle="--", linewidth=1, color="darkviolet")
            ax.set_xlim(frequencies[1]/2, max(frequencies)*2)
            ax.axhline(0, linestyle=":", linewidth=1, zorder=9)

            ax.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
            ax.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
            ax.yaxis.set_major_formatter(plt.FuncFormatter(corrfunc.format_func))
            ax.set_ylim([-3.8, 3.8])

            ax.set_xscale("log")
            ax.set_xlabel("Frequency (Hz)")
            ax.set_ylabel('Phase Lag (rad)')
            ax.grid(linestyle=":", zorder=5)
            ax.add_patch(phase_rect)

            ## --- Time Lags ---
            ##                   left, base, wide, high
            ax_time = fig.add_axes([ 0.0, 0.18,  0.48, 0.08], [])

            ax_time.plot(frequencies[1:len(frequencies)],  model_lags[1:len(model_lags)]/(2*np.pi*frequencies[1:len(frequencies)])
                , linestyle="-",  linewidth=2, color="darkviolet")
            ax_time.plot(frequencies[1:len(frequencies)], -model_lags[1:len(model_lags)]/(2*np.pi*frequencies[1:len(frequencies)])
                , linestyle="--", linewidth=2, color="darkviolet")
            ax_time.set_xlim(frequencies[1]/2, max(frequencies)*2)
            ax_time.set_xscale("log")
            ax_time.set_yscale("log")
            ax_time.set_xlabel("Frequency (Hz)")
            ax_time.set_ylabel('Time Lag (s)')
            ax_time.grid(linestyle=":", zorder=5)
            # ax_time.add_patch(time_rect)


        ## ----- Compilation: Output -----


        ## --- Correlation Function ---

        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.52, 0.57, 0.48, 0.10], [])

        ax.step(ccf_av[:,0], ccf_av[:,1], where='mid', lw=1, color="black")
        ax.set_xlim(-ccf_segment/3, ccf_segment/3)

        ax.minorticks_on()
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_xlabel('Lag (s)')
        ax.set_ylabel('CCF Coefficient')

        ax.grid(linestyle=":", zorder=5)

        ax.axhline(0, linestyle="--", linewidth=1, zorder=9)
        ax.axvline(0, linestyle="--", linewidth=1, zorder=9)

        ax.annotate("CCF Range: " + str(ccf_segment) + "s", xy=(0.65, 0.55), xycoords="figure fraction")


        ## --- Power Spectra ---
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.52,  0.42,  0.48, 0.08], [])
        ax.errorbar(fourier_freqs, fourier_freqs*powerspectra_A, yerr=fourier_freqs*powerspectra_A_err \
            , ms=4, lw=0.5, capsize=3, fmt='bo', zorder=9)
        ax.errorbar(fourier_freqs, fourier_freqs*powerspectra_B, yerr=fourier_freqs*powerspectra_B_err \
            , ms=4, lw=0.5, capsize=3, fmt='ro', zorder=9)
        ax.set_ylim(ps_ylimits)

        ax.xaxis.set_ticklabels([])
        ax.xaxis.tick_top()

        ax.set_xscale('log')
        ax.set_ylabel('Freq * Power')
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_yscale('log', nonpositive='clip')
        ax.grid(linestyle='--')


        ## --- Coherence ---
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.52, 0.34,  0.48, 0.08], [])
        ax.set_xscale('log')
        ax.set_yscale('log', nonpositive="clip")
        ax.set_ylim(coh_ylimits)

        ax.errorbar(fourier_freqs, coherence, yerr=coherence_err, lw=0.5, fmt='mo', zorder=9)
        ax.set_ylabel("Coherence")
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.grid(linestyle='--')


        ## --- Phase Lags ---
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.52,  0.26,  0.48, 0.08], [])
        ax.set_xscale('log')

        ax.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(corrfunc.format_func))
        ax.set_ylim([-3.8, 3.8])

        ax.axhline(0, linestyle=":", linewidth=1, zorder=9)
        ax.errorbar(fourier_freqs, phase_lags, yerr=phase_lags_err, lw=0.5, fmt='mo')

        ax.set_ylabel("Phase Lag (rad)")
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.grid(linestyle='--')


        ## --- Time Lags ---
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.52, 0.18,  0.48, 0.08], [])
        ax.set_xscale('log')
        ax.set_yscale('log', nonpositive="clip")

        ax.errorbar(fourier_freqs, time_lags,     yerr=time_lags_err, lw=0.5, fmt='mo', label = "+ve Lags")
        ax.errorbar(fourier_freqs, time_lags_inv, yerr=time_lags_err, lw=0.5, linestyle="none", marker='o', fillstyle = "none", color="darkviolet", label = "-ve Lags")

        #ax.legend()

        ax.set_ylabel("Time Lag (s)")
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_xlabel("Frequency (Hz)")
        ax.grid(linestyle='--')

        timelags_xlim = ax.get_xlim()
        timelags_ylim = ax.get_ylim()

        time_rect = patches.Rectangle(
            ( timelags_xlim[0]
            , timelags_ylim[0])
            , timelags_xlim[1]-timelags_xlim[0]
            , timelags_ylim[1]-timelags_ylim[0]
            , linewidth=1, linestyle='--', edgecolor='k', facecolor='none')

        ax_time.add_patch(time_rect)


        ## --- Log File ---
        ##                   left, base, wide, high
        ax = fig.add_axes([ 0.0, 0.04,  1.0, 0.09], [])
        ax.axis('off')

        ax.text(0.0, 0.8, "Observation Length: " + str(obs_length) + "s")
        ax.text(0.5, 0.8, "Time Resolution: " + str(time_res) + "s")
        ax.text(0.0, 0.6, "Mean Counts/s (A): " + str(mean_counts_per_sec_A))
        ax.text(0.5, 0.6, "Mean Counts/s (B): " + str(mean_counts_per_sec_B))
        ax.text(0.0, 0.4, "F_rms (A): " + str(F_rms_A_orig))
        ax.text(0.5, 0.4, "F_rms (B): " + str(F_rms_B_orig))

        ax.text(0.0, 0.1, "Noise Sources:   Red: " + apply_red + "   Poisson: " + apply_poisson + \
            "   Scintillation: " + apply_scintillation + "   Readout: " + apply_readout)


        ## ----- Save Figure----

        savename = "./" + output_dir + "/" + fileprefix + "Compilation.png"
        plt.savefig(savename, bbox_inches='tight')

        plt.close('all')

    print("--------------- CorrSim Complete! ---------------")

if __name__ == "__main__":
    CorrSim()
