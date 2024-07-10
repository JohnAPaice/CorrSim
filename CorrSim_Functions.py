## A selection of functions for use with CorrSim.

import warnings
## Ignore specific Stingray message:
warnings.filterwarnings("ignore", message="Large Datasets may not be processed efficiently due to computational constraints")

import copy
import math
import pandas
import numpy as np
import scipy.signal as ss
from stingray import Lightcurve, AveragedPowerspectrum, AveragedCrossspectrum
import time


##    ----------------------------------
##            General Functions
##    ----------------------------------

def largest_prime_factor(n):            ## Finds the largest prime factor of a number n
    """Find the largest prime factor of a number n"""
    i = 2
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
    return n

def place(value, list):                    ## Finds the (first) index of a given value in a list, or the closest to it.
    """Find the (first) index of a given value in a list, or the closest to it."""
    return min(range(len(list)), key=lambda i: abs(list[i]-value))

def ma(x, N):                            ## Creates a Moving Average of a list x, over N points.
    """Create a moving average of a list x, over N points"""
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def rev(lst):                            ## Reverses a list
    """Reverse a list"""
    return [ele for ele in reversed(lst)]

def ccf(x, y, lag=10000):                ## Creates an R-like Cross Correlation Function of two series (x, y)
    """
    Run an R-like Cross Correlation Function on two series, x and y.
    Uses scipy.signal as ss.

    Keyword parameters:
    x   - Series
    y   - Series
    lag - Integer. Maximum lag to investigate to. Default: 1000.
    """
    result = ss.correlate(y - np.mean(y), x - np.mean(x), method='direct') \
        / (np.std(y) * np.std(x) * len(y))
    length = (len(result) - 1) // 2
    lo = length - lag
    hi = length + (lag + 1)

    if lo < 0:        ## If the maximum lag goes negative, then reset to zero for intended behaviour.
        lo = 0

    return result[lo:hi]

def split_list(list, n):                ## Splits a list into segments of n length
                                        ## Usage: list(split_list(list,n))
    """
    Split a list into segments of n in length
    Usage: list(split_list(list,n))
    """
    for i in range(0, len(list), n):
        yield list[i:(i + n)]

def format_func(value, tick_number):    ## For plotting phases in units of pi/2
    """
    Plot phases in units of pi/2.
    Uses numpy as np.
    """
    N = int(np.round(2 * value / np.pi))
    if N == -2:
        return r"-$\pi$"
    elif N == -1:
        return r"-$\pi/2$"
    elif N == 0:
        return "0"
    elif N == 1:
        return r"$\pi/2$"
    elif N == 2:
        return r"$\pi$"
    elif N % 2 > 0:
        return r"${0}\pi/2$".format(N)
    else:
        return r"${0}\pi$".format(N // 2)

def statAv(rate,n):                        ## Binning a series every n points
    """Bin a series (rate) every n points"""

    binning_half    = int(np.floor(n/2))

    rate_binned     = []
    k               = binning_half + 1
    for i in range(int(np.floor(len(rate)/n))):
        while k < len(rate):
            rate_binned.append(np.mean(rate[(k - binning_half):(k + binning_half)]))
            k += n

    return rate_binned


def sort_list(list1, list2):            ## Sort list1 using an ordered list2
    """Sort list1 using an ordered list2"""
    zipped_pairs = zip(list2, list1)
    z = [x for _, x in sorted(zipped_pairs)]
    return z

def log_rebin(series, num=50):            ## Logarithmically rebin a series into a number of bins
    """Logarithmically rebin a series into a number of bins"""
    breaks = np.logspace(np.log10(series[1]), np.log10(series[-1]), num=num)
    cuts = pandas.cut(series, bins=breaks,labels=False)

    return breaks[0:-2], cuts[0:-2]

def log_rebin_average(series, cuts, breaks):    ## Bin a series into previously-create log bins
    """Bin a series into previously-create log bins"""
    series_binned = []
    for i in range(len(breaks)):
        k = np.where(cuts == i)
        if len(k[0]) == 0:
            if i > 0:
                series_binned.append(series_binned[i-1])
            else:
                series_binned.append(series[1])
        else:
            series_binned.append(np.mean(series[k]))
        if np.isfinite(series_binned[i]) == False:
            series_binned[i] = series_binned[i-1]

    return series_binned

def reverse_and_append(series):     ## Appends a reversed version of a series to the same series
    """Append a reversed version of a series to the same series"""
    series2 = np.append(series, rev(series[1:(len(series)-1)]))
    return series2

def complex_reverse_and_append(FTL_in):     ## Appends a reversed version of a series to the same series
    """Append the reversed complex conjugate of a complex series to that series"""
    FTL_out = np.append(FTL_in, rev(np.real(FTL_in)[1:(len(FTL_in)-1)] + (flipsign(np.imag(FTL_in)[1:(len(FTL_in)-1)]*1j))))
    return FTL_out

def flipsign(lst):                  ## Flips the sign of a list
    """Flip the sign of a list"""
    return [ -i for i in lst ]

##    ------------------------------------------
##            Log Distribution Functions
##    ------------------------------------------
## Find the given distributions in either log-log or semi-log space, from given coordinates.

def powerlaw(x1, x2, y1, y2):                                        ## --- Powerlaw
    """Find the normalisation and power of a powerlaw that connects two points."""
    slope   = math.log(y1/y2, x1/x2)                                ## Find the power that connects the two points,
    norm    = y1/(x1**slope)                                        ## and the normalisation.
    return [norm, slope]

def linear(x_1, y_1, x_2, y_2):                                        ## --- Linear
    """Find the slope and intercept of a straight line that connects two points."""
    slope        = (y_2-y_1) / (x_2-x_1)                                ## Find the slope between the two points,
    intercept    = y_2 - slope * x_2                                    ## and the intercept with the y axis.
    return [slope, intercept]

def semilog_linear(x_1, y_1, x_2, y_2, n):                            ## --- Linear (in semi-log space)
    """Find the slope and intercept of a straight line (in semi-log space) that connects two points."""
    slope        = (y_2-y_1) / ( math.log(x_2,n)-math.log(x_1,n) )    ## Find the slope between the two points,
    intercept    = y_2 - slope*math.log(x_2,n)                        ## and the intercept with the y axis.
    return [slope, intercept]

def log_linear(x_1, y_1, x_2, y_2, n):                                ## --- Linear (in log space)
    """Find the slope and intercept of a straight line (in semi-log space) that connects two points."""
    slope        = ( math.log(y_2,n)-math.log(y_1,n) ) / ( math.log(x_2,n)-math.log(x_1,n) )        ## Find the slope between the two points,
    intercept    = math.log(y_2,n) - slope*math.log(x_2,n)            ## and the intercept with the y axis.
    return [slope, intercept]

def semilog_polynomial(x_1, y_1, x_2, y_2, x_3, y_3, n):            ## --- Polynomial (in semi-log space)
    """Find the values of a polynomial (in semi-log space) that connects three points."""
    a    = [math.log(x_1,n), math.log(x_2,n), math.log(x_3,n)]        ## Find the three x-axis values...
    b    = [y_1, y_2, y_3]                                            ## Set the three y-axis values...
    solved_poly = np.polyfit(a,b,2)                                    ## And solve the polynomial for them.
    return solved_poly

def log_polynomial(x_1, y_1, x_2, y_2, x_3, y_3, n):                ## --- Polynomial (in log space)
    """Find the values of a polynomial (in log space) that connects three points."""
    a    = [math.log(x_1,n), math.log(x_2,n), math.log(x_3,n)]        ## Find the three x-axis values...
    b    = [math.log(y_1,n), math.log(y_2,n), math.log(y_3,n)]        ## Set the three y-axis values...
    solved_poly = np.polyfit(a,b,2)                                    ## And solve the polynomial for them.
    return solved_poly

## Lag Distribution function for time lags: Return the required distribution in the form Ax**2 + Bx**P + C
def lag_section_time(distribution, start_freq, start_lag, close_freq, \
                     close_lag=0, extra_freq=0, extra_lag=0):
    """
    For time lags.
    Takes input coordinates and a chosen distribution,
    and returns a required distribution in the form Ax**2 + Bx**P + C.

    Inputs:
    distribution    - The chosen shape of the line ('Constant', 'Linear', 'Polynomial', 'Power', 'Const. Phase')
    start_freq      - The starting frequency of the distribution (i.e. the first x coordinate).
    start_lag       - The starting lag of the distribution (i.e. the first y coordinate).
    close_freq      - The closing frequency of the distribution (i.e. the second x coordinate).
    close_lag       - The closing lag of the distribution (i.e. the second y coordinate). Default: 0
    extra_freq      - An extra frequency for Polynomial distributions (i.e. the third x coordinate).
    extra_lag       - An extra lag for Polynomial distributions (i.e. the third y coordinate).

    Outputs:
    start   - Starting frequency of the distribution.
    close   - Closing frequency of the distribution.
    A       - 'A' in Ax**2 + Bx**P + C
    B       - 'B' in Ax**2 + Bx**P + C
    C       - 'C' in Ax**2 + Bx**P + C
    P       - 'P' in Ax**2 + Bx**P + C
    log     - Boolean. Is this in log space? 1 (Yes) or 0 (No)

    """

    A = 0
    B = 0
    C = 0
    P = 1
    start = start_freq
    close = close_freq
    log = 0

    if distribution == "Constant" or distribution == "Const. Time":

        C = start_lag

    elif distribution == "Linear":
        phase_params = powerlaw(start_freq, close_freq, start_lag, close_lag)

        B = phase_params[0]
        P = phase_params[1]

    elif distribution == "Polynomial":

        phase_params = semilog_polynomial(
            start_freq, start_lag,
            extra_freq, extra_lag,
            close_freq, close_lag,
            10
        )

        A = phase_params[0]
        B = phase_params[1]
        C = phase_params[2]
        log = 1

        if extra_freq > close_freq:
            close = extra_freq


    elif distribution == "Power":

        phase_params = linear(start_freq, start_lag, close_freq, close_lag)

        B = phase_params[0]
        C = phase_params[1]

    elif distribution == "Const. Phase":

        B = start_lag/(2*np.pi)
        P = -1

    return [start, close, A, B, C, P, log]


## Lag Distribution function for phase lags: Return the required distribution in the form Ax**2 + Bx**P + C
def lag_section_phase(distribution, start_freq, start_lag, close_freq, \
                      close_lag=0, extra_freq=0, extra_lag=0):
    """
    For phase lags.
    Takes input coordinates and a chosen distribution,
    and returns a required distribution in the form Ax**2 + Bx**P + C.

    Inputs:
    distribution    - The chosen shape of the line ('Constant', 'Linear', 'Polynomial', 'Power', 'Const. Time')
    start_freq      - The starting frequency of the distribution (i.e. the first x coordinate).
    start_lag       - The starting lag of the distribution (i.e. the first y coordinate).
    close_freq      - The closing frequency of the distribution (i.e. the second x coordinate).
    close_lag       - The closing lag of the distribution (i.e. the second y coordinate). Default: 0
    extra_freq      - An extra frequency for Polynomial distributions (i.e. the third x coordinate).
    extra_lag       - An extra lag for Polynomial distributions (i.e. the third y coordinate).

    Outputs:
    start   - Starting frequency of the distribution.
    close   - Closing frequency of the distribution.
    A       - 'A' in Ax**2 + Bx**P + C
    B       - 'B' in Ax**2 + Bx**P + C
    C       - 'C' in Ax**2 + Bx**P + C
    P       - 'P' in Ax**2 + Bx**P + C
    log     - Boolean. Is this in log space? 1 (Yes) or 0 (No)

    """

    A = 0
    B = 0
    C = 0
    P = 1
    start = start_freq
    close = close_freq
    log = 0

    if distribution == "Constant" or distribution == "Const. Phase":

        C = start_lag

    elif distribution == "Linear":

        phase_params = semilog_linear(start_freq, start_lag, close_freq, close_lag, 10)

        B = phase_params[0]
        C = phase_params[1]
        log = 1

    elif distribution == "Polynomial":

        phase_params = semilog_polynomial(
            start_freq, start_lag,
            extra_freq, extra_lag,
            close_freq, close_lag,
            10
        )

        A = phase_params[0]
        B = phase_params[1]
        C = phase_params[2]
        log = 1

        if extra_freq > close_freq:
            close = extra_freq


    elif distribution == "Power":

        phase_params = linear(start_freq, start_lag, close_freq, close_lag)

        B = phase_params[0]
        C = phase_params[1]

    elif distribution == 'Const. Time':

        B = start_lag*2*np.pi

    return [start, close, A, B, C, P, log]



##    ----------------------------------
##            Main Functions
##    ----------------------------------


##            Create Power Spectra
##    ----------------------------------

## Creates power spectra using broken powerlaws
def power_model_broken_powerlaw(frequencies, ps_power, break_freq, \
                                coh_constant, coh_power, coh_break_freq):
    """
    Create two power spectra using a broken powerlaw model, combined with a broken coherence powerlaw.
    """

    coherence       = np.zeros(len(frequencies)) + coh_constant            ## Define constant section of Coherence
    coh_break_point = place(coh_break_freq, frequencies)                ## Find the break point for the Coherence
    coherence[coh_break_point:len(coherence)] =    \
        (frequencies[coh_break_point:len(coherence)] ** coh_power) * \
        (coh_constant / (frequencies[coh_break_point] ** coh_power))    ## Define and normalise the power law section of the Coherence

    unnorm_power_spectrum_A        = np.append(0, (frequencies[1:len((frequencies)+1)] ** ps_power))                    ## Define Band A PS...
    unnorm_power_spectrum_B_coh    = np.append(0, (frequencies[1:len((frequencies)+1)] ** ps_power)) * coherence        ## Define Coherent Band B PS...
    unnorm_power_spectrum_B_inc    = np.append(0, (frequencies[1:len((frequencies)+1)] ** ps_power)) * (1 - coherence)    ## Define Incoherent Band B PS...

    break_point = place(break_freq, frequencies)                        ## Define the break point for the power spectra
    constant    = frequencies[break_point] ** ps_power                    ## Find the value for the constant portion of the power spectrum

    unnorm_power_spectrum_A[    1:break_point]    = constant                ## Set the values below the break point to the constant (Band A)
    unnorm_power_spectrum_B_coh[1:break_point]    = (coherence[1:break_point]*0 + constant) * coherence[1:break_point]        ## Set the values below the break point as constant (based on coherence)
    unnorm_power_spectrum_B_inc[1:break_point]    = (coherence[1:break_point]*0 + constant) * (1 - coherence[1:break_point])    ## Set the values below the break point as constant (based on incoherence)

    return [unnorm_power_spectrum_A, unnorm_power_spectrum_B_coh, unnorm_power_spectrum_B_inc, coherence]

## Creates power spectra using lorentzians
def power_model_lorentzians(frequencies, Lorentz_params_A, Lorentz_params_B):
    """
    Create two power spectra using lists of Lorentzian parameters.

    NOTE:
    Lorentz_params_A must be divisible by 3 (Norm, Width, Mid), while
    Lorentz_params_B must be divisible by 4 (Norm, Width, Mid, Coherence)

    Inputs:
    - Frequencies       The range of frequencies
    - Lorentz_params_A  Parameters for Power Spectrum A.
    - Lorentz_params_B  Parameters for Power Spectrum B.
    """
    if (len(Lorentz_params_A)/3).is_integer() == False:
        raise Exception("Number of Lorentz Parameters in Series A must be divisible by 3.")
    if (len(Lorentz_params_B)/4).is_integer() == False:
        raise Exception("Number of Lorentz Parameters in Series B must be divisible by 4.")

    ## ========== Series A ==========

    ## Find the first instance of Normalisation = 0; this gives us how many Lorentzians we're dealing with.
    norms = np.arange(0, len(Lorentz_params_A), 3)
    # norms = [0, 3, 6, 9, 12, 15]
    lorentz_norms = [Lorentz_params_A[i] for i in norms]
    try:
        index_pos = lorentz_norms.index(0)
        num_lorentz_A = index_pos
    except ValueError as e:
        num_lorentz_A = int(len(Lorentz_params_A)/3)

    Lorentz_params_A = Lorentz_params_A[0:(num_lorentz_A*3)]

    # num_lorentz_A = int(len(Lorentz_params_A)/3)                                                                ## Find number of Lorentzians
    Lorentz_params_array_A = np.transpose(np.array(Lorentz_params_A).reshape(num_lorentz_A, 3))                    ## Read in the parameters
    Lorentzians_A = np.zeros((len(frequencies), num_lorentz_A))                                                    ## Create an array for the Lorentzians

    for i in range(num_lorentz_A):                                                                                ## For each Lorentzian...
        Lorentzians_A[:,i]    = Lorentz_params_array_A[0,i] * (1/np.pi) * (0.5 * Lorentz_params_array_A[1,i]) / \
            ((frequencies - Lorentz_params_array_A[2,i])**2 + (0.5 * Lorentz_params_array_A[1,i])**2)            ## Create it as a function of frequency

    unnorm_power_spectrum_A = np.zeros(len(frequencies))                                                        ## Set up the power spectrum...
    for i in range(1, len(frequencies)):
        unnorm_power_spectrum_A[i] = sum(Lorentzians_A[i,:])                                                    ## And sum all the Lorentzians to make it.


    ## ========== Series B ==========

    ## Find the first instance of Normalisation = 0; this gives us how many Lorentzians we're dealing with.
    norms = np.arange(0, len(Lorentz_params_B), 4)
    # norms = [0, 4, 8, 12, 16, 20]
    lorentz_norms = [Lorentz_params_B[i] for i in norms]
    try:
        index_pos = lorentz_norms.index(0)
        num_lorentz_B = index_pos
    except ValueError as e:
        num_lorentz_B = int(len(Lorentz_params_B)/4)

    Lorentz_params_B = Lorentz_params_B[0:(num_lorentz_B*4)]

    Lorentz_params_array_B = np.transpose(np.array(Lorentz_params_B).reshape(num_lorentz_B, 4))        ## Read in the parameters

    Lorentzians_B_coh = np.zeros((len(frequencies), num_lorentz_B))                                    ## Create an array for the coherent Lorentzians
    Lorentzians_B_inc = np.zeros((len(frequencies), num_lorentz_B))                                    ## Create an array for the incoherent Lorentzians

    for i in range(num_lorentz_B):                                                                    ## For each Lorentzian...
        Lorentzians_B_coh[:,i]          = ( Lorentz_params_array_B[0,i] * (1/np.pi) * (0.5 * Lorentz_params_array_B[1,i]) / \
            ((frequencies - Lorentz_params_array_B[2,i])**2 + (0.5 * Lorentz_params_array_B[1,i])**2) ) * Lorentz_params_array_B[3,i]        ## Create the coherent part...
        Lorentzians_B_inc[:,i]          = ( Lorentz_params_array_B[0,i] * (1/np.pi) * (0.5 * Lorentz_params_array_B[1,i]) / \
            ((frequencies - Lorentz_params_array_B[2,i])**2 + (0.5 * Lorentz_params_array_B[1,i])**2) ) * (1-Lorentz_params_array_B[3,i])    ## ...and the incoherent part.

    unnorm_power_spectrum_B_coh         = np.zeros(len(frequencies))                                        ## Set up the coherent power spectrum,
    unnorm_power_spectrum_B_inc         = np.zeros(len(frequencies))                                        ## and the incoherent one,
    coherence_model                     = np.zeros(len(frequencies))                                        ## as well as the model coherence.
    for i in range(1, len(frequencies)):
        unnorm_power_spectrum_B_coh[i]  = sum(Lorentzians_B_coh[i,:])                                ## Sum the coherent lorentzians,
        unnorm_power_spectrum_B_inc[i]  = sum(Lorentzians_B_inc[i,:])                                ## And the incoherent lorentzians
        # coherence_model[i]              = sum(Lorentzians_B_coh[i,:])/sum(Lorentzians_B_inc[i,:])    ## And divide for the coherence.
        coherence_model[i]              = sum(Lorentzians_B_coh[i,:])/(sum(Lorentzians_B_coh[i,:])+sum(Lorentzians_B_inc[i,:]))    ## And divide for the coherence.


    # print(coherence_model[(int(len(coherence_model)/2)-20):(int(len(coherence_model)/2)+20)])

    return [unnorm_power_spectrum_A, unnorm_power_spectrum_B_coh, unnorm_power_spectrum_B_inc, \
        coherence_model, Lorentzians_A, Lorentzians_B_coh, Lorentzians_B_inc, num_lorentz_A, num_lorentz_B]



##        Normalise Power Spectra
##    ----------------------------------

## Normalises power spectra; Uses F_rms normalisation.
def normalise_power_spectra(unnorm_power_spectrum_A, F_rms_A, \
                            unnorm_power_spectrum_B_coh, \
                            unnorm_power_spectrum_B_inc, F_rms_B, obs_length, \
                            time_res, num_bins, \
                            mean_counts_A, mean_counts_B):
    """Normalise both power spectra."""

    ## ========== Series A ==========

    normalisation_A         = ( F_rms_A**2 / sum(unnorm_power_spectrum_A / obs_length) )    ## Define normalisation
    power_spectrum_A        = unnorm_power_spectrum_A * normalisation_A                        ## Apply normalisation

    frac_rms_norm_A         = ((2 * time_res) / (num_bins * mean_counts_A**2))                ## Find F_rms normalisation (Vaughan 2003, Pg. 12)
    power_spectrum_A        = power_spectrum_A / frac_rms_norm_A                            ## And apply that too.


    ## ========== Series B ==========

    normalisation_B         = ( F_rms_B**2 /                                                ## Define normalisation
                                sum((unnorm_power_spectrum_B_coh + unnorm_power_spectrum_B_inc)
                                     / obs_length)
                                )
    power_spectrum_B_coh    = ( unnorm_power_spectrum_B_coh ) * normalisation_B                ## Apply normalisation to coherent power spectrum
    power_spectrum_B_inc    = ( unnorm_power_spectrum_B_inc ) * normalisation_B                ## and incoherent power spectrum

    frac_rms_norm_B            = ((2 * time_res) / (num_bins * mean_counts_B**2))                ## Find F_rms normalisation (Vaughan 2003, Pg. 12)
    power_spectrum_B_coh    = power_spectrum_B_coh / frac_rms_norm_B                        ## And apply that to the coherent...
    power_spectrum_B_inc    = power_spectrum_B_inc / frac_rms_norm_B                        ## ...and incoherent power spectra too.

    power_spectrum_B        = power_spectrum_B_coh + power_spectrum_B_inc                    ## Combine to make the final power spectrum for Series B.

    return [power_spectrum_A, power_spectrum_B, power_spectrum_B_coh, power_spectrum_B_inc, \
        normalisation_A, frac_rms_norm_A, normalisation_B, frac_rms_norm_B]

## Normalises power spectra; Uses F_rms normalisation. NEW: General version.
def normalise_single_power_spectra(unnorm_power_spectrum, F_rms, obs_length, \
                                   time_res, num_bins, mean_counts):
    """Normalise a single power spectra. Important for Red Noise calculations."""

    normalisation           = ( F_rms**2 / sum(unnorm_power_spectrum / obs_length) )    ## Define normalisation
    power_spectrum          = unnorm_power_spectrum * normalisation                        ## Apply normalisation

    frac_rms_norm           = ((2 * time_res) / (num_bins * mean_counts**2))                ## Find F_rms normalisation (Vaughan 2003, Pg. 12)
    power_spectrum          = power_spectrum / frac_rms_norm                            ## And apply that too.

    return [power_spectrum, normalisation, frac_rms_norm]



##         Apply Red Noise
##    ------------------------


def create_red_noise_real_imag_legacy(frequencies, power_spectrum, mean_counts, \
                               num_bins, norm, slope):
    """
    Create real and imaginary parts of a red noise power spectrum.
    LEGACY VERSION
    """

    rednoise_1      = np.random.normal(0, size=len(power_spectrum))                            ## First random number
    rednoise_2      = np.random.normal(0, size=len(power_spectrum))                            ## Second random number

    rednoise_1[1:]  = [norm * ((0.5*i)**(slope/2))*k for i, k in zip(frequencies[1:], rednoise_1[1:])]           ## Multiply for the first,
    rednoise_2[1:]  = [norm * ((0.5*i)**(slope/2))*k for i, k in zip(frequencies[1:], rednoise_2[1:])]           ## And for the second.
    rednoise_1[0]   = 0
    rednoise_2[0]   = 0

    # return(print("Woof"))
    FTL_real        = np.append(rednoise_1, rev(         rednoise_1[1:(len(rednoise_1)-1)]))    ## The first number is the real part of the Fourier-transformed data...
    FTL_imag        = np.append(rednoise_2, rev(flipsign(rednoise_2[1:(len(rednoise_2)-1)])))    ## and the second number is the imaginary part.

    # FTL_real = _ps + _rn

    FTL_imag[0]    = 0                                                                        ## Set the first imaginary bin as zero,
    if len(frequencies) % 2 == 1:                                                        ## And if it's an odd number of frequencies,
        FTL_imag[len(frequencies)-1]    = 0                                                ## the last freq. bin too.

    return [FTL_real, FTL_imag]

def create_red_noise_nums(power_spectrum):
    """Create two randomly-distributed numbers for the red noise"""
    rednoise_1      = np.random.normal(0, size=len(power_spectrum))                            ## First random number
    rednoise_2      = np.random.normal(0, size=len(power_spectrum))                            ## Second random number

    return [rednoise_1, rednoise_2]


def create_red_noise_real_imag(rednoise_1, rednoise_2, frequencies, mean_counts, \
                               num_bins, norm, slope):
    """Create real and imaginary parts of a red noise power spectrum."""

    rednoise_1a = copy.deepcopy(rednoise_1)
    rednoise_2a = copy.deepcopy(rednoise_2)

    rednoise_1a[1:]  = [norm * ((0.5*i)**(slope/2))*k for i, k in zip(frequencies[1:], rednoise_1a[1:])]           ## Multiply for the first,
    rednoise_2a[1:]  = [norm * ((0.5*i)**(slope/2))*k for i, k in zip(frequencies[1:], rednoise_2a[1:])]           ## And for the second.
    rednoise_1a[0]   = 0
    rednoise_2a[0]   = 0

    # return(print("Woof"))
    FTL_real        = np.append(rednoise_1a, rev(         rednoise_1a[1:(len(rednoise_1a)-1)]))    ## The first number is the real part of the Fourier-transformed data...
    FTL_imag        = np.append(rednoise_2a, rev(flipsign(rednoise_2a[1:(len(rednoise_2a)-1)])))    ## and the second number is the imaginary part.

    # FTL_real = _ps + _rn

    FTL_imag[0]    = 0                                                                        ## Set the first imaginary bin as zero,
    if len(frequencies) % 2 == 1:                                                        ## And if it's an odd number of frequencies,
        FTL_imag[len(frequencies)-1]    = 0                                                ## the last freq. bin too.

    return [FTL_real, FTL_imag]




##          Create Power Spectra from Fourier-Transformed Lightcurves
##    ---------------------------------------------------------------------

def create_ps_real_imag(frequencies, power_spectrum):
    """
    Create the real and imaginary parts of a complex
    Fourier-Transformed Lightcurve from a power spectrum.
    """

    Re = np.random.normal(0, scale = np.sqrt(power_spectrum/2), size=len(power_spectrum))                            ## First random number
    Im = np.random.normal(0, scale = np.sqrt(power_spectrum/2), size=len(power_spectrum))                            ## Second random number

    return [Re, Im]

def real_imag_to_ps(frequencies, FTL_real, FTL_imag):
    """
    Create a power spectrum from real/imaginary parts.
    """

    amplitude        = abs(FTL_real + FTL_imag*1j)                                        ## Find the amplitude by taking the absolute values,
    power_spectrum = amplitude[0:(len(frequencies))]**2                                    ## The power spectrum is the square of the amplitude!

    return(power_spectrum)


def ftl_to_lightcurve(FTL, correct_negative=True):
    """
    Create a lightcurve from a complex series.
    """

    ## Turn the Fourier-transformed lightcurves into actual lightcurves
    LC          = np.fft.ifft(FTL)         ## Inversely transform
    flux        = LC.real                  ## Find the flux

    ## Correct negative counts - assume nothing in that bin. WARNING: Messes up F_rms!
    if correct_negative == True:
        flux[flux<0] = 0

    return flux


def make_model_lags(frequencies, model_lag_array, overall_lag):
    """
    Make model lags from the input model_lag_array.
    """

    ## This will be calculating the lags from the lag parameters you set earlier.
    model_lags = np.zeros(len(frequencies)) + overall_lag            ## Define the basic model lags

    ## We'll be considering each section defined in turn.
    section_pointer = 0                                                ## Create a variable that points to a section
    for i in range(len(model_lags)):                                ## For each frequency bin...
        ## First, get into a defined section.
        if frequencies[i] < model_lag_array[section_pointer, 0]:    ## If we're not yet in the current section,
            next                                                    ## Continue until we reach it.
        elif frequencies[i] > model_lag_array[section_pointer, 1]:    ## If we're past the current section.
            section_pointer = section_pointer + 1                    ## Increment to the next section
        if section_pointer == len(model_lag_array[:,1])-1:            ## If we're out of defined sections.
            break                                                    ## Break.

        elif frequencies[i] >= model_lag_array[section_pointer, 0] and frequencies[i] <= model_lag_array[section_pointer, 1]:    ## If we're in a section,
            if (model_lag_array[section_pointer, 6] == 1):                                                                        ## If we need to use logs,
                model_lags[i] = model_lag_array[section_pointer, 2]*math.log(frequencies[i],10)**2 +    \
                    model_lag_array[section_pointer, 3]*math.log(frequencies[i],10)**model_lag_array[section_pointer, 5] + \
                    model_lag_array[section_pointer, 4]                                                                            ## Calculate the lag based on the section's parameters.
            else:                                                                                                                ## If we don't need to use logs,
                model_lags[i] = model_lag_array[section_pointer, 2]*frequencies[i]**2 + \
                    model_lag_array[section_pointer, 3]*frequencies[i]**model_lag_array[section_pointer, 5] + \
                    model_lag_array[section_pointer, 4]                                                                            ## Calculate the lag based on the section's parameters.

    return model_lags


##      Create Model CCF
##    ----------------------------------

def calc_model_CCF(coherence_model, power_spectrum_model_A, frac_rms_norm_A, \
                   power_spectrum_model_B, frac_rms_norm_B, model_lags, \
                   time_resolution, F_rms_A, F_rms_B):
    """
    Calculate model CCF.
    """

    ## We first make amplitudes and arguments...
    amp_1               = coherence_model * np.sqrt(power_spectrum_model_B*frac_rms_norm_B * power_spectrum_model_A*frac_rms_norm_A)
    amplitude_ccf       = np.sqrt(np.append(amp_1, rev(amp_1[1:(len(amp_1)-1)])))
    arguments_ccf       = np.append(model_lags, rev(-model_lags[1:(len(model_lags)-1)]))

    if len(arguments_ccf) % 2 == 1:
        arguments_ccf[int(arguments_ccf/2)]    = 0

    FT_ccf              = amplitude_ccf * np.cos(arguments_ccf) + 1j * amplitude_ccf * np.sin(arguments_ccf)    ## Find the Fourier Transformed version,
    IFT_ccf             = np.fft.ifft(FT_ccf)                                                                    ## Then inversely transform it...

    ccf_normalisation   = 1 / (2 * time_resolution * np.sqrt(F_rms_A * F_rms_B))                        ## Find the normalisation,

    model_ccf           = IFT_ccf.real * ccf_normalisation                                                        ## And calculate it!

    ## It's initially formatted wrongly about the middle - so cut in two, and stitch together
    model_ccf_former    = model_ccf[0:int(len(model_ccf)/2)]
    model_ccf_latter    = model_ccf[int(len(model_ccf)/2):len(model_ccf)]
    # model_ccf_former    = model_ccf[0:int(np.ceil(len(model_ccf)/2))]
    # model_ccf_latter    = model_ccf[int(np.ceil(len(model_ccf)/2)):len(model_ccf)]

    model_ccf           = rev(np.append(model_ccf_latter, model_ccf_former))

    ## Calculate the lag!
    lag = np.arange(1, len(model_ccf)+1)
    lag = lag - int(len(model_ccf)/2)
    lag = lag*time_resolution

    # lower = int(len(lag)/2)-20
    # upper = int(len(lag)/2)+20

    # print(lag[lower:upper])
    # print(amp_1[lower:upper])
    # print(amplitude_ccf[lower:upper])
    # print(arguments_ccf[lower:upper])
    # print(FT_ccf[lower:upper])
    # print(IFT_ccf[lower:upper])
    # print(ccf_normalisation)
    # print((IFT_ccf.real * ccf_normalisation)[lower:upper])
    # print(model_ccf_former[lower:upper])
    # print(model_ccf_latter[lower:upper])
    # print(model_ccf[lower:upper])

    return [lag, model_ccf]



##             Apply Noise
##    ---------------------------
## Various noise sources

def apply_poisson_noise(flux):
    """Add Poisson noise to a series."""
    return np.random.poisson(flux, len(flux))

def calculate_scintillation_noise(empirical_value, telescope_diameter, \
                                  exposure_time, target_altitude, \
                                  telescope_altitude, turbulence_height):
    """
    Calculate a value for the scintillation noise
    (To then use with apply_scintillation_noise)
    """
    zenith_distance = ( (90-target_altitude) * 2*np.pi) / 360    #Convert to Zenith Distance in Radians

    C_Y          = empirical_value
    D            = telescope_diameter
    t            = exposure_time
    gamma        = zenith_distance
    h_obs        = telescope_altitude
    H            = turbulence_height

    return 10 * 10**(-6) * C_Y**2 * D**(-4/3) * t**(-1) * np.cos(gamma)**(-3) * np.exp((-2*h_obs)/H)    ## Osborn+2015, Eqn. 7

def apply_scintillation_noise(scin_noise, flux):
    """
    Add Scintillation noise to a series.
    (Uses scin_noise from calculate_scintillation_noise)
    """
    return np.random.normal(loc = flux, scale = np.sqrt(scin_noise)*flux, size = len(flux))

def apply_readout_noise(readout_noise, flux):
    """Add normally-distributed noise to a series."""
    return flux + np.random.normal(loc = 0, scale = readout_noise, size = len(flux))


##             Calculate Average CCF
##    ----------------------------------------

def simulated_ccf(rateA, rateB, time_res, seconds, binning):
    """
    Calculate a Cross-Correlation Function between rateA and rateB,
    'seconds' in length, averaged over 'binning' bins.
    """

    ##            Binning (Optional)
    ## For if the data is high resolution, and you just want a quick look...
    if binning > 1:
        rateA_used        = statAv(rateA,binning)
        rateB_used        = statAv(rateB,binning)
        time_res_used    = time_res * binning
        print("Binning Complete.")

    else:
        rateA_used        = rateA
        rateB_used        = rateB
        time_res_used    = time_res


    ##            Segment Calculation
    ## Find the number of bins per segment
    bins_per_segment = int(np.floor(seconds/time_res_used))

    ## Check to make sure that the segment sizes are big enough
    if (bins_per_segment <= 10):
        print("WARNING: 10 or less points per segment in one band. Reduce binning or increase segment size.")
        #quit()

    ## Split up the data into segments of that size
    rateA_segments = list(split_list(rateA_used, bins_per_segment))
    rateB_segments = list(split_list(rateB_used, bins_per_segment))

    ## If the last segment does not have the right size, then discard
    if len(rateA_segments[0]) != len(rateA_segments[-1]):
        rateA_segments = rateA_segments[0:-1]
        rateB_segments = rateB_segments[0:-1]

    num_segments = len(rateA_segments)


    ##        Begin CCF of each section...
    ccf_arr = np.zeros(((bins_per_segment*2)-1,num_segments))
    for i in range(num_segments):

        ## From the full time range, select the data for this segment.
        ra = rateA_segments[i]
        rb = rateB_segments[i]

        ## Prewhiten the data - i.e. remove any long-term trends.
        pre_whiten = "(Pre-Whitened)"
        #pre_whiten = "n"
        if (pre_whiten == "(Pre-Whitened)"):
            av_range = 10
            rateA_start = np.mean(ra[0:av_range])
            rateA_close = np.mean(ra[(len(ra)-(av_range-1)):len(ra)])
            rateB_start = np.mean(rb[0:av_range])
            rateB_close = np.mean(rb[(len(rb)-(av_range-1)):len(rb)])

            ra2 = []
            rb2 = []
            fluxDiffA    = rateA_start-rateA_close
            fluxDiffB    = rateB_start-rateB_close
            for k in range(bins_per_segment):
                fraction    = (k)/(bins_per_segment)
                ra2.append(rateA_start-fluxDiffA*fraction)
                rb2.append(rateB_start-fluxDiffB*fraction)

            ra = (ra - ra2) + np.mean(ra)
            rb = (rb - rb2) + np.mean(rb)

        ## Calculate the CCF!
        ccf_arr[:,i] = ccf(ra, rb, bins_per_segment)

    ##     Post-Production - Average the results
    ccf_av = np.zeros(((bins_per_segment*2)-1,3))
    ccf_av[:,0] = np.arange(-bins_per_segment+1, bins_per_segment)*time_res_used
    for i in range((bins_per_segment*2)-1):
        ccf_av[i,1] = np.nanmean(ccf_arr[i,:])
        ccf_av[i,2] = np.nanstd( ccf_arr[i,:]) / np.sqrt(num_segments)

    return ccf_av


##             Calculate Fourier Data
##    --------------------------------------

def fourier(time, band_1, band_2, time_resolution, segment_size, rebin, \
            reference_freq, apply_poisson, remove_whitenoise):
    """
    Calculates all Fourier products:
Frequencies
    Power Spectra
    Coherence
    Phase
    Time Lags
    """

    ## --- Create Lightcurves ---
    lightcurve_A = Lightcurve(time, band_1, dt=time_resolution, skip_checks=True)
    lightcurve_B = Lightcurve(time, band_2, dt=time_resolution, skip_checks=True)


    ## --- Create Power Spectra ---
    avg_ps_A = AveragedPowerspectrum(lightcurve_A, segment_size, norm='frac')
    avg_ps_B = AveragedPowerspectrum(lightcurve_B, segment_size, norm='frac')


    ## --- Binning & White Noise ---
    powerspectra_binned_A = avg_ps_A.rebin_log(f=rebin)
    powerspectra_binned_B = avg_ps_B.rebin_log(f=rebin)

    if remove_whitenoise == 1:
        """
        Calculate white/poisson noise.
        This assumes that the last few points are white-noise dominated.
        The last point is usually based on only one bin, so we calculate it
        from the second-to-last point and assume the white noise is 99% of that.
        Note that this is a big approximation - proper fitting
        should be done for an actual calculation!
        """
        white_noise_A = np.mean(powerspectra_binned_A.power[-2])*0.99
        white_noise_B = np.mean(powerspectra_binned_B.power[-2])*0.99

    else:
        white_noise_A = 0
        white_noise_B = 0

    ## Remove the white noise (or not, if it's set to zero)
    powereadout_noiseless_binned_A = powerspectra_binned_A.power - white_noise_A
    powereadout_noiseless_binned_B = powerspectra_binned_B.power - white_noise_B

    ## --- Cross Spectra ---
    cross_spectrum = AveragedCrossspectrum(lightcurve_A, lightcurve_B, segment_size, norm='none')
    cross_spectrum_binned = cross_spectrum.rebin_log(f=rebin)

    ## --- Coherence ---
    coherence, err_coherence = cross_spectrum.coherence()                        ## Original
    coherence_binned, err_coherence_binned = cross_spectrum_binned.coherence()    ## Binned

    coherence_intrinsic_binned    = coherence_binned

    ## --- Coherence Errors ---
    ## Mainly from Vaughan & Nowak (1997) (Hereafter VN97) and racrs.pro
    ## Assuming High Powers and High Measured Coherence

    ## Use the same notation as VN97
    m   = powerspectra_binned_B.m    ## Number of Segments (# Error averaging from sections and rebinning)(?)
    N_1 = white_noise_A                                    ## White Noise (A); N^2 in VN97
    N_2 = white_noise_B                                    ## White Noise (B); N^2 in VN97
    S_1 = np.array(powereadout_noiseless_binned_A)        ## Noiseless Power Spectra (A);    S^2 in VN97
    S_2 = np.array(powereadout_noiseless_binned_B)        ## Noiseless Power Spectra (B); S^2 in VN97
    P_1 = np.array(powerspectra_binned_A.power)            ## Noisy Power Spectra (A)
    P_2 = np.array(powerspectra_binned_B.power)            ## Noisy Power Spectra (B)
    C   = np.array(coherence_binned * (S_1*S_2))        ## Complex-valued Cross Spectrum; C^2 in VN97
    n   = (S_1*N_2 + S_2*N_1 + N_1*N_2) / m                ## Simplification; n^2 in VN97

    ## --- New Error Calculations: (Done as close as possible to VN97 - 2021/08/04) ---
    err_1   = np.sqrt( 2 / m ) * (1-coherence_intrinsic_binned)/np.sqrt(coherence_intrinsic_binned)    ## Bendat & Piersol 1986
    err_1   = (m * err_1**2 ) / coherence_intrinsic_binned**2                                            ## Fourth term in the error equation

    err_coherence = ( (2 * n**2 * m) / ((C - n)**2) ) \
        + (N_1/S_1)**2 + (N_2/S_2)**2 \
        + err_1                                                                    ## Terms 1, 2+3,s 4

    err_coherence = np.sqrt(err_coherence)/np.sqrt(m)                            ## Deal with the m

    err_coherence_upper = ( (C - n) / (S_1 * S_2) ) * (1 + err_coherence)        ## Upper bound
    err_coherence_lower = ( (C - n) / (S_1 * S_2) ) * (1 - err_coherence)        ## Lower bound

    ## Average the upper and lower bounds for a simple (fractional) error value
    err_coherence_intrinsic_binned_frac = err_coherence_upper
    for i in range(len(err_coherence_upper)):
        err_coherence_intrinsic_binned_frac[i] = np.mean([abs(err_coherence_upper[i] - coherence_intrinsic_binned[i]), abs(err_coherence_lower[i] - coherence_intrinsic_binned[i])])

    err_coherence_intrinsic_binned = coherence_intrinsic_binned * err_coherence_intrinsic_binned_frac

    ## --- Freq Lags ---
    time_lags, time_lags_err = cross_spectrum_binned.time_lag()
    time_lags     = -time_lags
    time_lags_inv = -time_lags

    ## --- Phase Lags ---
    phase_lags                = []
    phase_lags_err            = []
    for lag_i in range(0, len(time_lags)):
        phase_lags.append(time_lags[lag_i] * 2 * np.pi * cross_spectrum_binned.freq[lag_i])
        phase_lags_err.append(time_lags_err[lag_i] * 2 * np.pi * cross_spectrum_binned.freq[lag_i])

    phase_lags2        = np.array(phase_lags)

    '''
    Phase lags can often go over the +/- pi boundary.
    In order to solve this problem, we assume a continuous distribution.
    First, a reference point is stated. From there, we work back and forwards.
    With each point, we look at its relation to the previous three points' mean.
    If it's closer to that value when adding/minusing 2pi,
        we do just that to that value and all others yet to be investigated.
    '''
    reference_point = place(reference_freq, powerspectra_binned_B.freq.real)

    ## To lower frequencies...
    for freq_i in range(reference_point-1,-1,-1):
        phase_diff = phase_lags2[freq_i] - np.mean([phase_lags2[(freq_i+1):(freq_i+3)]])
        if freq_i > 0:
            if phase_diff > np.pi:
                phase_lags2[0:freq_i+1] -= 2 * np.pi
            if phase_diff < -np.pi:
                phase_lags2[0:freq_i+1] += 2 * np.pi
        if freq_i == 0:
            if phase_diff > np.pi:
                phase_lags2[0] -= 2 * np.pi
            if phase_diff < -np.pi:
                phase_lags2[0] += 2 * np.pi

    ## To higher frequencies...
    for freq_i in range(reference_point+1,len(phase_lags2),1):
        phase_diff = phase_lags2[freq_i] - np.mean([phase_lags2[(freq_i-2):(freq_i)]])
        if phase_diff > np.pi:
            phase_lags2[freq_i:len(phase_lags2)] -= 2 * np.pi
        if phase_diff < -np.pi:
            phase_lags2[freq_i:len(phase_lags2)] += 2 * np.pi

    time_lags2 = []
    time_lags2_inv = []
    for freq_i in range(0, len(phase_lags2)):
        time_lags2.append(phase_lags2[freq_i] / ( 2 * np.pi * cross_spectrum_binned.freq[freq_i]))
        time_lags2_inv.append(-time_lags2[freq_i])

    time_lags_err[0] = segment_size

    ## Cut out the last bin; in testing, this always seems to be based on a single point after logarithmic binning.
    coherence_freqs            = powerspectra_binned_B.freq.real[0:-1]
    powerspectra_A            = powereadout_noiseless_binned_A[0:-1]
    powerspectra_A_err        = powerspectra_binned_A.power_err.real[0:-1]
    powerspectra_B            = powereadout_noiseless_binned_B[0:-1]
    powerspectra_B_err        = powerspectra_binned_B.power_err.real[0:-1]
    coherence                = coherence_intrinsic_binned[0:-1]
    coherence_err            = err_coherence_intrinsic_binned[0:-1]
    phase_lags                = phase_lags[0:-1]
    phase_lags_err            = phase_lags_err[0:-1]
    time_lags                = time_lags2[0:-1]
    time_lags_inv            = time_lags2_inv[0:-1]
    time_lags_err            = time_lags_err[0:-1]

    return [coherence_freqs, powerspectra_A, powerspectra_A_err, powerspectra_B, powerspectra_B_err,
        coherence, coherence_err, phase_lags, phase_lags_err, time_lags, time_lags_inv, time_lags_err]
