import math
import importlib
import numpy as np
import matplotlib.pyplot as plt
import CorrSim_Functions as corrfunc
import os
# import scipy

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
corrfunc.lag_section_phase(distribution, start_freq, start_lag, close_freq, close_lag=0, \
    extra_freq=0, extra_lag=0)
corrfunc.lag_section_time(distribution, start_freq, start_lag, close_freq, close_lag=0, \
    extra_freq=0, extra_lag=0)

    -Inputs:-

To define a lag, call corrfunc.lag_section_time or corrfunc.lag_section_time,
and define the distribution and the required parameters.

NOTE: Make sure to define your sections IN ORDER OF INCREASING FREQUENCY.
If you don't, it will create incorrect lags!


The following lag distributions are available:
    - Const. Time        = Constant lag in time
                            Creates phase lags that give a constant time (set by Start Lag)
                            Requires: start_freq, close_freq, start_lag

    - Const. Phase        = Constant lag in phase
                            Creates time lags that give a constant phase (set by Start Lag)
                            Requires: start_freq, close_freq, start_lag

    - Linear            = When time_or_phase="Phase": Linear in semi-log space.
                            When time_or_phase="Time": Linear in log-log space.
                            Requires: start_freq, close_freq, start_lag, close_lag

    - Power                = Exponentially increasing lag (Linear in normal space).
                            Requires: start_freq, close_freq, start_lag, close_lag

    - Polynomial        = When time_or_phase="Phase": Polynomial in semi-log space.
                            When time_or_phase="Time": Polynomial in log-log space.
                            Works by solving a second-order polynomial for the three coordinates given.
                            Requires: start_freq, close_freq, start_lag, close_lag, extra_freq, extra_lag



The output is seven numbers for each section. These are:
    - start        = Starting frequency of the section
    - close        = Ending frequency of the section
    - A            = Constant (see below)
    - B            = Constant (see below)
    - C            = Constant (see below)
    - P            = Constant (see below)
    - log        = Boolean; defines whether or not distributions have to be calculated using logs (1=Yes, 0=No)

    A, B, C, P are constants in the expression (Ax^2 + Bx^P + C)

'''



##    ===============================================================
##    ===============         Input Parameters        ===============


print("Define output directory (Default: 'CorrSim_Outputs'): ", end="")
output_dir = input()
if output_dir == "":
    output_dir = 'CorrSim_Outputs'

print("Define fileprefix directory (Default: 'CorrSim - '): ", end="")
fileprefix = input()
if fileprefix == "":
    fileprefix = 'CorrSim - '


print("Define length of observation (Default: 1000): ", end="")
length_of_observation = input()
if length_of_observation == "":
    length_of_observation = 1000

print("Define time resolution (Default: 0.1): ", end="")
time_resolution = input()
if time_resolution == "":
    time_resolution = 0.1

print("Which lags do you want to define, 'phase' or 'time'? (Default: 'phase'): ", end="")
time_or_phase = input()
if time_or_phase == "":
    time_or_phase = 'phase'

print("Define the overall lag (the lag for all non-defined sections) (Default: 0): ", end="")
overall_lag = input()
if overall_lag == "":
    overall_lag = 0
overall_lag = float(overall_lag)

print("How many lag sections do you want to define (max. 6)? ", end="")
num_sections = input()
while num_sections not in ["1","2","3","4","5","6"]:
    print("Please define a number between 1 and 6 inclusive:")
    num_sections = input()
num_sections = int(num_sections)

print("")
print("====================")
print("Beginning lag distributions.")
print("NOTE: Please define sections in order of ascending frequency!")
print("")


model_lag_array = np.zeros((6, 7))    ## Sets up the array. Don't Change!
distribution_list = ["Const. Time", "Const. Phase", "Linear", "Power", "Polynomial"]

for section in range(num_sections):
    print("----------")
    print("Section " + str(section+1) + ":")
    print("")
    print("Define Distribution")
    print("Const. Time / Const. Phase / Linear / Power / Polynomial")
    distribution = input()

    while distribution not in distribution_list:
        print("Distribution not found (or misttyped)! Please try again.")
        print("Const. Time / Const. Phase / Linear / Power / Polynomial")
        distribution = input()

    print("")
    if distribution == "Const. Time":
        print("Const. Time:")
        print("Defines a constant lag in time (set by Start Lag).")
        print("Requires: start_freq, close_freq, start_lag")

    if distribution == "Const. Phase":
        print("Const. Phase:")
        print("Defines a constant lag in phase (set by Start Lag).")
        print("Requires: start_freq, close_freq, start_lag")

    if distribution == "Linear":
        print("Linear:")
        print("A linear distribution between (start_freq, start_lag) and (close_freq, close_lag).")
        print("When time_or_phase='phase': Linear in semi-log space.")
        print("When time_or_phase='Time': Linear in log-log space.")
        print("Requires: start_freq, close_freq, start_lag, close_lag")

    if distribution == "Power":
        print("Power:")
        print("An exponentially increasing lag (Linear in normal space).")
        print("Requires: start_freq, close_freq, start_lag, close_lag")

    if distribution == "Polynomial":
        print("Polynomial:")
        print("A polynomial distribution.")
        print("When time_or_phase=\'Phase\': Polynomial in semi-log space.")
        print("When time_or_phase=\'Time\': Polynomial in log-log space.")
        print("Works by solving a second-order polynomial for the three coordinates given.")
        print("Requires: start_freq, close_freq, start_lag, close_lag, extra_freq, extra_lag")

    print("")

    ## Start Frequency (Freq_1)
    start_freq = None
    while type(start_freq) != float:
        try:
            print("Define starting frequency (start_freq): ", end="")
            start_freq = float(input())

            if start_freq <= 0:
                start_freq = None
                print("NOTE: The starting frequency must be positive.")
                continue

        except:
            print("NOTE: The starting frequency must be a positive, non-zero number.")
            pass

    ## Start Lag (Lag_1)
    start_lag = None
    while type(start_lag) != float:
        try:
            print("Define starting lag (start_lag): ", end="")
            start_lag = float(input())

            if time_or_phase == "time" and start_lag <= 0:
                start_lag = None
                print("NOTE: Lag distributions cannot go to, or across, 0 when defining time lags.")
                continue

        except:
            print("NOTE: The starting lag must be a number.")
            pass


    ## Start Frequency (Freq_2)
    close_freq = None
    while type(close_freq) != float:
        try:
            print("Define ending frequency (close_freq): ", end="")
            close_freq = float(input())

            if close_freq <= 0:
                close_freq = None
                print("NOTE: The ending frequency must be positive.")
                continue

            if close_freq <= start_freq:
                close_freq = None
                print("NOTE: The ending frequency must be greater than the starting frequency.")
                continue

        except:
            print("NOTE: The ending frequency must be a positive, non-zero number.")
            pass


    ## Close Lag (Lag_2)
    close_lag = 0
    if distribution != "Const. Time" and distribution != "Const. Phase":
        close_lag = None
        while type(close_lag) != float:
            try:
                print("Define ending lag (close_lag): ", end="")
                close_lag = float(input())

                if time_or_phase == "time" and start_lag <= 0:
                    close_lag = 0
                    print("NOTE: Lag distributions cannot go to, or across, 0 when defining time lags.")
                    continue

            except:
                print("NOTE: The ending lag must be a number.")
                pass

    extra_freq  = 0
    extra_lag   = 0
    if distribution == "Polynomial":

        ## Extra Frequency (Freq_3)
        extra_freq = None
        while type(extra_freq) != float:
            try:
                print("Define extra frequency (extra_freq): ", end="")
                extra_freq = float(input())
                if extra_freq <= 0:
                    extra_freq = None
                    print("NOTE: The extra frequency must be positive.")
                    continue
            except:
                print("NOTE: The extra frequency must be a positive, non-zero number.")
                pass

        ## Extra Lag (Lag_3)
        extra_lag = None
        while type(extra_lag) != float:
            try:
                print("Define extra lag (extra_lag): ", end="")
                extra_lag = float(input())
                if time_or_phase == "time" and start_lag <= 0:
                    extra_lag = None
                    print("NOTE: Lag distributions cannot go to, or across, 0 when defining time lags.")
                    continue
            except:
                print("NOTE: The extra lag must be a number.")
                pass

    print("")

    model_lag_array[section,] = corrfunc.lag_section_phase(distribution, start_freq, start_lag, close_freq, extra_freq, extra_lag)


print("====================")
print("")




## ----- Name -----

# output_dir                = 'CorrSim_Outputs'        ## Define output directory (Will be created if it doesn't exist)
# fileprefix                = 'CorrSim - Test - '    ## Define some prefix; this will be appended to all filenames.

## ----- Observation Parameters -----

## Set length of observation and time resolution (required to find the frequencies)
# length_of_observation    = 1027        ## Length of Observation in Seconds
# time_resolution            = 2**-5        ## Time Resolution in Seconds


## ========== Set Lag Parameters ==========

# ##    - Step 1: Choose whether to define time or phase lags.
# time_or_phase    = "phase"        ## "time"/"phase"
#
# model_lag_array = np.zeros((6, 7))    ## Sets up the array. Don't Change!
# if time_or_phase == "phase":
#     ##    - Step 2: Set the overall lag.
#     overall_lag = -4*np.pi/3            ## Choose the overall (constant) lag for all non-defined sections
#
#     ##    - Step 3: Set the lag model using corrfunc.lag_section_phase
#     ##                                                distribution    start_freq,    close_freq,     start_lag,     close_lag=0
#     model_lag_array[0,] = corrfunc.lag_section_phase("Const. Phase"    ,    0.001,  -4*np.pi/3,       0.02,   -4*np.pi/3)
#     model_lag_array[1,] = corrfunc.lag_section_phase("Power"            ,     0.02,  -4*np.pi/3,       0.25,    2*np.pi/5)
#     model_lag_array[2,] = corrfunc.lag_section_phase("Linear"        ,     0.25,   2*np.pi/5,        0.4,            0)
#     model_lag_array[3,] = corrfunc.lag_section_phase("Linear"        ,      0.4,           0,          5,        np.pi)
#     model_lag_array[4,] = corrfunc.lag_section_phase("Polynomial"    ,        5,       np.pi,        200,        np.pi, extra_freq = 28, extra_lag = 5*np.pi/2)
#
# elif time_or_phase == "time":
#     ##    - Step 2: Set the overall lag.
#     overall_lag        = -1E3        ## Choose the overall (constant) lag for all non-defined sections
#
#     ##    - Step 3: Set the lag model using corrfunc.lag_section_time.
#     ##    - NOTE: Time lags cannot go to or across zero. Please split distributions instead.
#     ##                                            distribution    start_freq,    close_freq,    start_lag,    close_lag=0
#     model_lag_array[0,] = corrfunc.lag_section_time("Linear"            ,    0.001,                -1E5,      0.01,         -1E1    )
#     model_lag_array[1,] = corrfunc.lag_section_time("Const. Time"    ,     0.01,          -1E1,         0.1                    )
#     model_lag_array[2,] = corrfunc.lag_section_time("Power"            ,      0.1,          -1E1,         0.4,         -1E-6    )
#     model_lag_array[3,] = corrfunc.lag_section_time("Linear"            ,      0.4,          1E-6,         0.8,          0.2    )
#     model_lag_array[4,] = corrfunc.lag_section_time("Linear"            ,      0.8,           0.2,           5,          0.1    )
#     model_lag_array[5,] = corrfunc.lag_section_time("Const. Phase"    ,        5,         np.pi,         200                 )
#


##    ===============            End of Input Parameters            ===============
##    =======================================================================


##                    Make Output Directory
##    --------------------------------------------------------

if not os.path.exists(output_dir):
    os.makedirs(output_dir)



## ===== Calculations =====

print("Creating model lags...")

## --- Define frequencies

obs_length    = length_of_observation
time_res    = time_resolution

freq_min    = 1 / obs_length                        ## Minimum Frequency
freq_max    = 1 / (2*time_res)                        ## Maximum Frequency
freq_res    = 1 / obs_length                        ## Frequency Resolution

frequencies    = np.arange(0, freq_max+(freq_res*0.1), freq_res)


## --- Calculate model lags

model_lags = np.zeros(len(frequencies)) + overall_lag            ## Define the basic model lags

section_pointer = 0                                                ## Create a variable that points to a section
for i in range(len(model_lags)):                                ## For each frequency bin...
    if frequencies[i] < model_lag_array[section_pointer, 0]:    ## If we're not yet in the current section,
        next                                                    ## Continue until we reach one.
    elif frequencies[i] > model_lag_array[section_pointer, 1]:    ## If we're past the current section.
        section_pointer = section_pointer + 1                    ## Increment to the next section

    if section_pointer == len(model_lag_array[:,1]):            ## If we're out of defined sections.
        break                                                    ## Break.
    elif frequencies[i] >= model_lag_array[section_pointer, 0] and frequencies[i] <= model_lag_array[section_pointer, 1]:    ## If we're in a section,
        if (model_lag_array[section_pointer, 6] == 1):                                                    ## If we need to use logs,
            model_lags[i] = model_lag_array[section_pointer, 2]*math.log(frequencies[i],10)**2 +    \
                model_lag_array[section_pointer, 3]*math.log(frequencies[i],10)**model_lag_array[section_pointer, 5] + \
                model_lag_array[section_pointer, 4]                                                        ## Calculate the lag based on the section's parameters.
        else:                                                                                    ## If we don't need to use logs,
            model_lags[i] = model_lag_array[section_pointer, 2]*frequencies[i]**2 + \
                model_lag_array[section_pointer, 3]*frequencies[i]**model_lag_array[section_pointer, 5] + \
                model_lag_array[section_pointer, 4]                                                        ## Calculate the lag based on the section's parameters.

if time_or_phase == 'time':
    print("Using: Time Lags")

    model_time_lags        = model_lags
    model_phase_lags    = np.array([])
    for j in range(len(model_lags)):
        model_phase_lags = np.append(model_phase_lags, model_lags[j] * 2 * np.pi * frequencies[j])

elif time_or_phase == 'phase':
    print("Using: Phase Lags")

    model_time_lags        = np.array([])
    model_phase_lags    = model_lags

    model_time_lags = np.append(model_time_lags, np.nan)
    for j in range(1, len(model_lags)):
        model_time_lags = np.append(model_time_lags, model_lags[j] / (2 * np.pi * frequencies[j]))

# model_time_lags_inv = [-x for x in model_time_lags]



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

ax.plot(frequencies, model_phase_lags-(4*np.pi),    lw=2, zorder=3, linestyle=':',  color="red", label = "-4π")
ax.plot(frequencies, model_phase_lags-(2*np.pi),    lw=2, zorder=3, linestyle='--', color="firebrick", label = "-2π")
ax.plot(frequencies, model_phase_lags,                 lw=2, zorder=3, linestyle='-',  color="black", label = "±0")
ax.plot(frequencies, model_phase_lags+(2*np.pi),    lw=2, zorder=3, linestyle='--', color="mediumblue", label = "+2π")
ax.plot(frequencies, model_phase_lags+(4*np.pi),    lw=2, zorder=3, linestyle=':',  color="dodgerblue", label = "+4π")

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

ax.plot(frequencies, model_time_lags,    linestyle="-",  lw=2, color="firebrick", label = "+ve Lags")
ax.plot(frequencies, -model_time_lags,    linestyle="--", lw=2, color="blue", label = "-ve Lags")

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
plt.savefig(fileOut, bbox_inches='tight')

plt.close('all')

print("Plotted model lags.")

print("")
print("model_lag_array (for use with CorrMain.CorrSim(); see Readme):")

# print(model_lag_array)
print("np." + repr(model_lag_array))
