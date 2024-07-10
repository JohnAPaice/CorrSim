## 'Least I got chicken.

import CorrSim_Functions as corrfunc
import CorrSim
import CorrSim_TestLags_Function as corrlags_test
import CorrSim_TestPSDs_Function as corrpsds_test
import importlib
import numpy as np
import tkinter as tk
import tkinter.filedialog
from tkinter import Label, LEFT, SOLID, Toplevel, ttk



print("\n")
print("--------------- Launching CorrSim GUI... ---------------")
print("\n")

importlib.reload(corrfunc)
importlib.reload(CorrSim)
importlib.reload(corrlags_test)
importlib.reload(corrpsds_test)

class ToolTip(object):

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() +27
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = Label(tw, text=self.text, justify=LEFT,
                      background="#ffffe0", relief=SOLID, borderwidth=1,
                      font=("Arial", "9", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)

class App(tk.Tk):
    def __init__(self):
        super().__init__()

        ### ===================================
        ### --------------- GUI ---------------

        self.title("CorrSim")
        ## Set Width and Height
        self.geometry("1115x850")

        ## ========== Title Image ==========

        self.title_image = tk.PhotoImage(file = "CorrSim_Title.png")
        # title_image = title_image.zoom(5)
        # title_image = title_image.subsample(10)
        canvas1 = tk.Canvas(self, width=900, height=100)
        canvas1.grid(row=0, column=0)
        canvas1.create_image(0,0, image=self.title_image, anchor="nw")


        ## ========== Menu Bar ==========

        self.menubar = tk.Menu(self)

        self.filemenu = tk.Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label="Load", command=self.loadfile)
        self.filemenu.add_command(label="Save", command=self.savefile)
        self.filemenu.add_command(label="Save as...", command=self.savefileas)

        self.menubar.add_cascade(label="File", menu=self.filemenu)

        over_row_num    = 1
        mid_row_num     = 0
        row_num         = 0


        ## =============================================================================
        ## ==================== Observational and Source Parameters ====================

        row_num              = 0
        self.frm_inputs      = tk.Frame(master=self, padx=5, pady=0)
        self.frm_inputs.grid(row=over_row_num, column=0, padx=0, pady=0, sticky="nsew")
        over_row_num         +=1

        self.frm_inputs_left = tk.Frame(master=self.frm_inputs, padx=0, pady=0)
        self.frm_inputs_left.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")

        row_num              = 0
        self.frm_main_params = tk.LabelFrame(master=self.frm_inputs_left, text="Observation and Source Parameters", font=('bitstream vera sans', 9, 'bold'), padx=3, pady=0)
        self.frm_main_params.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")
        over_row_num         +=1

        ## ========== Observation Length and Number of Bins ==========

        ## Set up the frame.
        row_num              = 0
        self.frm_numbins     = tk.Frame(master=self.frm_main_params)
        self.frm_numbins.grid(row=mid_row_num, column=0, padx=0, pady=1)
        mid_row_num          +=1


        ## Observation length
        self.obs_length      = tk.StringVar()
        self.lbl_obs_length  = tk.Label(master=self.frm_numbins, text="Observation Length (s):")
        self.lbl_obs_length.grid(row=row_num, column=0, sticky="e")
        self.ent_obs_length  = tk.Entry(master=self.frm_numbins, width=10, textvariable=self.obs_length)
        self.ent_obs_length.grid(row=row_num, column=1)
        row_num              += 1
        CreateToolTip(self.lbl_obs_length, text =
            "Length of simulated observation (seconds). \n"
            "Particularly relevant for the number of bins,\n"
            "which is the biggest factor in the size of\n"
            "the output files and the speed of CorrSim.")


        ## Time Resolution
        self.time_res        = tk.StringVar()
        self.lbl_time_res    = tk.Label(master=self.frm_numbins, text="Time Resolution (s):")
        self.lbl_time_res.grid(row=row_num, column=0, sticky="e")
        self.ent_time_res    = tk.Entry(master=self.frm_numbins, width=10, textvariable=self.time_res)
        self.ent_time_res.grid(row=row_num, column=1)
        row_num              += 1
        CreateToolTip(self.lbl_time_res, text =
            "Time resolution of simulated observation (seconds).\n"
            "Particularly relevant for the number of bins,\n"
            "which is the biggest factor in the size of\n"
            "the output files and the speed of CorrSim.")

        ## Optimise
        self.optimise          = tk.IntVar()    ## Initialise the optimise keyword
        self.chk_optimise      = tk.Checkbutton(master=self.frm_numbins, text="Optimise Bins", variable=self.optimise)
        self.chk_optimise.grid(row=0, column=2)
        row_num                += 1
        CreateToolTip(self.chk_optimise, text =
            "Optional; Adjusts the observation length so that\n"
            "the number of bins has no factor bigger than 7.\n"
            "This can vastly speed up calculations.")

        ## Calculate Bins
        self.btn_calculate        = tk.Button(master=self.frm_numbins, text="Calculate Bins", command=self.calculate_num_bins)
        self.btn_calculate.grid(row=row_num, column=0)

        ## Number of Bins Printout
        self.lbl_numbin_result    = tk.Label(master=self.frm_numbins, text="")
        self.lbl_numbin_result.grid(row=row_num, column=1)

        ## Number of Bins Warning
        self.lbl_numbin_wrning       = tk.Label(master=self.frm_numbins, text="")
        self.lbl_numbin_wrning["fg"] = "red"
        self.lbl_numbin_wrning.grid(row=row_num, column=2)
        row_num                      += 1
        CreateToolTip(self.lbl_numbin_wrning, text =
            "For large numbers of bins (>1e6), CorrSim can take\n"
            "a long time to run, and the output files can be \n"
            "very large (>100 MB). Be aware of the size of\n"
            "the files you are creating!")

        ttk.Separator(self.frm_main_params, orient='horizontal').grid(row=mid_row_num, column=0, columnspan=3, sticky='nsew', pady=0)
        mid_row_num += 1

        ## ========== Source Parameters ==========

        ## Set up the frame.
        row_num                = 0
        self.frm_source        = tk.Frame(master=self.frm_main_params)
        self.frm_source.grid(row=mid_row_num, column=0, padx=0, pady=1)
        mid_row_num_save = mid_row_num

        self.mean_counts_A     = tk.StringVar()
        self.mean_counts_B     = tk.StringVar()
        self.F_rms_A           = tk.StringVar()
        self.F_rms_B           = tk.StringVar()

        tk.Label(master=self.frm_source, text="A").grid(row=row_num, column=1, ipadx=10)
        tk.Label(master=self.frm_source, text="B").grid(row=row_num, column=2, ipadx=10)
        row_num                += 1

        ## Mean Counts/s
        self.lbl_mean_counts      = tk.Label(master=self.frm_source, text="Mean Counts/s:")
        self.lbl_mean_counts.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_mean_counts_A    = tk.Entry(master=self.frm_source, width=10, textvariable=self.mean_counts_A)
        self.ent_mean_counts_A.grid(row=row_num, column=1)
        self.ent_mean_counts_B    = tk.Entry(master=self.frm_source, width=10, textvariable=self.mean_counts_B)
        self.ent_mean_counts_B.grid(row=row_num, column=2)
        row_num                   += 1
        CreateToolTip(self.lbl_mean_counts, text = "Mean count rate (Counts/s)")

        ## Fractional RMS
        self.lbl_F_rms            = tk.Label(master=self.frm_source, text="F_rms:")
        self.lbl_F_rms.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_F_rms_A          = tk.Entry(master=self.frm_source, width=10, textvariable=self.F_rms_A)
        self.ent_F_rms_A.grid(row=row_num, column=1)
        self.ent_F_rms_B          = tk.Entry(master=self.frm_source, width=10, textvariable=self.F_rms_B)
        self.ent_F_rms_B.grid(row=row_num, column=2)
        row_num                   += 1
        CreateToolTip(self.lbl_F_rms, text = "Fractional RMS")

        ttk.Separator(self.frm_source, orient='horizontal').grid(row=row_num, column=1, columnspan=2, sticky='nsew', pady=1)
        row_num += 1


        ## ========== Noise Sources ==========

        self.apply_red_A            = tk.IntVar()
        self.apply_red_B            = tk.IntVar()
        self.apply_poisson_A        = tk.IntVar()
        self.apply_poisson_B        = tk.IntVar()
        self.F_rms_rednoise_A       = tk.StringVar()
        self.F_rms_rednoise_B       = tk.StringVar()
        self.rednoise_slope_A       = tk.StringVar()
        self.rednoise_slope_B       = tk.StringVar()


        ## ----- Red Noise -----

        self.chk_apply_red_A        = tk.Checkbutton(master=self.frm_source, variable=self.apply_red_A, command=self.red_enable_disable_A)
        self.chk_apply_red_B        = tk.Checkbutton(master=self.frm_source, variable=self.apply_red_B, command=self.red_enable_disable_B)
        self.lbl_red_noise          = tk.Label(master=self.frm_source, text="Apply Red Noise")
        self.lbl_red_noise.grid(row=row_num, column=0, sticky="e", ipadx=10, ipady=0)
        self.chk_apply_red_A.grid(row=row_num, column=1)
        self.chk_apply_red_B.grid(row=row_num, column=2)
        row_num                     += 1
        CreateToolTip(self.lbl_red_noise, text =
            "Red noise is also known as ‘Brownian noise’,\n"
            "and is proportional to frequency^slope.\n"
            "CorrSim uses the methods of Timmer & Koenig 1995\n"
            "(1995A&A...300..707T); For each Fourier Frequency,\n"
            "two Gaussian distributed random numbers are drawn, and\n"
            "multiplied by sqrt(0.5*(power_spectrum)). The result is\n"
            "the real & imaginary part of the Fourier transformed data.")


        ## Define fractional rms of the red noise (0-1)
        self.lbl_red_frms        = tk.Label(master=self.frm_source, text="Red Noise F_rms:")
        self.lbl_red_frms.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_red_frms_A      = tk.Entry(master=self.frm_source, width=10, textvariable=self.F_rms_rednoise_A)
        self.ent_red_frms_A.grid(row=row_num, column=1)
        self.ent_red_frms_B      = tk.Entry(master=self.frm_source, width=10, textvariable=self.F_rms_rednoise_B)
        self.ent_red_frms_B.grid(row=row_num, column=2)
        row_num                  +=1
        CreateToolTip(self.lbl_red_frms, text =
            "Fractional RMS of the Red Noise\n"
            "(Typical values are 0.03-0.2)\n"
            "Must be smaller than F_rms in each band.")

        ## Define slope of the red noise
        self.lbl_red_frms        = tk.Label(master=self.frm_source, text="Red Noise Slope:")
        self.lbl_red_frms.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_red_slope_A     = tk.Entry(master=self.frm_source, width=10, textvariable=self.rednoise_slope_A)
        self.ent_red_slope_A.grid(row=row_num, column=1)
        self.ent_red_slope_B     = tk.Entry(master=self.frm_source, width=10, textvariable=self.rednoise_slope_B)
        self.ent_red_slope_B.grid(row=row_num, column=2)
        row_num                  +=1
        CreateToolTip(self.lbl_red_frms, text =
            "Slope of the Red Noise\n"
            "(Typically -2)")


        ttk.Separator(self.frm_source, orient='horizontal').grid(row=row_num, column=1, columnspan=2, sticky='nsew', pady=1)
        row_num                += 1


        ## ----- Poisson Noise -----

        self.chk_apply_poisson_A  = tk.Checkbutton(master=self.frm_source, variable=self.apply_poisson_A)
        self.chk_apply_poisson_B  = tk.Checkbutton(master=self.frm_source, variable=self.apply_poisson_B)
        self.lbl_poisson_noise    = tk.Label(master=self.frm_source, text="Apply Poisson Noise")
        self.lbl_poisson_noise.grid(row=row_num, column=0, sticky="e", ipadx=10, ipady=0)
        self.chk_apply_poisson_A.grid(row=row_num, column=1)
        self.chk_apply_poisson_B.grid(row=row_num, column=2)
        row_num                   += 1
        CreateToolTip(self.lbl_poisson_noise, text =
            "Poisson noise, also known as ‘shot noise’, comes\n"
            "from the random nature of emitted photons.\n"
            "CorrSim applies poisson noise to the lightcurves by\n"
            "drawing a random number from a poisson distribution\n"
            "for each bin, with the mean centred on the flux.")


        ## ----- Read Noise -----

        ttk.Separator(self.frm_source, orient='horizontal').grid(row=row_num, column=1, columnspan=2, sticky='nsew', pady=1)
        row_num += 1

        self.apply_readout_A         = tk.IntVar()
        self.apply_readout_B         = tk.IntVar()
        self.read_noise_A            = tk.StringVar()
        self.read_noise_B            = tk.StringVar()

        ## Apply Readout Noise
        self.lbl_apply_readout       = tk.Label(master=self.frm_source, text="Apply Readout Noise")
        self.lbl_apply_readout.grid(row=row_num, column=0, sticky="e", ipadx=10, ipady=0)
        self.chk_apply_readout_A     = tk.Checkbutton(master=self.frm_source, variable=self.apply_readout_A, command=self.read_enable_disable_A)
        self.chk_apply_readout_A.grid(row=row_num, column=1)
        self.chk_apply_readout_B     = tk.Checkbutton(master=self.frm_source, variable=self.apply_readout_B, command=self.read_enable_disable_B)
        self.chk_apply_readout_B.grid(row=row_num, column=2)
        row_num                      += 1
        CreateToolTip(self.lbl_apply_readout, text =
            "Readout noise comes from fluctuations in\n"
            "reading out charge from a CCD. This is \n"
            "dependant on the Readout Noise parameters.")

        ## Define amount of readout noise (electrons per bin)
        self.lbl_readout_noise       = tk.Label(master=self.frm_source, text="Readout Noise (e-):")
        self.lbl_readout_noise.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_read_noise_A        = tk.Entry(master=self.frm_source, width=10, textvariable=self.read_noise_A)
        self.ent_read_noise_A.grid(row=row_num, column=1)
        self.ent_read_noise_B        = tk.Entry(master=self.frm_source, width=10, textvariable=self.read_noise_B)
        self.ent_read_noise_B.grid(row=row_num, column=2)
        row_num                      +=1
        CreateToolTip(self.lbl_readout_noise, text =
            "Readout noise from a CCD in electrons.")

        ## ----- Scintillation Noise -----

        ttk.Separator(self.frm_source, orient='horizontal').grid(row=row_num, column=1, columnspan=2, sticky='nsew', pady=1)
        row_num += 1

        self.apply_scintillation_A    = tk.IntVar()
        self.apply_scintillation_B    = tk.IntVar()


        ## Apply Scintillation Noise
        self.lbl_apply_scintillation      = tk.Label(master=self.frm_source, text="Apply Scintillation Noise")
        self.lbl_apply_scintillation.grid(row=row_num, column=0, sticky="e", ipadx=10, ipady=0)
        self.chk_apply_scintillation_A    = tk.Checkbutton(master=self.frm_source, variable=self.apply_scintillation_A, command=self.scin_enable_disable_A)
        self.chk_apply_scintillation_A.grid(row=row_num, column=1)
        self.chk_apply_scintillation_B    = tk.Checkbutton(master=self.frm_source, variable=self.apply_scintillation_B, command=self.scin_enable_disable_B)
        self.chk_apply_scintillation_B.grid(row=row_num, column=2)
        row_num                           += 1
        CreateToolTip(self.lbl_apply_scintillation, text =
            "This noise comes from atmospheric effects, and is thus\n"
            "only applicable to simulations from ground-based observatories.")

        ## -- Define Scintillation Noise Parameters --
        self.telescope_diameter_A        = tk.StringVar()
        self.telescope_altitude_A        = tk.StringVar()
        self.exposure_time_A             = tk.StringVar()
        self.turbulence_height_A         = tk.StringVar()
        self.target_altitude_A           = tk.StringVar()
        self.empirical_value_A           = tk.StringVar()
        self.telescope_diameter_B        = tk.StringVar()
        self.telescope_altitude_B        = tk.StringVar()
        self.exposure_time_B             = tk.StringVar()
        self.turbulence_height_B         = tk.StringVar()
        self.target_altitude_B           = tk.StringVar()
        self.empirical_value_B           = tk.StringVar()

        ## Telescope Diameter (m)
        self.lbl_telescope_diameter      = tk.Label(master=self.frm_source, text="Telescope Diameter (m):")
        self.lbl_telescope_diameter.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_telescope_diameter_A    = tk.Entry(master=self.frm_source, width=10, textvariable=self.telescope_diameter_A)
        self.ent_telescope_diameter_A.grid(row=row_num, column=1)
        self.ent_telescope_diameter_B    = tk.Entry(master=self.frm_source, width=10, textvariable=self.telescope_diameter_B)
        self.ent_telescope_diameter_B.grid(row=row_num, column=2)
        row_num                          +=1
        CreateToolTip(self.lbl_telescope_diameter, text =
            "Diameter of the telescope (metres).\n"
            "Increasing this value decreases the scintillation noise.")

        ## Telescope Altitude (m)
        self.lbl_telescope_altitude      = tk.Label(master=self.frm_source, text="Telescope Altitude (m):")
        self.lbl_telescope_altitude.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_telescope_altitude_A    = tk.Entry(master=self.frm_source, width=10, textvariable=self.telescope_altitude_A)
        self.ent_telescope_altitude_A.grid(row=row_num, column=1)
        self.ent_telescope_altitude_B    = tk.Entry(master=self.frm_source, width=10, textvariable=self.telescope_altitude_B)
        self.ent_telescope_altitude_B.grid(row=row_num, column=2)
        row_num                          +=1
        CreateToolTip(self.lbl_telescope_altitude, text =
            "Altitude of the observing telescope (metres).\n"
            "Increasing this value decreases the scintillation noise.")

        ## Exposure Time (s)
        self.lbl_exposure_time            = tk.Label(master=self.frm_source, text="Exposure Time (s):")
        self.lbl_exposure_time.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_exposure_time_A          = tk.Entry(master=self.frm_source, width=10, textvariable=self.exposure_time_A)
        self.ent_exposure_time_A.grid(row=row_num, column=1)
        self.ent_exposure_time_B          = tk.Entry(master=self.frm_source, width=10, textvariable=self.exposure_time_B)
        self.ent_exposure_time_B.grid(row=row_num, column=2)
        row_num                           +=1
        CreateToolTip(self.lbl_exposure_time, text =
            "Exposure time of the observing telescope. This is\n"
            "different from time resolution, as it is does NOT include\n"
            "the time taken to read out information from the CCD ('Deadtime').\n"
            "Increasing this value decreases the scintillation noise.")

        ## Turbulence Height (m)
        self.lbl_turbulence_height        = tk.Label(master=self.frm_source, text="Turbulence Height (m):")
        self.lbl_turbulence_height.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_turbulence_height_A      = tk.Entry(master=self.frm_source, width=10, textvariable=self.turbulence_height_A)
        self.ent_turbulence_height_A.grid(row=row_num, column=1)
        self.ent_turbulence_height_B      = tk.Entry(master=self.frm_source, width=10, textvariable=self.turbulence_height_B)
        self.ent_turbulence_height_B.grid(row=row_num, column=2)
        row_num                           +=1
        CreateToolTip(self.lbl_turbulence_height, text =
            "Height of turbulence in the atmosphere (metres).\n"
            "A typical value here is around 8000.\n"
            "Increasing this value increases the scintillation noise.")

        ## Target Altitude (Degrees)
        self.lbl_target_altitude          = tk.Label(master=self.frm_source, text="Target Altitude (º):")
        self.lbl_target_altitude.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_target_altitude_A        = tk.Entry(master=self.frm_source, width=10, textvariable=self.target_altitude_A)
        self.ent_target_altitude_A.grid(row=row_num, column=1)
        self.ent_target_altitude_B        = tk.Entry(master=self.frm_source, width=10, textvariable=self.target_altitude_B)
        self.ent_target_altitude_B.grid(row=row_num, column=2)
        row_num                           +=1
        CreateToolTip(self.lbl_target_altitude, text =
            "Altitude of the source (degrees).\n"
            "Increasing this value decreases the scintillation noise.")

        ## Empirical Value
        self.lbl_empirical_value          = tk.Label(master=self.frm_source, text="Empirical Value:")
        self.lbl_empirical_value.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_empirical_value_A        = tk.Entry(master=self.frm_source, width=10, textvariable=self.empirical_value_A)
        self.ent_empirical_value_A.grid(row=row_num, column=1)
        self.ent_empirical_value_B        = tk.Entry(master=self.frm_source, width=10, textvariable=self.empirical_value_B)
        self.ent_empirical_value_B.grid(row=row_num, column=2)
        row_num                           +=1
        CreateToolTip(self.lbl_empirical_value, text =
            "An empirical coefficient that varies depending on\n"
            "the site of the telescope. This is defined as C_Y in\n"
            "Osborn et al. 2015 (DOI: 10.1093/mnras/stv1400), where\n"
            "several sites are listed. The mean value is 1.5.\n"
            "Increasing this value increases the scintillation noise.")


        ## ==================== Observational and Source Parameters ====================
        ## =============================================================================


        ttk.Separator(self.frm_inputs, orient='vertical').grid(row=0, column=1, rowspan=3, sticky='nsew', padx=2)

        ## =============================================================================
        ## ============================ Fourier Parameters =============================

        ## Set up the frame for all Fourier parameters.
        over_row_num        = 0
        row_num             = 0
        self.frm_fourier    = tk.Frame(master=self.frm_inputs)
        self.frm_fourier.grid(row=over_row_num, column=2, padx=0, pady=0)
        over_row_num        +=1


        ## ========== Power Spectra ==========

        ## Set up the frame for Power Spectra.
        row_num             = 0
        mid_row_num         = 0
        self.frm_ps         = tk.LabelFrame(master=self.frm_fourier, text="Define Power Spectra", font=('bitstream vera sans', 9, 'bold'), padx=3, pady=2)
        self.frm_ps.grid(row=mid_row_num, column=0, padx=0, pady=0, sticky="nsew")
        mid_row_num         +=1


        ## ========== Power Spectral Choices ==========

        self.power_model_type            = tk.IntVar()

        ## Set up the frame for the broken powerlaw radio button
        ## (These frames are cut up in this way for layout purposes)
        row_num                   = 0
        self.frm_ps_0             = tk.Frame(master=self.frm_ps)
        self.frm_ps_0.grid(row=mid_row_num, column=0, padx=0, pady=0)

        ## Create the broken powerlaw radio button
        self.rad_bpl              = tk.Radiobutton(master=self.frm_ps_0, text="Broken Powerlaw",    variable = self.power_model_type, value = 1, command=self.psd_enable_disable)
        self.rad_bpl.grid(row = 0, column = 0, sticky='w', ipadx=5, ipady=0)
        self.rad_bpl.select()
        CreateToolTip(self.rad_bpl, text =
            "Uses a simple model of constant power that breaks at\n"
            "some frequency, and then becomes a power law component\n"
            "that decreases with increasing frequency.\n"
            "The coherence is described in the same way.")

        ## Set up the frame for the lorentzian radio button
        ## (These frames are cut up in this way for layout purposes)
        row_num                   = 0
        self.frm_ps_1             = tk.Frame(master=self.frm_ps)
        self.frm_ps_1.grid(row=mid_row_num, column=2, padx=0, pady=0)
        over_row_num        +=1

        ## Create the Lorentzian radio button
        self.rad_lor              = tk.Radiobutton(master=self.frm_ps_1, text="Lorentzians",        variable = self.power_model_type, value = 2, command=self.psd_enable_disable)
        self.rad_lor.grid(row = 0, column = 1, sticky='w', ipadx=5, ipady=0, rowspan=2)
        mid_row_num               +=1
        CreateToolTip(self.rad_lor, text =
            "Allows you to define a series of Lorentzians that\n"
            "sum to the power spectrum for each series. Series B\n"
            "includes an extra parameter for each Lorentzian\n"
            "to define how coherent it is with Series A.\n"
            "\n"
            "Lorentzians are in the form:\n"
            "L(f) = N/π (0.5*W)/((f-M)^2+(0.5*W)^2)\n"
            "Where f is frequency, and N, W, and M are the\n"
            "Normalisation, Width, and Midpoint respectively.")




        ## ========== Power Spectral Parameters ==========

        ## ----- Broken Powerlaw -----

        ## Set up the frame for the broken powerlaw values
        ## (These frames are cut up in this way for layout purposes)
        row_num                = 1
        self.frm_ps_2          = tk.Frame(master=self.frm_ps)
        self.frm_ps_2.grid(row=mid_row_num, column=0, padx=0, pady=2, sticky="ns")
        over_row_num           +=1
        mid_row_num            +=1


        self.ps_power          = tk.StringVar()
        self.break_freq        = tk.StringVar()
        self.coh_constant      = tk.StringVar()
        self.coh_power         = tk.StringVar()
        self.coh_break_freq    = tk.StringVar()

        ## Spacer
        tk.Label(master=self.frm_ps_2, text="Power Spectra...").grid(row=row_num, column=0, sticky="e")
        row_num                +=1

        ## Power Spectral Index
        self.lbl_ps_power      = tk.Label(master=self.frm_ps_2, text="Index:")
        self.lbl_ps_power.grid(row=row_num, column=0, sticky="e")
        self.ent_ps_power      = tk.Entry(master=self.frm_ps_2, width=6, textvariable=self.ps_power)
        self.ent_ps_power.grid(row=row_num, column=1)
        row_num                +=1
        CreateToolTip(self.lbl_ps_power, text =
            "Index of the power law component of the power spectra.")

        ## Power Spectral Break Frequency
        self.lbl_break_freq    = tk.Label(master=self.frm_ps_2, text="Break Freq.:")
        self.lbl_break_freq.grid(row=row_num, column=0, sticky="e")
        self.ent_break_freq    = tk.Entry(master=self.frm_ps_2, width=6, textvariable=self.break_freq)
        self.ent_break_freq.grid(row=row_num, column=1)
        row_num                +=1
        CreateToolTip(self.lbl_break_freq, text =
            "Break Frequency of the power law\n"
            "(the frequency at which it transitions\n"
            "from a constant to a power law).")

        tk.Label(master=self.frm_ps_2, text=" ").grid(row=row_num, column=0, sticky="e")
        row_num                +=1

        tk.Label(master=self.frm_ps_2, text="Coherence...").grid(row=row_num, column=0, sticky="e")
        row_num                +=1

        ## Coherence Constant
        self.lbl_coh_constant    = tk.Label(master=self.frm_ps_2, text="Constant:")
        self.lbl_coh_constant.grid(row=row_num, column=0, sticky="e")
        self.ent_coh_constant    = tk.Entry(master=self.frm_ps_2, width=6, textvariable=self.coh_constant)
        self.ent_coh_constant.grid(row=row_num, column=1)
        row_num                  +=1
        CreateToolTip(self.lbl_coh_constant, text =
            "The fraction that Series B is coherent with Series A\n"
            "during the constant component of the coherence.")

        ## Coherence Index
        self.lbl_coh_power        = tk.Label(master=self.frm_ps_2, text="Index:")
        self.lbl_coh_power.grid(row=row_num, column=0, sticky="e")
        self.ent_coh_power        = tk.Entry(master=self.frm_ps_2, width=6, textvariable=self.coh_power)
        self.ent_coh_power.grid(row=row_num, column=1)
        row_num                   +=1
        CreateToolTip(self.lbl_coh_power, text =
            "Index of the Power Law component of the coherence.")

        ## Coherence Break Frequency
        self.lbl_coh_break_freq    = tk.Label(master=self.frm_ps_2, text="Break Freq.:")
        self.lbl_coh_break_freq.grid(row=row_num, column=0, sticky="e")
        self.ent_coh_break_freq    = tk.Entry(master=self.frm_ps_2, width=6, textvariable=self.coh_break_freq)
        self.ent_coh_break_freq.grid(row=row_num, column=1)
        row_num                    +=1
        CreateToolTip(self.lbl_coh_break_freq, text =
            "Break Frequency of the coherence\n"
            "(the frequency at which it transitions\n"
            "from a constant to a power law)")


        ## ----- Lorentzians -----
        ## Separator between Broken Powerlaw and Lorentzian inputs
        ttk.Separator(self.frm_ps, orient='vertical').grid(row=0, column=1, rowspan=3, sticky='nsew', padx=5)


        ## Set up the frame for the lorentzian values
        ## (These frames are cut up in this way for layout purposes)
        mid_row_num        -=1
        row_num            = 0
        self.frm_ps_3      = tk.Frame(master=self.frm_ps)
        self.frm_ps_3.grid(row=mid_row_num, column=2, sticky="nsew")#, padx=20, pady=20)
        mid_row_num        +=1

        ## -- Series A: --
        row_num = 0

        self.lorA_1_Nrm    = tk.StringVar()
        self.lorA_2_Nrm    = tk.StringVar()
        self.lorA_3_Nrm    = tk.StringVar()
        self.lorA_4_Nrm    = tk.StringVar()
        self.lorA_5_Nrm    = tk.StringVar()
        self.lorA_6_Nrm    = tk.StringVar()
        self.lorA_1_Wdt    = tk.StringVar()
        self.lorA_2_Wdt    = tk.StringVar()
        self.lorA_3_Wdt    = tk.StringVar()
        self.lorA_4_Wdt    = tk.StringVar()
        self.lorA_5_Wdt    = tk.StringVar()
        self.lorA_6_Wdt    = tk.StringVar()
        self.lorA_1_Mid    = tk.StringVar()
        self.lorA_2_Mid    = tk.StringVar()
        self.lorA_3_Mid    = tk.StringVar()
        self.lorA_4_Mid    = tk.StringVar()
        self.lorA_5_Mid    = tk.StringVar()
        self.lorA_6_Mid    = tk.StringVar()

        col_num = 0

        ## Labels
        tk.Label(master=self.frm_ps_3, text="A:").grid(row=row_num,  column=col_num+2)
        row_num +=1

        self.lbl_normalisation_A        = tk.Label(master=self.frm_ps_3, text="Norm:")
        self.lbl_width_A                = tk.Label(master=self.frm_ps_3, text="Width:")
        self.lbl_midpoint_A             = tk.Label(master=self.frm_ps_3, text="Mid:")
        self.lbl_normalisation_A.grid(row=row_num,  column=col_num+1)
        self.lbl_width_A.grid(row=row_num,  column=col_num+2)
        self.lbl_midpoint_A.grid(row=row_num,  column=col_num+3)
        row_num                         +=1
        CreateToolTip(self.lbl_normalisation_A, text = "Normalisation of the Lorentzians")
        CreateToolTip(self.lbl_width_A, text = "Width (Full Width at Half Maximum) of the Lorentzians")
        CreateToolTip(self.lbl_midpoint_A, text = "Midpoint of the Lorentzians (0 = Zero-centered Lorentzian)")

        ## Entry Boxes
        tk.Label(master=self.frm_ps_3, text="1:").grid(row=row_num+0, column=col_num, sticky="e")
        tk.Label(master=self.frm_ps_3, text="2:").grid(row=row_num+1, column=col_num, sticky="e")
        tk.Label(master=self.frm_ps_3, text="3:").grid(row=row_num+2, column=col_num, sticky="e")
        tk.Label(master=self.frm_ps_3, text="4:").grid(row=row_num+3, column=col_num, sticky="e")
        tk.Label(master=self.frm_ps_3, text="5:").grid(row=row_num+4, column=col_num, sticky="e")
        tk.Label(master=self.frm_ps_3, text="6:").grid(row=row_num+5, column=col_num, sticky="e")

        self.ent_lorA_1_Nrm            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_1_Nrm)
        self.ent_lorA_1_Wdt            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_1_Wdt)
        self.ent_lorA_1_Mid            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_1_Mid)
        self.ent_lorA_2_Nrm            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_2_Nrm)
        self.ent_lorA_2_Wdt            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_2_Wdt)
        self.ent_lorA_2_Mid            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_2_Mid)
        self.ent_lorA_3_Nrm            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_3_Nrm)
        self.ent_lorA_3_Wdt            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_3_Wdt)
        self.ent_lorA_3_Mid            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_3_Mid)
        self.ent_lorA_4_Nrm            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_4_Nrm)
        self.ent_lorA_4_Wdt            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_4_Wdt)
        self.ent_lorA_4_Mid            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_4_Mid)
        self.ent_lorA_5_Nrm            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_5_Nrm)
        self.ent_lorA_5_Wdt            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_5_Wdt)
        self.ent_lorA_5_Mid            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_5_Mid)
        self.ent_lorA_6_Nrm            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_6_Nrm)
        self.ent_lorA_6_Wdt            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_6_Wdt)
        self.ent_lorA_6_Mid            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorA_6_Mid)
        self.ent_lorA_1_Nrm.grid(row=row_num+0, column=col_num+1)
        self.ent_lorA_2_Nrm.grid(row=row_num+1, column=col_num+1)
        self.ent_lorA_3_Nrm.grid(row=row_num+2, column=col_num+1)
        self.ent_lorA_4_Nrm.grid(row=row_num+3, column=col_num+1)
        self.ent_lorA_5_Nrm.grid(row=row_num+4, column=col_num+1)
        self.ent_lorA_6_Nrm.grid(row=row_num+5, column=col_num+1)
        self.ent_lorA_1_Wdt.grid(row=row_num+0, column=col_num+2)
        self.ent_lorA_2_Wdt.grid(row=row_num+1, column=col_num+2)
        self.ent_lorA_3_Wdt.grid(row=row_num+2, column=col_num+2)
        self.ent_lorA_4_Wdt.grid(row=row_num+3, column=col_num+2)
        self.ent_lorA_5_Wdt.grid(row=row_num+4, column=col_num+2)
        self.ent_lorA_6_Wdt.grid(row=row_num+5, column=col_num+2)
        self.ent_lorA_1_Mid.grid(row=row_num+0, column=col_num+3)
        self.ent_lorA_2_Mid.grid(row=row_num+1, column=col_num+3)
        self.ent_lorA_3_Mid.grid(row=row_num+2, column=col_num+3)
        self.ent_lorA_4_Mid.grid(row=row_num+3, column=col_num+3)
        self.ent_lorA_5_Mid.grid(row=row_num+4, column=col_num+3)
        self.ent_lorA_6_Mid.grid(row=row_num+5, column=col_num+3)


        ## -- Series B: --

        ttk.Separator(self.frm_ps_3, orient='vertical').grid(row=row_num, column=col_num+4, rowspan=6, sticky='nsew', padx=5)

        self.lorB_1_Nrm    = tk.StringVar()
        self.lorB_2_Nrm    = tk.StringVar()
        self.lorB_3_Nrm    = tk.StringVar()
        self.lorB_4_Nrm    = tk.StringVar()
        self.lorB_5_Nrm    = tk.StringVar()
        self.lorB_6_Nrm    = tk.StringVar()
        self.lorB_1_Wdt    = tk.StringVar()
        self.lorB_2_Wdt    = tk.StringVar()
        self.lorB_3_Wdt    = tk.StringVar()
        self.lorB_4_Wdt    = tk.StringVar()
        self.lorB_5_Wdt    = tk.StringVar()
        self.lorB_6_Wdt    = tk.StringVar()
        self.lorB_1_Mid    = tk.StringVar()
        self.lorB_2_Mid    = tk.StringVar()
        self.lorB_3_Mid    = tk.StringVar()
        self.lorB_4_Mid    = tk.StringVar()
        self.lorB_5_Mid    = tk.StringVar()
        self.lorB_6_Mid    = tk.StringVar()
        self.lorB_1_Coh    = tk.StringVar()
        self.lorB_2_Coh    = tk.StringVar()
        self.lorB_3_Coh    = tk.StringVar()
        self.lorB_4_Coh    = tk.StringVar()
        self.lorB_5_Coh    = tk.StringVar()
        self.lorB_6_Coh    = tk.StringVar()

        ## Labels
        row_num                    -=2
        tk.Label(master=self.frm_ps_3, text="B:").grid(row=row_num,  column=col_num+6)
        row_num                    +=1

        self.lbl_normalisation_B        = tk.Label(master=self.frm_ps_3, text="Norm:")
        self.lbl_width_B                = tk.Label(master=self.frm_ps_3, text="Width:")
        self.lbl_midpoint_B             = tk.Label(master=self.frm_ps_3, text="Mid:")
        self.lbl_coherence_B            = tk.Label(master=self.frm_ps_3, text="Coh:")
        self.lbl_normalisation_B.grid(row=row_num,    column=col_num+5)
        self.lbl_width_B.grid(row=row_num,            column=col_num+6)
        self.lbl_midpoint_B.grid(row=row_num,        column=col_num+7)
        self.lbl_coherence_B.grid(row=row_num,        column=col_num+8)
        row_num                    +=1
        CreateToolTip(self.lbl_normalisation_B, text = "Normalisation of the Lorentzians")
        CreateToolTip(self.lbl_width_B, text = "Width (Full Width at Half Maximum) of the Lorentzians")
        CreateToolTip(self.lbl_midpoint_B, text = "Midpoint of the Lorentzians (0 = Zero-centered Lorentzian)")
        CreateToolTip(self.lbl_coherence_B, text = "Fraction of the Lorentzian that is coherent with Series A (0 < Coh < 1)")

        ## Entry Boxes
        self.ent_lorB_1_Nrm            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_1_Nrm)
        self.ent_lorB_1_Wdt            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_1_Wdt)
        self.ent_lorB_1_Mid            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_1_Mid)
        self.ent_lorB_1_Coh            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_1_Coh)
        self.ent_lorB_2_Nrm            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_2_Nrm)
        self.ent_lorB_2_Wdt            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_2_Wdt)
        self.ent_lorB_2_Mid            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_2_Mid)
        self.ent_lorB_2_Coh            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_2_Coh)
        self.ent_lorB_3_Nrm            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_3_Nrm)
        self.ent_lorB_3_Wdt            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_3_Wdt)
        self.ent_lorB_3_Mid            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_3_Mid)
        self.ent_lorB_3_Coh            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_3_Coh)
        self.ent_lorB_4_Nrm            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_4_Nrm)
        self.ent_lorB_4_Wdt            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_4_Wdt)
        self.ent_lorB_4_Mid            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_4_Mid)
        self.ent_lorB_4_Coh            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_4_Coh)
        self.ent_lorB_5_Nrm            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_5_Nrm)
        self.ent_lorB_5_Wdt            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_5_Wdt)
        self.ent_lorB_5_Mid            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_5_Mid)
        self.ent_lorB_5_Coh            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_5_Coh)
        self.ent_lorB_6_Nrm            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_6_Nrm)
        self.ent_lorB_6_Wdt            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_6_Wdt)
        self.ent_lorB_6_Mid            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_6_Mid)
        self.ent_lorB_6_Coh            = tk.Entry(master=self.frm_ps_3, width=6, textvariable=self.lorB_6_Coh)
        self.ent_lorB_1_Nrm.grid(row=row_num+0, column=col_num+5)
        self.ent_lorB_2_Nrm.grid(row=row_num+1, column=col_num+5)
        self.ent_lorB_3_Nrm.grid(row=row_num+2, column=col_num+5)
        self.ent_lorB_4_Nrm.grid(row=row_num+3, column=col_num+5)
        self.ent_lorB_5_Nrm.grid(row=row_num+4, column=col_num+5)
        self.ent_lorB_6_Nrm.grid(row=row_num+5, column=col_num+5)
        self.ent_lorB_1_Wdt.grid(row=row_num+0, column=col_num+6)
        self.ent_lorB_2_Wdt.grid(row=row_num+1, column=col_num+6)
        self.ent_lorB_3_Wdt.grid(row=row_num+2, column=col_num+6)
        self.ent_lorB_4_Wdt.grid(row=row_num+3, column=col_num+6)
        self.ent_lorB_5_Wdt.grid(row=row_num+4, column=col_num+6)
        self.ent_lorB_6_Wdt.grid(row=row_num+5, column=col_num+6)
        self.ent_lorB_1_Mid.grid(row=row_num+0, column=col_num+7)
        self.ent_lorB_2_Mid.grid(row=row_num+1, column=col_num+7)
        self.ent_lorB_3_Mid.grid(row=row_num+2, column=col_num+7)
        self.ent_lorB_4_Mid.grid(row=row_num+3, column=col_num+7)
        self.ent_lorB_5_Mid.grid(row=row_num+4, column=col_num+7)
        self.ent_lorB_6_Mid.grid(row=row_num+5, column=col_num+7)
        self.ent_lorB_1_Coh.grid(row=row_num+0, column=col_num+8)
        self.ent_lorB_2_Coh.grid(row=row_num+1, column=col_num+8)
        self.ent_lorB_3_Coh.grid(row=row_num+2, column=col_num+8)
        self.ent_lorB_4_Coh.grid(row=row_num+3, column=col_num+8)
        self.ent_lorB_5_Coh.grid(row=row_num+4, column=col_num+8)
        self.ent_lorB_6_Coh.grid(row=row_num+5, column=col_num+8)



        ## ----- Test Power Spectra -----

        ## Set up the frame
        self.frm_ps_4         = tk.Frame(master=self.frm_ps)
        self.frm_ps_4.grid(row=mid_row_num, column=2, pady=0, sticky="e")

        ## Test Power Spectra Button
        self.btn_testpsds    = tk.Button(master=self.frm_ps_4, text="Test Power Spectra", command= lambda: self.test_PSDs( \
            self.output_dir.get(), self.fileprefix.get(), float(self.obs_length.get()), float(self.time_res.get()), float(self.power_model_type.get()), \
            float(self.ps_power.get()),      float(self.break_freq.get()),    float(self.coh_constant.get()),  float(self.coh_power.get()),     float(self.coh_break_freq.get()), \
            float(self.lorA_1_Nrm.get()),    float(self.lorA_2_Nrm.get()),    float(self.lorA_3_Nrm.get()),    float(self.lorA_4_Nrm.get()),    float(self.lorA_5_Nrm.get()),    float(self.lorA_6_Nrm.get()), \
            float(self.lorA_1_Wdt.get()),    float(self.lorA_2_Wdt.get()),    float(self.lorA_3_Wdt.get()),    float(self.lorA_4_Wdt.get()),    float(self.lorA_5_Wdt.get()),    float(self.lorA_6_Wdt.get()), \
            float(self.lorA_1_Mid.get()),    float(self.lorA_2_Mid.get()),    float(self.lorA_3_Mid.get()),    float(self.lorA_4_Mid.get()),    float(self.lorA_5_Mid.get()),    float(self.lorA_6_Mid.get()), \
            float(self.lorB_1_Nrm.get()),    float(self.lorB_2_Nrm.get()),    float(self.lorB_3_Nrm.get()),    float(self.lorB_4_Nrm.get()),    float(self.lorB_5_Nrm.get()),    float(self.lorB_6_Nrm.get()), \
            float(self.lorB_1_Wdt.get()),    float(self.lorB_2_Wdt.get()),    float(self.lorB_3_Wdt.get()),    float(self.lorB_4_Wdt.get()),    float(self.lorB_5_Wdt.get()),    float(self.lorB_6_Wdt.get()), \
            float(self.lorB_1_Mid.get()),    float(self.lorB_2_Mid.get()),    float(self.lorB_3_Mid.get()),    float(self.lorB_4_Mid.get()),    float(self.lorB_5_Mid.get()),    float(self.lorB_6_Mid.get()), \
            float(self.lorB_1_Coh.get()),    float(self.lorB_2_Coh.get()),    float(self.lorB_3_Coh.get()),    float(self.lorB_4_Coh.get()),    float(self.lorB_5_Coh.get()),    float(self.lorB_6_Coh.get())))
        self.btn_testpsds.grid(row=0, column=0)





        # ## ========== Phase/Time Lags ==========


        ## ----- Main -----

        ## Set up the frame
        col_num                = 0
        row_num                = 0
        self.frm_lags          = tk.LabelFrame(master=self.frm_fourier, text="Define Lags", font=('bitstream vera sans', 9, 'bold'), padx=3, pady=3)
        self.frm_lags.grid(row=mid_row_num, column=0, padx=0, pady=0, sticky="nsew")
        over_row_num           +=1
        mid_row_num            +=1

        ## Set up the frame for time_or_phase and overall_lag
        self.frm_lags1         = tk.Frame(master=self.frm_lags)
        self.frm_lags1.grid(row=0, column=0, padx=0, pady=0)

        self.time_or_phase     = tk.IntVar()

        ## Time and Phase radio buttons
        self.rad_time          = tk.Radiobutton(master=self.frm_lags1, text="Time",    variable = self.time_or_phase, value = 1)
        self.rad_time.grid(row = row_num, column = col_num, sticky='w', ipady=0)
        row_num                +=1
        self.rad_phase         = tk.Radiobutton(master=self.frm_lags1, text="Phase",    variable = self.time_or_phase, value = 2)
        self.rad_phase.grid(row = row_num, column = col_num, sticky='w', ipady=0)
        self.rad_phase.select()
        row_num                +=1
        CreateToolTip(self.rad_time, text =
            "Define the lags in the Time domain.\n"
            "All distributions will be in log-log space.")
        CreateToolTip(self.rad_phase, text =
            "Define the lags in the Frequency domain.\n"
            "All distributions will be in semi-log space.")

        ## Overall Lag
        self.overall_lag            = tk.StringVar()

        self.lbl_overall_lag        = tk.Label(master=self.frm_lags1, text="Overall Lag:")
        self.lbl_overall_lag.grid(row=row_num, column=0, sticky='e')
        row_num                     +=1
        self.ent_overall_lag        = tk.Entry(master=self.frm_lags1, width=10, textvariable=self.overall_lag)
        self.ent_overall_lag.grid(row=row_num, column=col_num)
        row_num                     = 0
        CreateToolTip(self.lbl_overall_lag, text =
            "The time/phase lag to use for all frequencies\n"
            "not otherwise defined to the right.")

        col_num                +=1

        ## Lags
        ttk.Separator(self.frm_lags1, orient='vertical').grid(row=row_num, column=col_num, rowspan=7, sticky='nsew', padx=5)

        col_num                +=1

        self.lag1_start_frq        = tk.StringVar()
        self.lag1_start_lag        = tk.StringVar()
        self.lag1_close_frq        = tk.StringVar()
        self.lag1_close_lag        = tk.StringVar()
        self.lag1_extra_frq        = tk.StringVar()
        self.lag1_extra_lag        = tk.StringVar()
        self.lag2_distribtn        = tk.StringVar()
        self.lag2_start_frq        = tk.StringVar()
        self.lag2_start_lag        = tk.StringVar()
        self.lag2_close_frq        = tk.StringVar()
        self.lag2_close_lag        = tk.StringVar()
        self.lag2_extra_frq        = tk.StringVar()
        self.lag2_extra_lag        = tk.StringVar()
        self.lag3_distribtn        = tk.StringVar()
        self.lag3_start_frq        = tk.StringVar()
        self.lag3_start_lag        = tk.StringVar()
        self.lag3_close_frq        = tk.StringVar()
        self.lag3_close_lag        = tk.StringVar()
        self.lag3_extra_frq        = tk.StringVar()
        self.lag3_extra_lag        = tk.StringVar()
        self.lag4_distribtn        = tk.StringVar()
        self.lag4_start_frq        = tk.StringVar()
        self.lag4_start_lag        = tk.StringVar()
        self.lag4_close_frq        = tk.StringVar()
        self.lag4_close_lag        = tk.StringVar()
        self.lag4_extra_frq        = tk.StringVar()
        self.lag4_extra_lag        = tk.StringVar()
        self.lag5_distribtn        = tk.StringVar()
        self.lag5_start_frq        = tk.StringVar()
        self.lag5_start_lag        = tk.StringVar()
        self.lag5_close_frq        = tk.StringVar()
        self.lag5_close_lag        = tk.StringVar()
        self.lag5_extra_frq        = tk.StringVar()
        self.lag5_extra_lag        = tk.StringVar()
        self.lag6_distribtn        = tk.StringVar()
        self.lag6_start_frq        = tk.StringVar()
        self.lag6_start_lag        = tk.StringVar()
        self.lag6_close_frq        = tk.StringVar()
        self.lag6_close_lag        = tk.StringVar()
        self.lag6_extra_frq        = tk.StringVar()
        self.lag6_extra_lag        = tk.StringVar()

        ## Distributions
        distribution_list = [
        "(Unused)",
        "Const. Time",
        "Const. Phase",
        "Linear",
        "Power",
        "Polynomial"
        ]

        self.lag1_distribtn = tk.StringVar()
        self.lag2_distribtn = tk.StringVar()
        self.lag3_distribtn = tk.StringVar()
        self.lag4_distribtn = tk.StringVar()
        self.lag5_distribtn = tk.StringVar()
        self.lag6_distribtn = tk.StringVar()
        self.lag1_distribtn.set(distribution_list[0]) # default value
        self.lag2_distribtn.set(distribution_list[0]) # default value
        self.lag3_distribtn.set(distribution_list[0]) # default value
        self.lag4_distribtn.set(distribution_list[0]) # default value
        self.lag5_distribtn.set(distribution_list[0]) # default value
        self.lag6_distribtn.set(distribution_list[0]) # default value

        ## Labels
        self.lbl_distribution = tk.Label(master=self.frm_lags1, text="Distribution:")
        self.lbl_distribution.grid(row=row_num, column=col_num+1)
        self.lbl_frq_1        = tk.Label(master=self.frm_lags1, text="Freq 1:"    )
        self.lbl_frq_1.grid(row=row_num, column=col_num+2)
        self.lbl_lag_1        = tk.Label(master=self.frm_lags1, text="Lag 1:"    )
        self.lbl_lag_1.grid(row=row_num, column=col_num+3)
        self.lbl_frq_2        = tk.Label(master=self.frm_lags1, text="Freq 2:"    )
        self.lbl_frq_2.grid(row=row_num, column=col_num+5)
        self.lbl_lag_2        = tk.Label(master=self.frm_lags1, text="Lag 2:"    )
        self.lbl_lag_2.grid(row=row_num, column=col_num+6)
        self.lbl_frq_3        = tk.Label(master=self.frm_lags1, text="Freq 3:"    )
        self.lbl_frq_3.grid(row=row_num, column=col_num+8)
        self.lbl_lag_3        = tk.Label(master=self.frm_lags1, text="Lag 3:"    )
        self.lbl_lag_3.grid(row=row_num, column=col_num+9)
        row_num               +=1

        CreateToolTip(self.lbl_distribution, text =
            "Select Distribution:\n"
            "\n"
            "(Unused):        Unused line.                        (Keep all entries as 0)\n"
            "Const. Time:\tConstant lag (L_1) in time.                    (Requires: Freq_1, Lag_1, Freq_2)\n"
            "Const. Phase:\tConstant lag (L_1) in phase.                    (Requires: Freq_1, Lag_1, Freq_2)\n"
            "Linear:        Linear distribution between (F_1, L_1) and (F_2, L_2).        (Requires: Freq_1, Lag_1, Freq_2, Lag_2)\n"
            "Power:        Exponentially increasing lag between (F_1, L_1) and (F_2, L_2).    (Requires: Freq_1, Lag_1, Freq_2, Lag_2)\n"
            "Polynomial:    Solves a second-order polynomial for the three coordinates given.    (Requires: Freq_1, Lag_1, Freq_2, Lag_2, Freq_3, Lag_3)\n"
        )

        tk.Label(master=self.frm_lags1, text="1:").grid(row=row_num+0, column=col_num+0, sticky='e')
        tk.Label(master=self.frm_lags1, text="2:").grid(row=row_num+1, column=col_num+0, sticky='e')
        tk.Label(master=self.frm_lags1, text="3:").grid(row=row_num+2, column=col_num+0, sticky='e')
        tk.Label(master=self.frm_lags1, text="4:").grid(row=row_num+3, column=col_num+0, sticky='e')
        tk.Label(master=self.frm_lags1, text="5:").grid(row=row_num+4, column=col_num+0, sticky='e')
        tk.Label(master=self.frm_lags1, text="6:").grid(row=row_num+5, column=col_num+0, sticky='e')



        ## Entry Boxes
        self.opt_lag1_distribtn = tk.OptionMenu(self.frm_lags1, self.lag1_distribtn, *distribution_list, command=self.lag_enable_disable_1)
        self.opt_lag1_distribtn.config(width=10)
        self.ent_lag1_start_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag1_start_frq)
        self.ent_lag1_start_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag1_start_lag)
        self.ent_lag1_close_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag1_close_frq)
        self.ent_lag1_close_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag1_close_lag)
        self.ent_lag1_extra_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag1_extra_frq)
        self.ent_lag1_extra_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag1_extra_lag)
        self.opt_lag2_distribtn = tk.OptionMenu(self.frm_lags1, self.lag2_distribtn, *distribution_list, command=self.lag_enable_disable_2)
        self.opt_lag2_distribtn.config(width=10)
        self.ent_lag2_start_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag2_start_frq)
        self.ent_lag2_start_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag2_start_lag)
        self.ent_lag2_close_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag2_close_frq)
        self.ent_lag2_close_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag2_close_lag)
        self.ent_lag2_extra_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag2_extra_frq)
        self.ent_lag2_extra_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag2_extra_lag)
        self.opt_lag3_distribtn = tk.OptionMenu(self.frm_lags1, self.lag3_distribtn, *distribution_list, command=self.lag_enable_disable_3)
        self.opt_lag3_distribtn.config(width=10)
        self.ent_lag3_start_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag3_start_frq)
        self.ent_lag3_start_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag3_start_lag)
        self.ent_lag3_close_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag3_close_frq)
        self.ent_lag3_close_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag3_close_lag)
        self.ent_lag3_extra_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag3_extra_frq)
        self.ent_lag3_extra_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag3_extra_lag)
        self.opt_lag4_distribtn = tk.OptionMenu(self.frm_lags1, self.lag4_distribtn, *distribution_list, command=self.lag_enable_disable_4)
        self.opt_lag4_distribtn.config(width=10)
        self.ent_lag4_start_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag4_start_frq)
        self.ent_lag4_start_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag4_start_lag)
        self.ent_lag4_close_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag4_close_frq)
        self.ent_lag4_close_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag4_close_lag)
        self.ent_lag4_extra_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag4_extra_frq)
        self.ent_lag4_extra_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag4_extra_lag)
        self.opt_lag5_distribtn = tk.OptionMenu(self.frm_lags1, self.lag5_distribtn, *distribution_list, command=self.lag_enable_disable_5)
        self.opt_lag5_distribtn.config(width=10)
        self.ent_lag5_start_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag5_start_frq)
        self.ent_lag5_start_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag5_start_lag)
        self.ent_lag5_close_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag5_close_frq)
        self.ent_lag5_close_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag5_close_lag)
        self.ent_lag5_extra_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag5_extra_frq)
        self.ent_lag5_extra_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag5_extra_lag)
        self.opt_lag6_distribtn = tk.OptionMenu(self.frm_lags1, self.lag6_distribtn, *distribution_list, command=self.lag_enable_disable_6)
        self.opt_lag6_distribtn.config(width=10)
        self.ent_lag6_start_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag6_start_frq)
        self.ent_lag6_start_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag6_start_lag)
        self.ent_lag6_close_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag6_close_frq)
        self.ent_lag6_close_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag6_close_lag)
        self.ent_lag6_extra_frq        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag6_extra_frq)
        self.ent_lag6_extra_lag        = tk.Entry(master=self.frm_lags1, width=6, textvariable=self.lag6_extra_lag)

        ## Distribution options
        self.opt_lag1_distribtn.grid(row=row_num+0, column=col_num+1)
        self.opt_lag2_distribtn.grid(row=row_num+1, column=col_num+1)
        self.opt_lag3_distribtn.grid(row=row_num+2, column=col_num+1)
        self.opt_lag4_distribtn.grid(row=row_num+3, column=col_num+1)
        self.opt_lag5_distribtn.grid(row=row_num+4, column=col_num+1)
        self.opt_lag6_distribtn.grid(row=row_num+5, column=col_num+1)

        ## Freq 1 and Lag 1
        self.ent_lag1_start_frq.grid(row=row_num+0, column=col_num+2)
        self.ent_lag2_start_frq.grid(row=row_num+1, column=col_num+2)
        self.ent_lag3_start_frq.grid(row=row_num+2, column=col_num+2)
        self.ent_lag4_start_frq.grid(row=row_num+3, column=col_num+2)
        self.ent_lag5_start_frq.grid(row=row_num+4, column=col_num+2)
        self.ent_lag6_start_frq.grid(row=row_num+5, column=col_num+2)
        self.ent_lag1_start_lag.grid(row=row_num+0, column=col_num+3)
        self.ent_lag2_start_lag.grid(row=row_num+1, column=col_num+3)
        self.ent_lag3_start_lag.grid(row=row_num+2, column=col_num+3)
        self.ent_lag4_start_lag.grid(row=row_num+3, column=col_num+3)
        self.ent_lag5_start_lag.grid(row=row_num+4, column=col_num+3)
        self.ent_lag6_start_lag.grid(row=row_num+5, column=col_num+3)

        ttk.Separator(self.frm_lags1, orient='vertical').grid(row=row_num, column=col_num+4, rowspan=6, sticky='nsew', padx=5)

        ## Freq 2 and Lag 2
        self.ent_lag1_close_frq.grid(row=row_num+0, column=col_num+5)
        self.ent_lag2_close_frq.grid(row=row_num+1, column=col_num+5)
        self.ent_lag3_close_frq.grid(row=row_num+2, column=col_num+5)
        self.ent_lag4_close_frq.grid(row=row_num+3, column=col_num+5)
        self.ent_lag5_close_frq.grid(row=row_num+4, column=col_num+5)
        self.ent_lag6_close_frq.grid(row=row_num+5, column=col_num+5)
        self.ent_lag1_close_lag.grid(row=row_num+0, column=col_num+6)
        self.ent_lag2_close_lag.grid(row=row_num+1, column=col_num+6)
        self.ent_lag3_close_lag.grid(row=row_num+2, column=col_num+6)
        self.ent_lag4_close_lag.grid(row=row_num+3, column=col_num+6)
        self.ent_lag5_close_lag.grid(row=row_num+4, column=col_num+6)
        self.ent_lag6_close_lag.grid(row=row_num+5, column=col_num+6)

        ttk.Separator(self.frm_lags1, orient='vertical').grid(row=row_num, column=col_num+7, rowspan=6, sticky='nsew', padx=5)

        ## Freq 3 and Lag 3
        self.ent_lag1_extra_frq.grid(row=row_num+0, column=col_num+8)
        self.ent_lag2_extra_frq.grid(row=row_num+1, column=col_num+8)
        self.ent_lag3_extra_frq.grid(row=row_num+2, column=col_num+8)
        self.ent_lag4_extra_frq.grid(row=row_num+3, column=col_num+8)
        self.ent_lag5_extra_frq.grid(row=row_num+4, column=col_num+8)
        self.ent_lag6_extra_frq.grid(row=row_num+5, column=col_num+8)
        self.ent_lag1_extra_lag.grid(row=row_num+0, column=col_num+9)
        self.ent_lag2_extra_lag.grid(row=row_num+1, column=col_num+9)
        self.ent_lag3_extra_lag.grid(row=row_num+2, column=col_num+9)
        self.ent_lag4_extra_lag.grid(row=row_num+3, column=col_num+9)
        self.ent_lag5_extra_lag.grid(row=row_num+4, column=col_num+9)
        self.ent_lag6_extra_lag.grid(row=row_num+5, column=col_num+9)


        ## ----- Test Lags -----

        self.frm_lags2         = tk.Frame(master=self.frm_lags)
        self.frm_lags2.grid(row=1, column=0, pady=0, sticky="e")

        self.btn_testlags      = tk.Button(master=self.frm_lags2, text="Test Lags", command= lambda: self.test_lags( \
            self.output_dir.get(), self.fileprefix.get(), float(self.obs_length.get()), float(self.time_res.get()), float(self.time_or_phase.get()), float(self.overall_lag.get()), \
            self.lag1_distribtn.get(), self.lag2_distribtn.get(), self.lag3_distribtn.get(), self.lag4_distribtn.get(), self.lag5_distribtn.get(), self.lag6_distribtn.get(), \
            float(self.lag1_start_frq.get()), float(self.lag2_start_frq.get()), float(self.lag3_start_frq.get()), float(self.lag4_start_frq.get()), float(self.lag5_start_frq.get()), float(self.lag6_start_frq.get()), \
            float(self.lag1_start_lag.get()), float(self.lag2_start_lag.get()), float(self.lag3_start_lag.get()), float(self.lag4_start_lag.get()), float(self.lag5_start_lag.get()), float(self.lag6_start_lag.get()), \
            float(self.lag1_close_frq.get()), float(self.lag2_close_frq.get()), float(self.lag3_close_frq.get()), float(self.lag4_close_frq.get()), float(self.lag5_close_frq.get()), float(self.lag6_close_frq.get()), \
            float(self.lag1_close_lag.get()), float(self.lag2_close_lag.get()), float(self.lag3_close_lag.get()), float(self.lag4_close_lag.get()), float(self.lag5_close_lag.get()), float(self.lag6_close_lag.get()), \
            float(self.lag1_extra_frq.get()), float(self.lag2_extra_frq.get()), float(self.lag3_extra_frq.get()), float(self.lag4_extra_frq.get()), float(self.lag5_extra_frq.get()), float(self.lag6_extra_frq.get()), \
            float(self.lag1_extra_lag.get()), float(self.lag2_extra_lag.get()), float(self.lag3_extra_lag.get()), float(self.lag4_extra_lag.get()), float(self.lag5_extra_lag.get()), float(self.lag6_extra_lag.get()), ))
        self.btn_testlags.grid(row=0, column=0)


        ## ============================ Fourier Parameters =============================
        ## =============================================================================


        ## =============================================================================
        ## ================================= Plotting ==================================

        # ttk.Separator(self.frm_inputs_left, orient='horizontal').grid(row=1, column=0, columnspan=3, sticky='nsew', pady=5)

        ## Set up the frame for plotting parameters.
        row_num                   = 0
        self.frm_plotting         = tk.LabelFrame(master=self.frm_inputs_left, text="Plot Options", font=('bitstream vera sans', 9, 'bold'), padx=3, pady=0)
        self.frm_plotting.grid(row=2, column=0, padx=0, pady=0)

        ## ----- Red Noise Plotting -----

        self.plot_rednoise        = tk.IntVar()
        self.lbl_plot_rednoise    = tk.Label(master=self.frm_plotting, text="Plot Red Noise")
        self.chk_plot_rednoise    = tk.Checkbutton(master=self.frm_plotting, variable=self.plot_rednoise)

        self.lbl_plot_rednoise.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.chk_plot_rednoise.grid(row=row_num, column=1)
        row_num                   += 1
        CreateToolTip(self.lbl_plot_rednoise, text =
            "Plot the input red noise.")

        ttk.Separator(self.frm_plotting, orient='horizontal').grid(row=row_num, column=0, columnspan=2, sticky='ew', pady=1)
        row_num                   += 1

        ## ----- CCF Plotting -----

        self.plot_model_ccf       = tk.IntVar()
        self.calculate_ccf        = tk.IntVar()
        self.ccf_segment          = tk.StringVar()
        self.ccf_binning          = tk.StringVar()

        self.chk_plot_model_ccf   = tk.Checkbutton(master=self.frm_plotting, variable=self.plot_model_ccf)
        self.chk_calculate_ccf    = tk.Checkbutton(master=self.frm_plotting, variable=self.calculate_ccf, command=self.ccf_enable_disable)

        ## Plot Model CCF
        self.lbl_plot_model_ccf   = tk.Label(master=self.frm_plotting, text="Plot Model CCF:")
        self.lbl_plot_model_ccf.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.chk_plot_model_ccf.grid(row=row_num, column=1)
        row_num                   += 1
        CreateToolTip(self.lbl_plot_model_ccf, text =
            "Plot the model Cross-Correlation Function.")

        ## Calculate CCF
        self.lbl_calculate_ccf    = tk.Label(master=self.frm_plotting, text="Calculate CCF:")
        self.lbl_calculate_ccf.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.chk_calculate_ccf.grid(row=row_num, column=1)
        row_num                   += 1
        CreateToolTip(self.lbl_calculate_ccf, text =
            "Calculate the Cross-Correlation Function\n"
            "from the simulated lightcurves.")

        ## Segment Size (s)
        self.lbl_ccf_segment        = tk.Label(master=self.frm_plotting, text="Segment Size (s):")
        self.lbl_ccf_segment.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_ccf_segment        = tk.Entry(master=self.frm_plotting, width=4, textvariable=self.ccf_segment)
        self.ent_ccf_segment.grid(row=row_num, column=1)
        row_num                     += 1
        CreateToolTip(self.lbl_ccf_segment, text =
            "Define the length (in seconds) the lightcurves are cut up\n"
            "into before the Cross-Correlation Function is calculated")

        ## Binning (Number of Bins)
        self.lbl_ccf_binning        = tk.Label(master=self.frm_plotting, text="Binning (Bins):")
        self.lbl_ccf_binning.grid(row=row_num, column=0, sticky="e", ipadx=10)
        self.ent_ccf_binning        = tk.Entry(master=self.frm_plotting, width=4, textvariable=self.ccf_binning)
        self.ent_ccf_binning.grid(row=row_num, column=1)
        row_num                     += 1
        CreateToolTip(self.lbl_ccf_binning, text =
            "Optional; if >1, average the lightcurves over this many bins\n"
            "before the Cross-Correlation Function is calculated.\n"
            "This can be done to speed up calculations.")

        ttk.Separator(self.frm_plotting, orient='vertical').grid(row=0, column=2, rowspan=6, sticky='nsew', padx=5)
        row_num += 1


        ## ----- Fourier Plotting -----

        self.plot_model_fourier      = tk.IntVar()
        self.calculate_fourier       = tk.IntVar()
        self.fourier_bins            = tk.StringVar()
        self.fourier_rebins          = tk.StringVar()
        self.reference_freq          = tk.StringVar()

        self.chk_plot_model_fourier  = tk.Checkbutton(master=self.frm_plotting, variable=self.plot_model_fourier)
        self.chk_calculate_fourier   = tk.Checkbutton(master=self.frm_plotting, variable=self.calculate_fourier, command=self.fourier_enable_disable)

        row_num                      = 0

        ## Plot Fourier Models
        self.lbl_plot_model_fourier  = tk.Label(master=self.frm_plotting, text="Plot Fourier Models:")
        self.lbl_plot_model_fourier.grid(row=row_num, column=3, sticky="e", ipadx=10)
        self.chk_plot_model_fourier.grid(row=row_num, column=4)
        row_num                      += 1
        CreateToolTip(self.lbl_plot_model_fourier, text =
            "Plot the input models of the Power Spectra,\n"
            "Coherence, Phase Lags, and Time Lags.")

        ## Calculate Fourier Products (Such as Power Spectra, Coherence, Phase Lags, and Time Lags)
        self.lbl_calculate_fourier   = tk.Label(master=self.frm_plotting, text="Calculate Fourier:")
        self.lbl_calculate_fourier.grid(row=row_num, column=3, sticky="e", ipadx=10)
        self.chk_calculate_fourier.grid(row=row_num, column=4)
        row_num                      += 1
        CreateToolTip(self.lbl_calculate_fourier, text =
            "Calculate the Power Spectra, Coherence, Phase Lags,\n"
            "and Time Lags from the simulated lightcurves.\n"
            "Uses the Stingray analysis software.")

        ## Number of Bins for Fourier calculations
        self.lbl_fourier_bins       = tk.Label(master=self.frm_plotting, text="Segment Size (bins):")
        self.lbl_fourier_bins.grid(row=row_num, column=3, sticky="e", ipadx=10)
        self.ent_fourier_bins       = tk.Entry(master=self.frm_plotting, width=4, textvariable=self.fourier_bins)
        self.ent_fourier_bins.grid(row=row_num, column=4)
        row_num                     += 1
        CreateToolTip(self.lbl_fourier_bins, text =
            "Define the length (in bins) the lightcurves are cut up\n"
            "into before the Fourier analysis is carried out.\n"
            "Ideally, this should be a power of two.")

        ## Amount of Rebinning for Fourier Calulations
        self.lbl_fourier_rebins        = tk.Label(master=self.frm_plotting, text="Rebinning:")
        self.lbl_fourier_rebins.grid(row=row_num, column=3, sticky="e", ipadx=10)
        self.ent_fourier_rebins        = tk.Entry(master=self.frm_plotting, width=4, textvariable=self.fourier_rebins)
        self.ent_fourier_rebins.grid(row=row_num, column=4)
        row_num                        += 1
        CreateToolTip(self.lbl_fourier_rebins, text =
            ">0. Define the amount of logarithmic rebinning to do\n"
            "when calculating and plotting the Fourier analysis.")

        ## Reference Frequency for Time Lags
        self.lbl_reference_freq        = tk.Label(master=self.frm_plotting, text="Reference Freq. (Hz):")
        self.lbl_reference_freq.grid(row=row_num, column=3, sticky="e", ipadx=10)
        self.ent_reference_freq        = tk.Entry(master=self.frm_plotting, width=4, textvariable=self.reference_freq)
        self.ent_reference_freq.grid(row=row_num, column=4)
        row_num                        += 1
        CreateToolTip(self.lbl_reference_freq, text =
            "This is used when calculating the phase lags. Phase lags are\n"
            "inherently constrained between +/-pi. If the actual phase lags\n"
            "are outside of this range, e.g. between pi and 2pi, then analysis\n"
            "will show that they are between -pi and 0 radians due to the\n"
            "cyclical nature of sine waves. A ‘reference frequency’ is thus\n"
            "defined to be the frequency at which the code calibrates the rest\n"
            "of the phase lags; this should be the frequency at which the\n"
            "measured phase lag is between +/-pi.")


        ## ----- White Noise -----
        ## Whether or not to (very simply) calculate and remove the white noise
        self.remove_whitenoise        = tk.IntVar()

        self.lbl_remove_whitenoise    = tk.Label(master=self.frm_plotting, text="Remove White Noise:")
        self.lbl_remove_whitenoise.grid(row=row_num, column=3, sticky="e")
        self.chk_remove_whitenoise    = tk.Checkbutton(master=self.frm_plotting, variable=self.remove_whitenoise)
        self.chk_remove_whitenoise.grid(row=row_num, column=4)
        row_num                       += 1
        CreateToolTip(self.lbl_remove_whitenoise, text =
            "This is used for calculating the power spectra and coherence.\n"
            "It attempts to remove any white noise by assuming that\n"
            "the last points in both power spectra are white-noise dominated,\n"
            "and removes 99% of those value from the power spectra.\n"
            "\n"
            "Note: This is ONLY a good approximation when the spectrum\n"
            "is white-noise dominated at the highest Fourier frequency.")



        ## ================================= Plotting ==================================
        ## =============================================================================


        ## =============================================================================
        ## ========================= Names and Code Execution ==========================

        ## Set up the frame.
        row_num                    = 0
        self.frm_final             = tk.LabelFrame(master=self.frm_fourier, text="Final", font=('bitstream vera sans', 9, 'bold'), padx=3, pady=5)
        self.frm_final.grid(row=mid_row_num, column=0, sticky="nsew")

        # frm_names             = tk.LabelFrame(master=frm_final, text="Output Directory and Filenames", font=('bitstream vera sans', 9, 'bold'), padx=10, pady=10)
        self.frm_names             = tk.Frame(master=self.frm_final)
        self.frm_names.grid(row=0, column=0, padx=20)

        ## Row 1: Output Directory
        self.output_dir            = tk.StringVar()
        self.lbl_output_dir        = tk.Label(master=self.frm_names, text="Output Directory:")
        self.lbl_output_dir.grid(row=row_num, column=0, sticky="e")
        self.ent_output_dir        = tk.Entry(master=self.frm_names, width=10, textvariable=self.output_dir)
        self.ent_output_dir.grid(row=row_num, column=1, ipadx=50)
        row_num                    += 1
        CreateToolTip(self.lbl_output_dir, text =
            "Defines the output directory; e.g. ‘CorrSim_Outputs’ will\n"
            "save all outputs to the folder ./CorrSim_Outputs.\n"
            "This folder will be created if it doesn't already exist.")

        ## Row 2: Fileprefix
        self.fileprefix            = tk.StringVar()
        self.lbl_fileprefix        = tk.Label(master=self.frm_names, text="File prefix:")
        self.lbl_fileprefix.grid(row=row_num, column=0, sticky="e")
        self.ent_fileprefix        = tk.Entry(master=self.frm_names, width=10, textvariable=self.fileprefix)
        self.ent_fileprefix.grid(row=row_num, column=1, ipadx=50)
        row_num                    += 1
        CreateToolTip(self.lbl_fileprefix, text =
            "Defines a prefix for the names of all outputted files.")



        ## ----- Execute CorrSim -----


        row_num                  = 0
        # frm_execute         = tk.LabelFrame(master=frm_final, text="Execute CorrSim", font=('bitstream vera sans', 9, 'bold'), padx=10, pady=10)
        self.frm_execute         = tk.Frame(master=self.frm_final)
        self.frm_execute.grid(row=0, column=1, padx=10)
        row_num                  += 1

        self.btn_corrsim    = tk.Button(master=self.frm_execute, text="Execute CorrSim", command= lambda: self.execute_CorrSim( \
            self.output_dir.get(), self.fileprefix.get(), float(self.obs_length.get()), float(self.time_res.get()), \
            float(self.mean_counts_A.get()), float(self.mean_counts_B.get()), float(self.F_rms_A.get()), float(self.F_rms_B.get()),
            self.apply_red_A.get(), self.apply_red_B.get(), \
            float(self.F_rms_rednoise_A.get()), float(self.F_rms_rednoise_B.get()), float(self.rednoise_slope_A.get()), float(self.rednoise_slope_B.get()), \
            self.apply_poisson_A.get(), self.apply_poisson_B.get(), \
            self.apply_readout_A.get(), self.apply_readout_B.get(), float(self.read_noise_A.get()), float(self.read_noise_B.get()), \
            self.apply_scintillation_A.get(), self.apply_scintillation_B.get(), \
            float(self.empirical_value_A.get()), float(self.telescope_diameter_A.get()), float(self.exposure_time_A.get()), \
            float(self.target_altitude_A.get()), float(self.telescope_altitude_A.get()), float(self.turbulence_height_A.get()), \
            float(self.empirical_value_B.get()), float(self.telescope_diameter_B.get()), float(self.exposure_time_B.get()), \
            float(self.target_altitude_B.get()), float(self.telescope_altitude_B.get()), float(self.turbulence_height_B.get()), \
            float(self.power_model_type.get()),
            float(self.ps_power.get()), float(self.break_freq.get()), float(self.coh_constant.get()), float(self.coh_power.get()), float(self.coh_break_freq.get()), \
            float(self.lorA_1_Nrm.get()), float(self.lorA_2_Nrm.get()), float(self.lorA_3_Nrm.get()), float(self.lorA_4_Nrm.get()), float(self.lorA_5_Nrm.get()), float(self.lorA_6_Nrm.get()), \
            float(self.lorA_1_Wdt.get()), float(self.lorA_2_Wdt.get()), float(self.lorA_3_Wdt.get()), float(self.lorA_4_Wdt.get()), float(self.lorA_5_Wdt.get()), float(self.lorA_6_Wdt.get()), \
            float(self.lorA_1_Mid.get()), float(self.lorA_2_Mid.get()), float(self.lorA_3_Mid.get()), float(self.lorA_4_Mid.get()), float(self.lorA_5_Mid.get()), float(self.lorA_6_Mid.get()), \
            float(self.lorB_1_Nrm.get()), float(self.lorB_2_Nrm.get()), float(self.lorB_3_Nrm.get()), float(self.lorB_4_Nrm.get()), float(self.lorB_5_Nrm.get()), float(self.lorB_6_Nrm.get()), \
            float(self.lorB_1_Wdt.get()), float(self.lorB_2_Wdt.get()), float(self.lorB_3_Wdt.get()), float(self.lorB_4_Wdt.get()), float(self.lorB_5_Wdt.get()), float(self.lorB_6_Wdt.get()), \
            float(self.lorB_1_Mid.get()), float(self.lorB_2_Mid.get()), float(self.lorB_3_Mid.get()), float(self.lorB_4_Mid.get()), float(self.lorB_5_Mid.get()), float(self.lorB_6_Mid.get()), \
            float(self.lorB_1_Coh.get()), float(self.lorB_2_Coh.get()), float(self.lorB_3_Coh.get()), float(self.lorB_4_Coh.get()), float(self.lorB_5_Coh.get()), float(self.lorB_6_Coh.get()), \
            float(self.time_or_phase.get()), float(self.overall_lag.get()), \
            self.lag1_distribtn.get(), self.lag2_distribtn.get(), self.lag3_distribtn.get(), self.lag4_distribtn.get(), self.lag5_distribtn.get(), self.lag6_distribtn.get(), \
            float(self.lag1_start_frq.get()), float(self.lag2_start_frq.get()), float(self.lag3_start_frq.get()), float(self.lag4_start_frq.get()), float(self.lag5_start_frq.get()), float(self.lag6_start_frq.get()), \
            float(self.lag1_start_lag.get()), float(self.lag2_start_lag.get()), float(self.lag3_start_lag.get()), float(self.lag4_start_lag.get()), float(self.lag5_start_lag.get()), float(self.lag6_start_lag.get()), \
            float(self.lag1_close_frq.get()), float(self.lag2_close_frq.get()), float(self.lag3_close_frq.get()), float(self.lag4_close_frq.get()), float(self.lag5_close_frq.get()), float(self.lag6_close_frq.get()), \
            float(self.lag1_close_lag.get()), float(self.lag2_close_lag.get()), float(self.lag3_close_lag.get()), float(self.lag4_close_lag.get()), float(self.lag5_close_lag.get()), float(self.lag6_close_lag.get()), \
            float(self.lag1_extra_frq.get()), float(self.lag2_extra_frq.get()), float(self.lag3_extra_frq.get()), float(self.lag4_extra_frq.get()), float(self.lag5_extra_frq.get()), float(self.lag6_extra_frq.get()), \
            float(self.lag1_extra_lag.get()), float(self.lag2_extra_lag.get()), float(self.lag3_extra_lag.get()), float(self.lag4_extra_lag.get()), float(self.lag5_extra_lag.get()), float(self.lag6_extra_lag.get()), \
            self.plot_rednoise.get(), self.plot_model_ccf.get(), self.calculate_ccf.get(), float(self.ccf_segment.get()), float(self.ccf_binning.get()), \
            self.plot_model_fourier.get(), self.calculate_fourier.get(), \
            float(self.fourier_bins.get()), float(self.fourier_rebins.get()), float(self.reference_freq.get()), self.remove_whitenoise.get() \
            ))
        self.btn_corrsim.grid(row=row_num, column=1, padx=35, pady=5)


        ## ========================= Names and Code Execution ==========================
        ## =============================================================================




        ## =============================================================================
        ## ============================ Set Default Values =============================

        self.ent_output_dir.insert(0, "CorrSim_Outputs")
        self.ent_fileprefix.insert(0, "CorrSim - TestA - ")
        self.ent_obs_length.insert(0, 1000)             ## Length of observation (seconds)
        self.ent_time_res.insert(0, 0.1)                ## Time resolution of observation (s)

        self.ent_mean_counts_A.insert(0, 1000)          ## Mean Count Rate in A (Counts/s)
        self.ent_mean_counts_B.insert(0, 5000)          ## Mean Count Rate in B (Counts/s)
        self.ent_F_rms_A.insert(0, 0.3)                 ## Fractional RMS in A (For XRBs in X-rays: Hard state ~0.3, Soft State ~0.04)
        self.ent_F_rms_B.insert(0, 0.1)                 ## Fractional RMS in B (For XRBs in Optical: Hard state ~0.1)

        self.chk_apply_red_A.select()                   ## Apply red noise to Series A
        self.chk_apply_red_B.select()                   ## Apply red noise to Series B
        self.ent_red_frms_A.insert(0, 0.2)              ## Readout Noise in Electrons - HiPERCAM
        self.ent_red_frms_B.insert(0, 0.03)             ## Readout Noise in Electrons - HiPERCAM
        self.ent_red_slope_A.insert(0, -2)              ## Readout Noise in Electrons - HiPERCAM
        self.ent_red_slope_B.insert(0, -2)              ## Readout Noise in Electrons - HiPERCAM

        self.chk_apply_poisson_A.select()               ## Apply poisson noise to Series A
        self.chk_apply_poisson_B.select()               ## Apply poisson noise to Series B

        self.chk_apply_readout_A.select()               ## Apply readout noise to Series A
        self.ent_read_noise_A.insert(0, 4.5)            ## Readout Noise in Electrons - HiPERCAM
        self.ent_read_noise_B.insert(0, 4.5)            ## Readout Noise in Electrons - HiPERCAM

        self.chk_apply_scintillation_A.select()         ## Apply scintillation noise to Series A
        self.ent_empirical_value_A.insert(0, 1.5)       ## Empirical Value from each telescope (Present: Mean)
        self.ent_telescope_diameter_A.insert(0, 3.58)   ## Telescope Diameter (m) (Present: NTT, La Silla)
        self.ent_exposure_time_A.insert(0, 0.09)        ## Exposure time (s)
        self.ent_target_altitude_A.insert(0, 40)        ## Altitude of Source (degrees)
        self.ent_telescope_altitude_A.insert(0, 2375)   ## Altitude of Observatory - (Present: La Silla)
        self.ent_turbulence_height_A.insert(0, 8000)    ## Scale height of Atmospheric Turbulence
        self.ent_empirical_value_B.insert(0, 1.5)       ## Empirical Value from each telescope (Present: Mean)
        self.ent_telescope_diameter_B.insert(0, 3.58)   ## Telescope Diameter (m) (Present: NTT, La Silla)
        self.ent_exposure_time_B.insert(0, 0.09)        ## Exposure time (s)
        self.ent_target_altitude_B.insert(0, 40)        ## Altitude of Source (degrees)
        self.ent_telescope_altitude_B.insert(0, 2375)   ## Altitude of Observatory - (Present: La Silla)
        self.ent_turbulence_height_B.insert(0, 8000)    ## Scale height of Atmospheric Turbulence


        self.ent_ps_power.insert(0,-1.5)                ## Power Law component
        self.ent_break_freq.insert(0,1)                 ## Break Frequency of the power law
        self.ent_coh_constant.insert(0,0.1)             ## Constant Coherence
        self.ent_coh_power.insert(0,-0.5)               ## Power Law component for Coherence
        self.ent_coh_break_freq.insert(0,0.5)           ## Break Frequency of the power law for the Coherence

        ## Lorentzian Values
        self.ent_lorA_1_Nrm.insert(0,0)
        self.ent_lorA_2_Nrm.insert(0,0)
        self.ent_lorA_3_Nrm.insert(0,0)
        self.ent_lorA_4_Nrm.insert(0,0)
        self.ent_lorA_5_Nrm.insert(0,0)
        self.ent_lorA_6_Nrm.insert(0,0)
        self.ent_lorA_1_Wdt.insert(0,0)
        self.ent_lorA_2_Wdt.insert(0,0)
        self.ent_lorA_3_Wdt.insert(0,0)
        self.ent_lorA_4_Wdt.insert(0,0)
        self.ent_lorA_5_Wdt.insert(0,0)
        self.ent_lorA_6_Wdt.insert(0,0)
        self.ent_lorA_1_Mid.insert(0,0)
        self.ent_lorA_2_Mid.insert(0,0)
        self.ent_lorA_3_Mid.insert(0,0)
        self.ent_lorA_4_Mid.insert(0,0)
        self.ent_lorA_5_Mid.insert(0,0)
        self.ent_lorA_6_Mid.insert(0,0)
        self.ent_lorB_1_Nrm.insert(0,0)
        self.ent_lorB_2_Nrm.insert(0,0)
        self.ent_lorB_3_Nrm.insert(0,0)
        self.ent_lorB_4_Nrm.insert(0,0)
        self.ent_lorB_5_Nrm.insert(0,0)
        self.ent_lorB_6_Nrm.insert(0,0)
        self.ent_lorB_1_Wdt.insert(0,0)
        self.ent_lorB_2_Wdt.insert(0,0)
        self.ent_lorB_3_Wdt.insert(0,0)
        self.ent_lorB_4_Wdt.insert(0,0)
        self.ent_lorB_5_Wdt.insert(0,0)
        self.ent_lorB_6_Wdt.insert(0,0)
        self.ent_lorB_1_Mid.insert(0,0)
        self.ent_lorB_2_Mid.insert(0,0)
        self.ent_lorB_3_Mid.insert(0,0)
        self.ent_lorB_4_Mid.insert(0,0)
        self.ent_lorB_5_Mid.insert(0,0)
        self.ent_lorB_6_Mid.insert(0,0)
        self.ent_lorB_1_Coh.insert(0,0)
        self.ent_lorB_2_Coh.insert(0,0)
        self.ent_lorB_3_Coh.insert(0,0)
        self.ent_lorB_4_Coh.insert(0,0)
        self.ent_lorB_5_Coh.insert(0,0)
        self.ent_lorB_6_Coh.insert(0,0)

        ## Phase/Time Lag Values
        self.ent_overall_lag.insert(0,0)                    ## Overall Lag

        self.ent_lag1_start_frq.insert(0,0)
        self.ent_lag2_start_frq.insert(0,0)
        self.ent_lag3_start_frq.insert(0,0)
        self.ent_lag4_start_frq.insert(0,0)
        self.ent_lag5_start_frq.insert(0,0)
        self.ent_lag6_start_frq.insert(0,0)
        self.ent_lag1_start_lag.insert(0,0)
        self.ent_lag2_start_lag.insert(0,0)
        self.ent_lag3_start_lag.insert(0,0)
        self.ent_lag4_start_lag.insert(0,0)
        self.ent_lag5_start_lag.insert(0,0)
        self.ent_lag6_start_lag.insert(0,0)
        self.ent_lag1_close_frq.insert(0,0)
        self.ent_lag2_close_frq.insert(0,0)
        self.ent_lag3_close_frq.insert(0,0)
        self.ent_lag4_close_frq.insert(0,0)
        self.ent_lag5_close_frq.insert(0,0)
        self.ent_lag6_close_frq.insert(0,0)
        self.ent_lag1_close_lag.insert(0,0)
        self.ent_lag2_close_lag.insert(0,0)
        self.ent_lag3_close_lag.insert(0,0)
        self.ent_lag4_close_lag.insert(0,0)
        self.ent_lag5_close_lag.insert(0,0)
        self.ent_lag6_close_lag.insert(0,0)
        self.ent_lag1_extra_frq.insert(0,0)
        self.ent_lag2_extra_frq.insert(0,0)
        self.ent_lag3_extra_frq.insert(0,0)
        self.ent_lag4_extra_frq.insert(0,0)
        self.ent_lag5_extra_frq.insert(0,0)
        self.ent_lag6_extra_frq.insert(0,0)
        self.ent_lag1_extra_lag.insert(0,0)
        self.ent_lag2_extra_lag.insert(0,0)
        self.ent_lag3_extra_lag.insert(0,0)
        self.ent_lag4_extra_lag.insert(0,0)
        self.ent_lag5_extra_lag.insert(0,0)
        self.ent_lag6_extra_lag.insert(0,0)

        self.ent_ccf_segment.insert(0,30)                ## Segment size of averaged CCF, in seconds.
        self.ent_ccf_binning.insert(0,0)                 ## Optional: bins the data (0 or 1 = no binning)

        self.ent_fourier_bins.insert(0,512)              ## Number of bins per segment in Fourier analysis (Ideally a power of 2)
        self.ent_fourier_rebins.insert(0,0.2)            ## Amount of logarithmic rebinning
        self.ent_reference_freq.insert(0,1)              ## Frequency at which the phase lag can be assumed to be correct (Not shifted +/-2pi)

        ## Disable all noise fields until used
        self.read_enable_disable_A()
        self.read_enable_disable_B()
        self.scin_enable_disable_A()
        self.scin_enable_disable_B()

        ## Disable Power Spectral fields until used
        self.psd_enable_disable()

        ## Disable all lag fields until used
        self.lag_enable_disable_1("(Unused)")
        self.lag_enable_disable_2("(Unused)")
        self.lag_enable_disable_3("(Unused)")
        self.lag_enable_disable_4("(Unused)")
        self.lag_enable_disable_5("(Unused)")
        self.lag_enable_disable_6("(Unused)")

        self.ccf_enable_disable()
        self.fourier_enable_disable()

        ## Add Menu bar
        self.config(menu=self.menubar)

        ## When closing, properly shutdown
        self.protocol("WM_DELETE_WINDOW", self.shutdown_ttk_repeat)


    ### =========================================
    ### --------------- Functions ---------------

    ## ----- Enable/Disable entry fields -----
    def red_enable_disable_A(self):
        self.ent_red_frms_A.configure(state='disable')
        self.ent_red_slope_A.configure(state='disable')
        if self.apply_red_B.get() == 0:
            self.chk_plot_rednoise.configure(state='disable')
        if self.apply_red_A.get() == 1:
            self.ent_red_frms_A.configure(state='normal')
            self.ent_red_slope_A.configure(state='normal')
            self.chk_plot_rednoise.configure(state='normal')

    def red_enable_disable_B(self):
        self.ent_red_frms_B.configure(state='disable')
        self.ent_red_slope_B.configure(state='disable')
        if self.apply_red_A.get() == 0:
            self.chk_plot_rednoise.configure(state='disable')
        if self.apply_red_B.get() == 1:
            self.ent_red_frms_B.configure(state='normal')
            self.ent_red_slope_B.configure(state='normal')
            self.chk_plot_rednoise.configure(state='normal')

    def read_enable_disable_A(self):
        self.ent_read_noise_A.configure(state='disable')
        if self.apply_readout_A.get() == 1:
            self.ent_read_noise_A.configure(state='normal')

    def read_enable_disable_B(self):
        self.ent_read_noise_B.configure(state='disable')
        if self.apply_readout_B.get() == 1:
            self.ent_read_noise_B.configure(state='normal')

    def scin_enable_disable_A(self):
        scin_list_A = [self.ent_telescope_diameter_A, self.ent_telescope_altitude_A, self.ent_exposure_time_A, self.ent_turbulence_height_A, self.ent_target_altitude_A, self.ent_empirical_value_A]
        [x.configure(state='disable') for x in scin_list_A]
        if self.apply_scintillation_A.get() == 1:
            [x.configure(state='normal') for x in scin_list_A]

    def scin_enable_disable_B(self):
        scin_list_B = [self.ent_telescope_diameter_B, self.ent_telescope_altitude_B, self.ent_exposure_time_B, self.ent_turbulence_height_B, self.ent_target_altitude_B, self.ent_empirical_value_B]
        [x.configure(state='disable') for x in scin_list_B]
        if self.apply_scintillation_B.get() == 1:
            [x.configure(state='normal') for x in scin_list_B]

    def ccf_enable_disable(self):
        ccf_list = [self.ent_ccf_segment, self.ent_ccf_binning]
        [x.configure(state='disable') for x in ccf_list]
        if self.calculate_ccf.get() == 1:
            [x.configure(state='normal') for x in ccf_list]

    def fourier_enable_disable(self):
        fourier_list = [self.ent_fourier_bins, self.ent_fourier_rebins, self.ent_reference_freq]
        [x.configure(state='disable') for x in fourier_list]
        self.chk_remove_whitenoise.configure(state='disable')
        if self.calculate_fourier.get() == 1:
            [x.configure(state='normal') for x in fourier_list]
            self.chk_remove_whitenoise.configure(state='normal')


    def psd_enable_disable(self):
        bpl_list = [self.ent_ps_power, self.ent_break_freq, self.ent_coh_constant, self.ent_coh_power, self.ent_coh_break_freq]
        lor_list = [
            self.ent_lorA_1_Nrm, self.ent_lorA_1_Wdt, self.ent_lorA_1_Mid, self.ent_lorA_2_Nrm, self.ent_lorA_2_Wdt, self.ent_lorA_2_Mid,
            self.ent_lorA_3_Nrm, self.ent_lorA_3_Wdt, self.ent_lorA_3_Mid, self.ent_lorA_4_Nrm, self.ent_lorA_4_Wdt, self.ent_lorA_4_Mid,
            self.ent_lorA_5_Nrm, self.ent_lorA_5_Wdt, self.ent_lorA_5_Mid, self.ent_lorA_6_Nrm, self.ent_lorA_6_Wdt, self.ent_lorA_6_Mid,
            self.ent_lorB_1_Nrm, self.ent_lorB_1_Wdt, self.ent_lorB_1_Mid, self.ent_lorB_1_Coh, self.ent_lorB_2_Nrm, self.ent_lorB_2_Wdt, self.ent_lorB_2_Mid, self.ent_lorB_2_Coh,
            self.ent_lorB_3_Nrm, self.ent_lorB_3_Wdt, self.ent_lorB_3_Mid, self.ent_lorB_3_Coh, self.ent_lorB_4_Nrm, self.ent_lorB_4_Wdt, self.ent_lorB_4_Mid, self.ent_lorB_4_Coh,
            self.ent_lorB_5_Nrm, self.ent_lorB_5_Wdt, self.ent_lorB_5_Mid, self.ent_lorB_5_Coh, self.ent_lorB_6_Nrm, self.ent_lorB_6_Wdt, self.ent_lorB_6_Mid, self.ent_lorB_6_Coh]
        [x.configure(state='disable') for x in bpl_list]
        [x.configure(state='disable') for x in lor_list]
        if self.power_model_type.get() == 1:
            [x.configure(state='normal') for x in bpl_list]
        elif self.power_model_type.get() == 2:
            [x.configure(state='normal') for x in lor_list]

    def lag_enable_disable_1(self, distribution):
        ent_lag1_list = [self.ent_lag1_start_frq, self.ent_lag1_start_lag, self.ent_lag1_close_frq, self.ent_lag1_close_lag, self.ent_lag1_extra_frq, self.ent_lag1_extra_lag]
        [x.configure(state='disable') for x in ent_lag1_list]
        if distribution != "(Unused)":
            self.ent_lag1_start_frq.configure(state='normal')
            self.ent_lag1_start_lag.configure(state='normal')
            self.ent_lag1_close_frq.configure(state='normal')
        if distribution == "Linear" or distribution == "Power" or distribution == "Polynomial":
            self.ent_lag1_close_lag.configure(state='normal')
        if distribution == "Polynomial":
            self.ent_lag1_extra_frq.configure(state='normal')
            self.ent_lag1_extra_lag.configure(state='normal')

    def lag_enable_disable_2(self, distribution):
        ent_lag2_list = [self.ent_lag2_start_frq, self.ent_lag2_start_lag, self.ent_lag2_close_frq, self.ent_lag2_close_lag, self.ent_lag2_extra_frq, self.ent_lag2_extra_lag]
        [x.configure(state='disable') for x in ent_lag2_list]
        if distribution != "(Unused)":
            self.ent_lag2_start_frq.configure(state='normal')
            self.ent_lag2_start_lag.configure(state='normal')
            self.ent_lag2_close_frq.configure(state='normal')
        if distribution == "Linear" or distribution == "Power" or distribution == "Polynomial":
            self.ent_lag2_close_lag.configure(state='normal')
        if distribution == "Polynomial":
            self.ent_lag2_extra_frq.configure(state='normal')
            self.ent_lag2_extra_lag.configure(state='normal')

    def lag_enable_disable_3(self, distribution):
        ent_lag3_list = [self.ent_lag3_start_frq, self.ent_lag3_start_lag, self.ent_lag3_close_frq, self.ent_lag3_close_lag, self.ent_lag3_extra_frq, self.ent_lag3_extra_lag]
        [x.configure(state='disable') for x in ent_lag3_list]
        if distribution != "(Unused)":
            self.ent_lag3_start_frq.configure(state='normal')
            self.ent_lag3_start_lag.configure(state='normal')
            self.ent_lag3_close_frq.configure(state='normal')
        if distribution == "Linear" or distribution == "Power" or distribution == "Polynomial":
            self.ent_lag3_close_lag.configure(state='normal')
        if distribution == "Polynomial":
            self.ent_lag3_extra_frq.configure(state='normal')
            self.ent_lag3_extra_lag.configure(state='normal')

    def lag_enable_disable_4(self, distribution):
        ent_lag4_list = [self.ent_lag4_start_frq, self.ent_lag4_start_lag, self.ent_lag4_close_frq, self.ent_lag4_close_lag, self.ent_lag4_extra_frq, self.ent_lag4_extra_lag]
        [x.configure(state='disable') for x in ent_lag4_list]
        if distribution != "(Unused)":
            self.ent_lag4_start_frq.configure(state='normal')
            self.ent_lag4_start_lag.configure(state='normal')
            self.ent_lag4_close_frq.configure(state='normal')
        if distribution == "Linear" or distribution == "Power" or distribution == "Polynomial":
            self.ent_lag4_close_lag.configure(state='normal')
        if distribution == "Polynomial":
            self.ent_lag4_extra_frq.configure(state='normal')
            self.ent_lag4_extra_lag.configure(state='normal')

    def lag_enable_disable_5(self, distribution):
        ent_lag5_list = [self.ent_lag5_start_frq, self.ent_lag5_start_lag, self.ent_lag5_close_frq, self.ent_lag5_close_lag, self.ent_lag5_extra_frq, self.ent_lag5_extra_lag]
        [x.configure(state='disable') for x in ent_lag5_list]
        if distribution != "(Unused)":
            self.ent_lag5_start_frq.configure(state='normal')
            self.ent_lag5_start_lag.configure(state='normal')
            self.ent_lag5_close_frq.configure(state='normal')
        if distribution == "Linear" or distribution == "Power" or distribution == "Polynomial":
            self.ent_lag5_close_lag.configure(state='normal')
        if distribution == "Polynomial":
            self.ent_lag5_extra_frq.configure(state='normal')
            self.ent_lag5_extra_lag.configure(state='normal')

    def lag_enable_disable_6(self, distribution):
        ent_lag6_list = [self.ent_lag6_start_frq, self.ent_lag6_start_lag, self.ent_lag6_close_frq, self.ent_lag6_close_lag, self.ent_lag6_extra_frq, self.ent_lag6_extra_lag]
        [x.configure(state='disable') for x in ent_lag6_list]
        if distribution != "(Unused)":
            self.ent_lag6_start_frq.configure(state='normal')
            self.ent_lag6_start_lag.configure(state='normal')
            self.ent_lag6_close_frq.configure(state='normal')
        if distribution == "Linear" or distribution == "Power" or distribution == "Polynomial":
            self.ent_lag6_close_lag.configure(state='normal')
        if distribution == "Polynomial":
            self.ent_lag6_extra_frq.configure(state='normal')
            self.ent_lag6_extra_lag.configure(state='normal')

    ## ----- Menu Functions -----
    ## 'Load' Dialog:
    def loadfile(self):
        path = tkinter.filedialog.askopenfile(filetypes = (("Text files", "*.txt"), ("All files", "*.*"))).name

        self.title('CorrSim - ' + path)

        with open(path, 'r') as f:
            content = f.readlines()

            ## Delete all current entries
            [widget.delete(0, tk.END)        for widget in self.frm_numbins.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.deselect()               for widget in self.frm_numbins.winfo_children() if isinstance(widget, tk.Checkbutton)]
            [widget.config(state='normal')   for widget in self.frm_source.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.delete(0, tk.END)        for widget in self.frm_source.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.deselect()               for widget in self.frm_source.winfo_children() if isinstance(widget, tk.Checkbutton)]
            [widget.delete(0, tk.END)        for widget in self.frm_fourier.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.delete(0, tk.END)        for widget in self.frm_ps.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.delete(0, tk.END)        for widget in self.frm_ps_1.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.config(state='normal')   for widget in self.frm_ps_2.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.delete(0, tk.END)        for widget in self.frm_ps_2.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.config(state='normal')   for widget in self.frm_ps_3.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.delete(0, tk.END)        for widget in self.frm_ps_3.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.config(state='normal')   for widget in self.frm_lags1.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.delete(0, tk.END)        for widget in self.frm_lags1.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.config(state='normal')   for widget in self.frm_plotting.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.delete(0, tk.END)        for widget in self.frm_plotting.winfo_children() if isinstance(widget, tk.Entry)]
            [widget.deselect()               for widget in self.frm_plotting.winfo_children() if isinstance(widget, tk.Checkbutton)]
            [widget.delete(0, tk.END)        for widget in self.frm_names.winfo_children() if isinstance(widget, tk.Entry)]

            ## Read in loaded file
            self.ent_output_dir.insert(0, content[0][0:-1])
            self.ent_fileprefix.insert(0, content[1][0:-1])
            self.ent_obs_length.insert(0, content[2][0:-1])
            self.ent_time_res.insert(0, content[3][0:-1])
            if content[4] == "1\n": self.chk_optimise.select()

            self.ent_mean_counts_A.insert(0, content[5][0:-1])
            self.ent_mean_counts_B.insert(0, content[6][0:-1])
            self.ent_F_rms_A.insert(0, content[7][0:-1])
            self.ent_F_rms_B.insert(0, content[8][0:-1])

            if content[9]  == "1\n": self.chk_apply_red_A.select()
            if content[10] == "1\n": self.chk_apply_red_B.select()
            self.ent_red_frms_A.insert(0,            content[11][0:-1])
            self.ent_red_frms_B.insert(0,            content[12][0:-1])
            self.ent_red_slope_A.insert(0,           content[13][0:-1])
            self.ent_red_slope_B.insert(0,           content[14][0:-1])

            if content[15] == "1\n": self.chk_apply_poisson_A.select()
            if content[16] == "1\n": self.chk_apply_poisson_B.select()

            if content[17] == "1\n": self.chk_apply_readout_A.select()
            if content[18] == "1\n": self.chk_apply_readout_B.select()
            self.ent_read_noise_A.insert(0,            content[19][0:-1])
            self.ent_read_noise_B.insert(0,            content[20][0:-1])

            if content[21] == "1\n": self.chk_apply_scintillation_A.select()
            if content[22] == "1\n": self.chk_apply_scintillation_B.select()

            self.ent_telescope_diameter_A.insert(0,    content[23][0:-1])
            self.ent_telescope_altitude_A.insert(0,    content[24][0:-1])
            self.ent_exposure_time_A.insert(0,         content[25][0:-1])
            self.ent_turbulence_height_A.insert(0,     content[26][0:-1])
            self.ent_target_altitude_A.insert(0,       content[27][0:-1])
            self.ent_empirical_value_A.insert(0,       content[28][0:-1])
            self.ent_telescope_diameter_B.insert(0,    content[29][0:-1])
            self.ent_telescope_altitude_B.insert(0,    content[30][0:-1])
            self.ent_exposure_time_B.insert(0,         content[31][0:-1])
            self.ent_turbulence_height_B.insert(0,     content[32][0:-1])
            self.ent_target_altitude_B.insert(0,       content[33][0:-1])
            self.ent_empirical_value_B.insert(0,       content[34][0:-1])


            if content[35][0:-1]  == "2": self.rad_lor.select()
            else: self.rad_bpl.select()
            self.ent_ps_power.insert(0,              content[36][0:-1])
            self.ent_break_freq.insert(0,            content[37][0:-1])
            self.ent_coh_constant.insert(0,          content[38][0:-1])
            self.ent_coh_power.insert(0,             content[39][0:-1])
            self.ent_coh_break_freq.insert(0,        content[40][0:-1])

            self.ent_lorA_1_Nrm.insert(0,            content[41][0:-1])
            self.ent_lorA_2_Nrm.insert(0,            content[42][0:-1])
            self.ent_lorA_3_Nrm.insert(0,            content[43][0:-1])
            self.ent_lorA_4_Nrm.insert(0,            content[44][0:-1])
            self.ent_lorA_5_Nrm.insert(0,            content[45][0:-1])
            self.ent_lorA_6_Nrm.insert(0,            content[46][0:-1])
            self.ent_lorA_1_Wdt.insert(0,            content[47][0:-1])
            self.ent_lorA_2_Wdt.insert(0,            content[48][0:-1])
            self.ent_lorA_3_Wdt.insert(0,            content[49][0:-1])
            self.ent_lorA_4_Wdt.insert(0,            content[50][0:-1])
            self.ent_lorA_5_Wdt.insert(0,            content[51][0:-1])
            self.ent_lorA_6_Wdt.insert(0,            content[52][0:-1])
            self.ent_lorA_1_Mid.insert(0,            content[53][0:-1])
            self.ent_lorA_2_Mid.insert(0,            content[54][0:-1])
            self.ent_lorA_3_Mid.insert(0,            content[55][0:-1])
            self.ent_lorA_4_Mid.insert(0,            content[56][0:-1])
            self.ent_lorA_5_Mid.insert(0,            content[57][0:-1])
            self.ent_lorA_6_Mid.insert(0,            content[58][0:-1])
            self.ent_lorB_1_Nrm.insert(0,            content[59][0:-1])
            self.ent_lorB_2_Nrm.insert(0,            content[60][0:-1])
            self.ent_lorB_3_Nrm.insert(0,            content[61][0:-1])
            self.ent_lorB_4_Nrm.insert(0,            content[62][0:-1])
            self.ent_lorB_5_Nrm.insert(0,            content[63][0:-1])
            self.ent_lorB_6_Nrm.insert(0,            content[64][0:-1])
            self.ent_lorB_1_Wdt.insert(0,            content[65][0:-1])
            self.ent_lorB_2_Wdt.insert(0,            content[66][0:-1])
            self.ent_lorB_3_Wdt.insert(0,            content[67][0:-1])
            self.ent_lorB_4_Wdt.insert(0,            content[68][0:-1])
            self.ent_lorB_5_Wdt.insert(0,            content[69][0:-1])
            self.ent_lorB_6_Wdt.insert(0,            content[70][0:-1])
            self.ent_lorB_1_Mid.insert(0,            content[71][0:-1])
            self.ent_lorB_2_Mid.insert(0,            content[72][0:-1])
            self.ent_lorB_3_Mid.insert(0,            content[73][0:-1])
            self.ent_lorB_4_Mid.insert(0,            content[74][0:-1])
            self.ent_lorB_5_Mid.insert(0,            content[75][0:-1])
            self.ent_lorB_6_Mid.insert(0,            content[76][0:-1])
            self.ent_lorB_1_Coh.insert(0,            content[77][0:-1])
            self.ent_lorB_2_Coh.insert(0,            content[78][0:-1])
            self.ent_lorB_3_Coh.insert(0,            content[79][0:-1])
            self.ent_lorB_4_Coh.insert(0,            content[80][0:-1])
            self.ent_lorB_5_Coh.insert(0,            content[81][0:-1])
            self.ent_lorB_6_Coh.insert(0,            content[82][0:-1])

            if content[83][0:-1]  == "2": self.rad_phase.select()
            else: self.rad_time.select()

            self.ent_overall_lag.insert(0,           content[84][0:-1])

            self.lag1_distribtn.set(content[85][0:-1])
            self.lag2_distribtn.set(content[86][0:-1])
            self.lag3_distribtn.set(content[87][0:-1])
            self.lag4_distribtn.set(content[88][0:-1])
            self.lag5_distribtn.set(content[89][0:-1])
            self.lag6_distribtn.set(content[90][0:-1])

            self.ent_lag1_start_frq.insert(0,        content[91][0:-1])
            self.ent_lag2_start_frq.insert(0,        content[92][0:-1])
            self.ent_lag3_start_frq.insert(0,        content[93][0:-1])
            self.ent_lag4_start_frq.insert(0,        content[94][0:-1])
            self.ent_lag5_start_frq.insert(0,        content[95][0:-1])
            self.ent_lag6_start_frq.insert(0,        content[96][0:-1])
            self.ent_lag1_start_lag.insert(0,        content[97][0:-1])
            self.ent_lag2_start_lag.insert(0,        content[98][0:-1])
            self.ent_lag3_start_lag.insert(0,        content[99][0:-1])
            self.ent_lag4_start_lag.insert(0,        content[100][0:-1])
            self.ent_lag5_start_lag.insert(0,        content[101][0:-1])
            self.ent_lag6_start_lag.insert(0,        content[102][0:-1])
            self.ent_lag1_close_frq.insert(0,        content[103][0:-1])
            self.ent_lag2_close_frq.insert(0,        content[104][0:-1])
            self.ent_lag3_close_frq.insert(0,        content[105][0:-1])
            self.ent_lag4_close_frq.insert(0,        content[106][0:-1])
            self.ent_lag5_close_frq.insert(0,        content[107][0:-1])
            self.ent_lag6_close_frq.insert(0,        content[108][0:-1])
            self.ent_lag1_close_lag.insert(0,        content[109][0:-1])
            self.ent_lag2_close_lag.insert(0,        content[110][0:-1])
            self.ent_lag3_close_lag.insert(0,        content[111][0:-1])
            self.ent_lag4_close_lag.insert(0,        content[112][0:-1])
            self.ent_lag5_close_lag.insert(0,        content[113][0:-1])
            self.ent_lag6_close_lag.insert(0,        content[114][0:-1])
            self.ent_lag1_extra_frq.insert(0,        content[115][0:-1])
            self.ent_lag2_extra_frq.insert(0,        content[116][0:-1])
            self.ent_lag3_extra_frq.insert(0,        content[117][0:-1])
            self.ent_lag4_extra_frq.insert(0,        content[118][0:-1])
            self.ent_lag5_extra_frq.insert(0,        content[119][0:-1])
            self.ent_lag6_extra_frq.insert(0,        content[120][0:-1])
            self.ent_lag1_extra_lag.insert(0,        content[121][0:-1])
            self.ent_lag2_extra_lag.insert(0,        content[122][0:-1])
            self.ent_lag3_extra_lag.insert(0,        content[123][0:-1])
            self.ent_lag4_extra_lag.insert(0,        content[124][0:-1])
            self.ent_lag5_extra_lag.insert(0,        content[125][0:-1])
            self.ent_lag6_extra_lag.insert(0,        content[126][0:-1])


            if content[127] == "1\n": self.chk_plot_rednoise.select()

            if content[128] == "1\n": self.chk_plot_model_ccf.select()
            if content[129] == "1\n": self.chk_calculate_ccf.select()
            self.ent_ccf_segment.insert(0,           content[130][0:-1])
            self.ent_ccf_binning.insert(0,           content[131][0:-1])

            if content[132] == "1\n": self.chk_plot_model_fourier.select()
            if content[133] == "1\n": self.chk_calculate_fourier.select()
            self.ent_fourier_bins.insert(0,          content[134][0:-1])
            self.ent_fourier_rebins.insert(0,        content[135][0:-1])
            self.ent_reference_freq.insert(0,        content[136][0:-1])

            if content[137] == "1\n": self.chk_remove_whitenoise.select()


            ## Check which fields need to be enabled/disabled
            self.red_enable_disable_A()
            self.red_enable_disable_B()

            self.read_enable_disable_A()
            self.read_enable_disable_B()

            self.scin_enable_disable_A()
            self.scin_enable_disable_B()

            self.psd_enable_disable()

            self.lag_enable_disable_1(content[85][0:-1])
            self.lag_enable_disable_2(content[86][0:-1])
            self.lag_enable_disable_3(content[87][0:-1])
            self.lag_enable_disable_4(content[88][0:-1])
            self.lag_enable_disable_5(content[89][0:-1])
            self.lag_enable_disable_6(content[90][0:-1])

            self.ccf_enable_disable()
            self.fourier_enable_disable()




    ## Save files:
    def savefile_prime(self, path):
        with open(path, 'w') as f:
            ## Save all current entries
            f.write(self.output_dir.get()+"\n")
            f.write(self.fileprefix.get()+"\n")
            f.write(self.obs_length.get()+"\n")
            f.write(self.time_res.get()+"\n")
            f.write(str(self.optimise.get())+"\n")

            f.write(self.mean_counts_A.get()+"\n")
            f.write(self.mean_counts_B.get()+"\n")
            f.write(self.F_rms_A.get()+"\n")
            f.write(self.F_rms_B.get()+"\n")

            f.write(str(self.apply_red_A.get())+"\n")
            f.write(str(self.apply_red_B.get())+"\n")
            f.write(self.F_rms_rednoise_A.get()+"\n")
            f.write(self.F_rms_rednoise_B.get()+"\n")
            f.write(self.rednoise_slope_A.get()+"\n")
            f.write(self.rednoise_slope_B.get()+"\n")

            f.write(str(self.apply_poisson_A.get())+"\n")
            f.write(str(self.apply_poisson_B.get())+"\n")

            f.write(str(self.apply_readout_A.get())+"\n")
            f.write(str(self.apply_readout_B.get())+"\n")
            f.write(self.read_noise_A.get()+"\n")
            f.write(self.read_noise_B.get()+"\n")

            f.write(str(self.apply_scintillation_A.get())+"\n")
            f.write(str(self.apply_scintillation_B.get())+"\n")
            f.write(self.telescope_diameter_A.get()+"\n")
            f.write(self.telescope_altitude_A.get()+"\n")
            f.write(self.exposure_time_A.get()+"\n")
            f.write(self.turbulence_height_A.get()+"\n")
            f.write(self.target_altitude_A.get()+"\n")
            f.write(self.empirical_value_A.get()+"\n")
            f.write(self.telescope_diameter_B.get()+"\n")
            f.write(self.telescope_altitude_B.get()+"\n")
            f.write(self.exposure_time_B.get()+"\n")
            f.write(self.turbulence_height_B.get()+"\n")
            f.write(self.target_altitude_B.get()+"\n")
            f.write(self.empirical_value_B.get()+"\n")

            f.write(str(self.power_model_type.get())+"\n")

            f.write(self.ps_power.get()+"\n")
            f.write(self.break_freq.get()+"\n")
            f.write(self.coh_constant.get()+"\n")
            f.write(self.coh_power.get()+"\n")
            f.write(self.coh_break_freq.get()+"\n")

            f.write(self.lorA_1_Nrm.get()+"\n")
            f.write(self.lorA_2_Nrm.get()+"\n")
            f.write(self.lorA_3_Nrm.get()+"\n")
            f.write(self.lorA_4_Nrm.get()+"\n")
            f.write(self.lorA_5_Nrm.get()+"\n")
            f.write(self.lorA_6_Nrm.get()+"\n")
            f.write(self.lorA_1_Wdt.get()+"\n")
            f.write(self.lorA_2_Wdt.get()+"\n")
            f.write(self.lorA_3_Wdt.get()+"\n")
            f.write(self.lorA_4_Wdt.get()+"\n")
            f.write(self.lorA_5_Wdt.get()+"\n")
            f.write(self.lorA_6_Wdt.get()+"\n")
            f.write(self.lorA_1_Mid.get()+"\n")
            f.write(self.lorA_2_Mid.get()+"\n")
            f.write(self.lorA_3_Mid.get()+"\n")
            f.write(self.lorA_4_Mid.get()+"\n")
            f.write(self.lorA_5_Mid.get()+"\n")
            f.write(self.lorA_6_Mid.get()+"\n")
            f.write(self.lorB_1_Nrm.get()+"\n")
            f.write(self.lorB_2_Nrm.get()+"\n")
            f.write(self.lorB_3_Nrm.get()+"\n")
            f.write(self.lorB_4_Nrm.get()+"\n")
            f.write(self.lorB_5_Nrm.get()+"\n")
            f.write(self.lorB_6_Nrm.get()+"\n")
            f.write(self.lorB_1_Wdt.get()+"\n")
            f.write(self.lorB_2_Wdt.get()+"\n")
            f.write(self.lorB_3_Wdt.get()+"\n")
            f.write(self.lorB_4_Wdt.get()+"\n")
            f.write(self.lorB_5_Wdt.get()+"\n")
            f.write(self.lorB_6_Wdt.get()+"\n")
            f.write(self.lorB_1_Mid.get()+"\n")
            f.write(self.lorB_2_Mid.get()+"\n")
            f.write(self.lorB_3_Mid.get()+"\n")
            f.write(self.lorB_4_Mid.get()+"\n")
            f.write(self.lorB_5_Mid.get()+"\n")
            f.write(self.lorB_6_Mid.get()+"\n")
            f.write(self.lorB_1_Coh.get()+"\n")
            f.write(self.lorB_2_Coh.get()+"\n")
            f.write(self.lorB_3_Coh.get()+"\n")
            f.write(self.lorB_4_Coh.get()+"\n")
            f.write(self.lorB_5_Coh.get()+"\n")
            f.write(self.lorB_6_Coh.get()+"\n")

            f.write(str(self.time_or_phase.get())+"\n")

            f.write(self.overall_lag.get()+"\n")

            f.write(self.lag1_distribtn.get()+"\n")
            f.write(self.lag2_distribtn.get()+"\n")
            f.write(self.lag3_distribtn.get()+"\n")
            f.write(self.lag4_distribtn.get()+"\n")
            f.write(self.lag5_distribtn.get()+"\n")
            f.write(self.lag6_distribtn.get()+"\n")
            f.write(self.lag1_start_frq.get()+"\n")
            f.write(self.lag2_start_frq.get()+"\n")
            f.write(self.lag3_start_frq.get()+"\n")
            f.write(self.lag4_start_frq.get()+"\n")
            f.write(self.lag5_start_frq.get()+"\n")
            f.write(self.lag6_start_frq.get()+"\n")
            f.write(self.lag1_start_lag.get()+"\n")
            f.write(self.lag2_start_lag.get()+"\n")
            f.write(self.lag3_start_lag.get()+"\n")
            f.write(self.lag4_start_lag.get()+"\n")
            f.write(self.lag5_start_lag.get()+"\n")
            f.write(self.lag6_start_lag.get()+"\n")
            f.write(self.lag1_close_frq.get()+"\n")
            f.write(self.lag2_close_frq.get()+"\n")
            f.write(self.lag3_close_frq.get()+"\n")
            f.write(self.lag4_close_frq.get()+"\n")
            f.write(self.lag5_close_frq.get()+"\n")
            f.write(self.lag6_close_frq.get()+"\n")
            f.write(self.lag1_close_lag.get()+"\n")
            f.write(self.lag2_close_lag.get()+"\n")
            f.write(self.lag3_close_lag.get()+"\n")
            f.write(self.lag4_close_lag.get()+"\n")
            f.write(self.lag5_close_lag.get()+"\n")
            f.write(self.lag6_close_lag.get()+"\n")
            f.write(self.lag1_extra_frq.get()+"\n")
            f.write(self.lag2_extra_frq.get()+"\n")
            f.write(self.lag3_extra_frq.get()+"\n")
            f.write(self.lag4_extra_frq.get()+"\n")
            f.write(self.lag5_extra_frq.get()+"\n")
            f.write(self.lag6_extra_frq.get()+"\n")
            f.write(self.lag1_extra_lag.get()+"\n")
            f.write(self.lag2_extra_lag.get()+"\n")
            f.write(self.lag3_extra_lag.get()+"\n")
            f.write(self.lag4_extra_lag.get()+"\n")
            f.write(self.lag5_extra_lag.get()+"\n")
            f.write(self.lag6_extra_lag.get()+"\n")


            f.write(str(self.plot_rednoise.get())+"\n")

            f.write(str(self.plot_model_ccf.get())+"\n")
            f.write(str(self.calculate_ccf.get())+"\n")
            f.write(self.ccf_segment.get()+"\n")
            f.write(self.ccf_binning.get()+"\n")

            f.write(str(self.plot_model_fourier.get())+"\n")
            f.write(str(self.calculate_fourier.get())+"\n")
            f.write(self.fourier_bins.get()+"\n")
            f.write(self.fourier_rebins.get()+"\n")
            f.write(self.reference_freq.get()+"\n")

            f.write(str(self.remove_whitenoise.get())+"\n")




    ## 'Save' Dialog:
    def savefile(self):
        ## If not yet saved, then default to 'Save As'
        try:
            path = self.title().split(' - ')[1][0:]
        except:
            path = ''
        if path != '':
            self.savefile_prime(path)
        else:
            self.savefileas()

    ## 'Save As' Dialog:
    def savefileas(self):
        try:
            path = tkinter.filedialog.asksaveasfile(filetypes = (("Text files", "*.txt"), ("All files", "*.*"))).name
            self.title('CorrSim - ' + path)
        except:
            return
        self.savefile_prime(path)


    ## Optimise bins:
    def optimize_bins(self, num_bins, time_res):
        direction = "d"
        iterations = 1
        while corrfunc.largest_prime_factor(num_bins) > 7:
            if direction == "d":
                num_bins -= iterations
                direction = "u"
                iterations += 1
            else:
                num_bins += iterations
                direction = "d"
                iterations += 1

        obs_length2 = num_bins*float(time_res)

        return [num_bins, obs_length2]

    ## Calculate the number of bins:
    def calculate_num_bins(self):
        obs_length    = self.ent_obs_length.get()
        time_res      = self.ent_time_res.get()

        num_bins      = int(np.floor(float(obs_length)/float(time_res)))

        optimise_yn   = self.optimise.get()
        if optimise_yn == 1:
            optimised  = self.optimize_bins(num_bins, time_res)
            num_bins   = optimised[0]
            obs_length = optimised[1]

        if num_bins < 1E6:
            self.lbl_numbin_result["text"]    = f"{num_bins} bins"
            self.lbl_numbin_result["fg"]      = "black"
            self.lbl_numbin_wrning["text"]    = ""
        else:
            self.lbl_numbin_result["text"]    = f"{num_bins} bins"
            self.lbl_numbin_result["fg"]      = "red"
            self.lbl_numbin_wrning["text"]    = "WARNING!"

        self.ent_obs_length.delete(0, tk.END)
        self.ent_obs_length.insert(0, obs_length)

    ## Test Model Power Spectra
    def test_PSDs(self, output_dir, fileprefix, obs_length, time_res, power_model_type, \
        ps_power, break_freq, coh_constant, coh_power, coh_break_freq, \
        lorA_1_Nrm, lorA_2_Nrm, lorA_3_Nrm, lorA_4_Nrm, lorA_5_Nrm, lorA_6_Nrm, \
        lorA_1_Wdt, lorA_2_Wdt, lorA_3_Wdt, lorA_4_Wdt, lorA_5_Wdt, lorA_6_Wdt, \
        lorA_1_Mid, lorA_2_Mid, lorA_3_Mid, lorA_4_Mid, lorA_5_Mid, lorA_6_Mid, \
        lorB_1_Nrm, lorB_2_Nrm, lorB_3_Nrm, lorB_4_Nrm, lorB_5_Nrm, lorB_6_Nrm, \
        lorB_1_Wdt, lorB_2_Wdt, lorB_3_Wdt, lorB_4_Wdt, lorB_5_Wdt, lorB_6_Wdt, \
        lorB_1_Mid, lorB_2_Mid, lorB_3_Mid, lorB_4_Mid, lorB_5_Mid, lorB_6_Mid, \
        lorB_1_Coh, lorB_2_Coh, lorB_3_Coh, lorB_4_Coh, lorB_5_Coh, lorB_6_Coh):

        self.handle_powspec_errors(power_model_type,\
            ps_power, break_freq, coh_constant, coh_power, coh_break_freq, \
            lorA_1_Nrm, lorA_2_Nrm, lorA_3_Nrm, lorA_4_Nrm, lorA_5_Nrm, lorA_6_Nrm, \
            lorA_1_Wdt, lorA_2_Wdt, lorA_3_Wdt, lorA_4_Wdt, lorA_5_Wdt, lorA_6_Wdt, \
            lorB_1_Nrm, lorB_2_Nrm, lorB_3_Nrm, lorB_4_Nrm, lorB_5_Nrm, lorB_6_Nrm, \
            lorB_1_Wdt, lorB_2_Wdt, lorB_3_Wdt, lorB_4_Wdt, lorB_5_Wdt, lorB_6_Wdt, \
            lorB_1_Coh, lorB_2_Coh, lorB_3_Coh, lorB_4_Coh, lorB_5_Coh, lorB_6_Coh)

        ## Power Spectra

        print("")
        print("-----------")

        if power_model_type == 2: power_model_type = "lorentzians"
        else: power_model_type = "broken_powerlaw"

        ## Construct Lorentzian Array

        Lorentz_params_A    = [
            ##        Norm        Width        Midpoint
                    lorA_1_Nrm,    lorA_1_Wdt,    lorA_1_Mid
                ,    lorA_2_Nrm,    lorA_2_Wdt,    lorA_2_Mid
                ,    lorA_3_Nrm,    lorA_3_Wdt,    lorA_3_Mid
                ,    lorA_4_Nrm,    lorA_4_Wdt,    lorA_4_Mid
                ,    lorA_5_Nrm,    lorA_5_Wdt,    lorA_5_Mid
                ,    lorA_6_Nrm,    lorA_6_Wdt,    lorA_6_Mid
        ]

        Lorentz_params_B    = [
            ##        Norm        Width        Midpoint    Coherence
                    lorB_1_Nrm,    lorB_1_Wdt,    lorB_1_Mid, lorB_1_Coh
                ,    lorB_2_Nrm,    lorB_2_Wdt,    lorB_2_Mid, lorB_2_Coh
                ,    lorB_3_Nrm,    lorB_3_Wdt,    lorB_3_Mid, lorB_3_Coh
                ,    lorB_4_Nrm,    lorB_4_Wdt,    lorB_4_Mid, lorB_4_Coh
                ,    lorB_5_Nrm,    lorB_5_Wdt,    lorB_5_Mid, lorB_5_Coh
                ,    lorB_6_Nrm,    lorB_6_Wdt,    lorB_6_Mid, lorB_6_Coh
        ]

        corrpsds_test.TestPSDs(output_dir, fileprefix, obs_length, time_res, power_model_type,
            ps_power, break_freq, coh_constant, coh_power, coh_break_freq, Lorentz_params_A, Lorentz_params_B)

        print("-----------")

    def handle_powspec_errors(self, power_model_type,\
        ps_power, break_freq, coh_constant, coh_power, coh_break_freq, \
        lorA_1_Nrm, lorA_2_Nrm, lorA_3_Nrm, lorA_4_Nrm, lorA_5_Nrm, lorA_6_Nrm, \
        lorA_1_Wdt, lorA_2_Wdt, lorA_3_Wdt, lorA_4_Wdt, lorA_5_Wdt, lorA_6_Wdt, \
        lorB_1_Nrm, lorB_2_Nrm, lorB_3_Nrm, lorB_4_Nrm, lorB_5_Nrm, lorB_6_Nrm, \
        lorB_1_Wdt, lorB_2_Wdt, lorB_3_Wdt, lorB_4_Wdt, lorB_5_Wdt, lorB_6_Wdt, \
        lorB_1_Coh, lorB_2_Coh, lorB_3_Coh, lorB_4_Coh, lorB_5_Coh, lorB_6_Coh):

        if power_model_type == 1:
            if break_freq <= 0:        raise Exception("Pow. Spec. Break Freq. must be greater than 0")
            if coh_constant < 0:    raise Exception("Coherence Constant must be greater than or equal to 0")
            if coh_break_freq <= 0:    raise Exception("Coherence Break Freq. must be greater than 0")
        if power_model_type == 2:
            if lorA_1_Nrm == 0:        raise Exception("Please fully define at least one Lorentzian in each band")
            if lorA_1_Wdt == 0:        raise Exception("Please fully define at least one Lorentzian in each band")
            if lorB_1_Nrm == 0:        raise Exception("Please fully define at least one Lorentzian in each band")
            if lorB_1_Wdt == 0:        raise Exception("Please fully define at least one Lorentzian in each band")
            if  lorA_1_Nrm < 0 or lorB_1_Nrm < 0 or \
                lorA_2_Nrm < 0 or lorB_2_Nrm < 0 or \
                lorA_3_Nrm < 0 or lorB_3_Nrm < 0 or \
                lorA_4_Nrm < 0 or lorB_4_Nrm < 0 or \
                lorA_5_Nrm < 0 or lorB_5_Nrm < 0 or \
                lorA_6_Nrm < 0 or lorB_6_Nrm < 0:    raise Exception("Lorentzian Normalisations must be >= 0")
            if  lorA_1_Wdt < 0 or lorB_1_Wdt < 0 or \
                lorA_2_Wdt < 0 or lorB_2_Wdt < 0 or \
                lorA_3_Wdt < 0 or lorB_3_Wdt < 0 or \
                lorA_4_Wdt < 0 or lorB_4_Wdt < 0 or \
                lorA_5_Wdt < 0 or lorB_5_Wdt < 0 or \
                lorA_6_Wdt < 0 or lorB_6_Wdt < 0:    raise Exception("Lorentzian Widths must be >= 0")
            if  lorB_1_Coh < 0 or lorB_1_Coh > 1 or \
                lorB_2_Coh < 0 or lorB_2_Coh > 1 or \
                lorB_3_Coh < 0 or lorB_3_Coh > 1 or \
                lorB_4_Coh < 0 or lorB_4_Coh > 1 or \
                lorB_5_Coh < 0 or lorB_5_Coh > 1 or \
                lorB_6_Coh < 0 or lorB_6_Coh > 1:    raise Exception("Lorentzian Coherences must be: 0 <= Coh <= 1")
            if  (lorA_1_Nrm != 0 and lorA_1_Wdt == 0) or \
                (lorA_2_Nrm != 0 and lorA_2_Wdt == 0) or \
                (lorA_3_Nrm != 0 and lorA_3_Wdt == 0) or \
                (lorA_4_Nrm != 0 and lorA_4_Wdt == 0) or \
                (lorA_5_Nrm != 0 and lorA_5_Wdt == 0) or \
                (lorA_6_Nrm != 0 and lorA_6_Wdt == 0) or \
                (lorB_1_Nrm != 0 and lorB_1_Wdt == 0) or \
                (lorB_2_Nrm != 0 and lorB_2_Wdt == 0) or \
                (lorB_3_Nrm != 0 and lorB_3_Wdt == 0) or \
                (lorB_4_Nrm != 0 and lorB_4_Wdt == 0) or \
                (lorB_5_Nrm != 0 and lorB_5_Wdt == 0) or \
                (lorB_6_Nrm != 0 and lorB_6_Wdt == 0):    raise Exception("All Lorentzians with Normalisation > 0 must have Width > 0")
            if  (lorA_1_Nrm == 0 and lorA_1_Wdt != 0) or \
                (lorA_2_Nrm == 0 and lorA_2_Wdt != 0) or \
                (lorA_3_Nrm == 0 and lorA_3_Wdt != 0) or \
                (lorA_4_Nrm == 0 and lorA_4_Wdt != 0) or \
                (lorA_5_Nrm == 0 and lorA_5_Wdt != 0) or \
                (lorA_6_Nrm == 0 and lorA_6_Wdt != 0) or \
                (lorB_1_Nrm == 0 and lorB_1_Wdt != 0) or \
                (lorB_2_Nrm == 0 and lorB_2_Wdt != 0) or \
                (lorB_3_Nrm == 0 and lorB_3_Wdt != 0) or \
                (lorB_4_Nrm == 0 and lorB_4_Wdt != 0) or \
                (lorB_5_Nrm == 0 and lorB_5_Wdt != 0) or \
                (lorB_6_Nrm == 0 and lorB_6_Wdt != 0):    raise Exception("All Lorentzians with Width > 0 must have Normalisation > 0")


    ## Test Model Lags
    def test_lags(self, output_dir, fileprefix, obs_length, time_res, time_or_phase, overall_lag,
        lag1_distribtn, lag2_distribtn, lag3_distribtn, lag4_distribtn, lag5_distribtn, lag6_distribtn, \
        lag1_start_frq, lag2_start_frq, lag3_start_frq, lag4_start_frq, lag5_start_frq, lag6_start_frq, \
        lag1_start_lag, lag2_start_lag, lag3_start_lag, lag4_start_lag, lag5_start_lag, lag6_start_lag, \
        lag1_close_frq, lag2_close_frq, lag3_close_frq, lag4_close_frq, lag5_close_frq, lag6_close_frq, \
        lag1_close_lag, lag2_close_lag, lag3_close_lag, lag4_close_lag, lag5_close_lag, lag6_close_lag, \
        lag1_extra_frq, lag2_extra_frq, lag3_extra_frq, lag4_extra_frq, lag5_extra_frq, lag6_extra_frq, \
        lag1_extra_lag, lag2_extra_lag, lag3_extra_lag, lag4_extra_lag, lag5_extra_lag, lag6_extra_lag):

        print("")
        print("-----------")

        model_lags_full = self.make_model_lags(time_or_phase, overall_lag,
            lag1_distribtn, lag2_distribtn, lag3_distribtn, lag4_distribtn, lag5_distribtn, lag6_distribtn, \
            lag1_start_frq, lag2_start_frq, lag3_start_frq, lag4_start_frq, lag5_start_frq, lag6_start_frq, \
            lag1_start_lag, lag2_start_lag, lag3_start_lag, lag4_start_lag, lag5_start_lag, lag6_start_lag, \
            lag1_close_frq, lag2_close_frq, lag3_close_frq, lag4_close_frq, lag5_close_frq, lag6_close_frq, \
            lag1_close_lag, lag2_close_lag, lag3_close_lag, lag4_close_lag, lag5_close_lag, lag6_close_lag, \
            lag1_extra_frq, lag2_extra_frq, lag3_extra_frq, lag4_extra_frq, lag5_extra_frq, lag6_extra_frq, \
            lag1_extra_lag, lag2_extra_lag, lag3_extra_lag, lag4_extra_lag, lag5_extra_lag, lag6_extra_lag)

        time_or_phase    = model_lags_full[0]
        model_lag_array    = model_lags_full[1]

        corrlags_test.TestLags(output_dir, fileprefix, obs_length, time_res, time_or_phase, overall_lag, model_lag_array)

        print("-----------")



    ## Make model_lag_Array
    def make_model_lags(self, time_or_phase, overall_lag,
        lag1_distribtn, lag2_distribtn, lag3_distribtn, lag4_distribtn, lag5_distribtn, lag6_distribtn, \
        lag1_start_frq, lag2_start_frq, lag3_start_frq, lag4_start_frq, lag5_start_frq, lag6_start_frq, \
        lag1_start_lag, lag2_start_lag, lag3_start_lag, lag4_start_lag, lag5_start_lag, lag6_start_lag, \
        lag1_close_frq, lag2_close_frq, lag3_close_frq, lag4_close_frq, lag5_close_frq, lag6_close_frq, \
        lag1_close_lag, lag2_close_lag, lag3_close_lag, lag4_close_lag, lag5_close_lag, lag6_close_lag, \
        lag1_extra_frq, lag2_extra_frq, lag3_extra_frq, lag4_extra_frq, lag5_extra_frq, lag6_extra_frq, \
        lag1_extra_lag, lag2_extra_lag, lag3_extra_lag, lag4_extra_lag, lag5_extra_lag, lag6_extra_lag):

        if time_or_phase == 2: time_or_phase = "phase"
        else: time_or_phase = "time"

        ## --- Check to make sure all the distributions are okay ---

        ## Are the lags out of order?
        if  (lag2_start_frq != 0 and lag1_start_frq > lag2_start_frq and lag2_distribtn != "(Unused)") or \
            (lag3_start_frq != 0 and lag2_start_frq > lag3_start_frq and lag3_distribtn != "(Unused)") or \
            (lag4_start_frq != 0 and lag3_start_frq > lag4_start_frq and lag4_distribtn != "(Unused)") or \
            (lag5_start_frq != 0 and lag4_start_frq > lag5_start_frq and lag5_distribtn != "(Unused)") or \
            (lag6_start_frq != 0 and lag5_start_frq > lag6_start_frq and lag6_distribtn != "(Unused)"):
            raise Exception("Please define lag distributions in order of increasing frequency (Freq. 1)")


        ## Are any frequencies in the wrong order?
        if  (lag2_start_frq < lag1_close_frq and lag2_close_frq > lag1_start_frq and lag2_distribtn != "(Unused)") or \
            (lag3_start_frq < lag2_close_frq and lag3_close_frq > lag2_start_frq and lag3_distribtn != "(Unused)") or \
            (lag4_start_frq < lag3_close_frq and lag4_close_frq > lag3_start_frq and lag4_distribtn != "(Unused)") or \
            (lag5_start_frq < lag4_close_frq and lag5_close_frq > lag4_start_frq and lag5_distribtn != "(Unused)") or \
            (lag6_start_frq < lag5_close_frq and lag6_close_frq > lag5_start_frq and lag6_distribtn != "(Unused)"):
            print("WARNING: Model lags overlap. Earlier distributions will be favoured.")

        ## Are any frequencies < 0?
        if  (lag1_start_frq < 0 and lag1_distribtn != "(Unused)") or \
            (lag1_close_frq < 0 and lag1_distribtn != "(Unused)") or \
            (lag1_extra_frq < 0 and lag1_distribtn != "(Unused)") or \
            (lag2_start_frq < 0 and lag2_distribtn != "(Unused)") or \
            (lag2_close_frq < 0 and lag2_distribtn != "(Unused)") or \
            (lag2_extra_frq < 0 and lag2_distribtn != "(Unused)") or \
            (lag3_start_frq < 0 and lag3_distribtn != "(Unused)") or \
            (lag3_close_frq < 0 and lag3_distribtn != "(Unused)") or \
            (lag3_extra_frq < 0 and lag3_distribtn != "(Unused)") or \
            (lag4_start_frq < 0 and lag4_distribtn != "(Unused)") or \
            (lag4_close_frq < 0 and lag4_distribtn != "(Unused)") or \
            (lag4_extra_frq < 0 and lag4_distribtn != "(Unused)") or \
            (lag5_start_frq < 0 and lag5_distribtn != "(Unused)") or \
            (lag5_close_frq < 0 and lag5_distribtn != "(Unused)") or \
            (lag5_extra_frq < 0 and lag5_distribtn != "(Unused)") or \
            (lag6_start_frq < 0 and lag6_distribtn != "(Unused)") or \
            (lag6_close_frq < 0 and lag6_distribtn != "(Unused)") or \
            (lag6_extra_frq < 0 and lag6_distribtn != "(Unused)"):
            raise Exception("All lag frequencies must be greater than or equal to 0")

        ## When using time lags, do any lag distributions go between positive and negative numbers, or have a zero?
        if time_or_phase == "time":
            # print(lag1_distribtn, lag1_start_lag, lag1_close_lag, lag1_extra_lag)
            # print(lag2_distribtn, lag2_start_lag, lag2_close_lag, lag2_extra_lag)
            # print(lag3_distribtn, lag3_start_lag, lag3_close_lag, lag3_extra_lag)
            # print(lag4_distribtn, lag4_start_lag, lag4_close_lag, lag4_extra_lag)
            # print(lag5_distribtn, lag5_start_lag, lag5_close_lag, lag5_extra_lag)
            # print(lag6_distribtn, lag6_start_lag, lag6_close_lag, lag6_extra_lag)
            if  lag1_start_lag * lag1_close_lag <= 0 and lag1_distribtn != "Const. Time" and lag1_distribtn != "Const. Phase" and lag1_distribtn != "(Unused)" or \
                lag1_start_lag * lag1_extra_lag <= 0 and lag1_distribtn == "Polynomial" or \
                lag2_start_lag * lag2_close_lag <= 0 and lag2_distribtn != "Const. Time" and lag2_distribtn != "Const. Phase" and lag2_distribtn != "(Unused)" or \
                lag2_start_lag * lag2_extra_lag <= 0 and lag2_distribtn == "Polynomial" or \
                lag3_start_lag * lag3_close_lag <= 0 and lag3_distribtn != "Const. Time" and lag3_distribtn != "Const. Phase" and lag3_distribtn != "(Unused)" or \
                lag3_start_lag * lag3_extra_lag <= 0 and lag3_distribtn == "Polynomial" or \
                lag4_start_lag * lag4_close_lag <= 0 and lag4_distribtn != "Const. Time" and lag4_distribtn != "Const. Phase" and lag4_distribtn != "(Unused)" or \
                lag4_start_lag * lag4_extra_lag <= 0 and lag4_distribtn == "Polynomial" or \
                lag5_start_lag * lag5_close_lag <= 0 and lag5_distribtn != "Const. Time" and lag5_distribtn != "Const. Phase" and lag5_distribtn != "(Unused)" or \
                lag5_start_lag * lag5_extra_lag <= 0 and lag5_distribtn == "Polynomial" or \
                lag6_start_lag * lag6_close_lag <= 0 and lag6_distribtn != "Const. Time" and lag6_distribtn != "Const. Phase" and lag6_distribtn != "(Unused)" or \
                lag6_start_lag * lag6_extra_lag <= 0 and lag6_distribtn == "Polynomial":

                raise Exception("When using time lags, distributions cannot go to or across zero.")

        ## Do any linear distributions have freq=0?
        if  (lag1_distribtn == "Linear" and lag1_start_frq <= 0) or \
            (lag2_distribtn == "Linear" and lag2_start_frq <= 0) or \
            (lag3_distribtn == "Linear" and lag3_start_frq <= 0) or \
            (lag4_distribtn == "Linear" and lag4_start_frq <= 0) or \
            (lag5_distribtn == "Linear" and lag5_start_frq <= 0) or \
            (lag6_distribtn == "Linear" and lag6_start_frq <= 0):
            raise Exception("For Linear distributions, Freqs. 1 and 2 must both be greater than 0")

        ## Do any linear distributions have freq=0?
        if  (lag1_distribtn == "Power" and ((lag1_start_frq >= lag1_close_frq) or (lag1_start_frq == 0))) or \
            (lag2_distribtn == "Power" and ((lag2_start_frq >= lag2_close_frq) or (lag2_start_frq == 0))) or \
            (lag3_distribtn == "Power" and ((lag3_start_frq >= lag3_close_frq) or (lag3_start_frq == 0))) or \
            (lag4_distribtn == "Power" and ((lag4_start_frq >= lag4_close_frq) or (lag4_start_frq == 0))) or \
            (lag5_distribtn == "Power" and ((lag5_start_frq >= lag5_close_frq) or (lag5_start_frq == 0))) or \
            (lag6_distribtn == "Power" and ((lag6_start_frq >= lag6_close_frq) or (lag6_start_frq == 0))):
            raise Exception("For Power distributions, Freqs. 1 and 2 must both be greater than 0")


        ## Are any polynomial frequencies set to zero?
        if  (lag1_distribtn == "Polynomial" and (lag1_start_frq <= 0 or lag1_close_frq <= 0 or lag1_extra_frq <= 0)) or \
            (lag2_distribtn == "Polynomial" and (lag2_start_frq <= 0 or lag2_close_frq <= 0 or lag2_extra_frq <= 0)) or \
            (lag3_distribtn == "Polynomial" and (lag3_start_frq <= 0 or lag3_close_frq <= 0 or lag3_extra_frq <= 0)) or \
            (lag4_distribtn == "Polynomial" and (lag4_start_frq <= 0 or lag4_close_frq <= 0 or lag4_extra_frq <= 0)) or \
            (lag5_distribtn == "Polynomial" and (lag5_start_frq <= 0 or lag5_close_frq <= 0 or lag5_extra_frq <= 0)) or \
            (lag6_distribtn == "Polynomial" and (lag6_start_frq <= 0 or lag6_close_frq <= 0 or lag6_extra_frq <= 0)):
            raise Exception("For Polynomial distributions, all frequencies must be greater than 0")

        ## Are any frequencies in the wrong order?
        if  (lag1_start_frq > lag1_close_frq and lag1_distribtn != "(Unused)") or \
            (lag2_start_frq > lag2_close_frq and lag2_distribtn != "(Unused)") or \
            (lag3_start_frq > lag3_close_frq and lag3_distribtn != "(Unused)") or \
            (lag4_start_frq > lag4_close_frq and lag4_distribtn != "(Unused)") or \
            (lag5_start_frq > lag5_close_frq and lag5_distribtn != "(Unused)") or \
            (lag6_start_frq > lag6_close_frq and lag6_distribtn != "(Unused)"):
            raise Exception("In all cases, Freq. 2 must be greater than or equal to Freq 1")

        ## --- End of error handling ---

        ## Construct Lag Array
        model_lag_array = np.zeros((6, 7))    ## Sets up the array.
        if time_or_phase == "phase":
            ##    - Set the lag model using corrfunc.lag_section_phase
            ##                                                distribution    start_freq,    close_freq,    start_lag,    close_lag=0
            model_lag_array[0,] = corrfunc.lag_section_phase(lag1_distribtn,    lag1_start_frq, lag1_start_lag, lag1_close_frq, lag1_close_lag, lag1_extra_frq, lag1_extra_lag)
            model_lag_array[1,] = corrfunc.lag_section_phase(lag2_distribtn,    lag2_start_frq, lag2_start_lag, lag2_close_frq, lag2_close_lag, lag2_extra_frq, lag2_extra_lag)
            model_lag_array[2,] = corrfunc.lag_section_phase(lag3_distribtn,    lag3_start_frq, lag3_start_lag, lag3_close_frq, lag3_close_lag, lag3_extra_frq, lag3_extra_lag)
            model_lag_array[3,] = corrfunc.lag_section_phase(lag4_distribtn,    lag4_start_frq, lag4_start_lag, lag4_close_frq, lag4_close_lag, lag4_extra_frq, lag4_extra_lag)
            model_lag_array[4,] = corrfunc.lag_section_phase(lag5_distribtn,    lag5_start_frq, lag5_start_lag, lag5_close_frq, lag5_close_lag, lag5_extra_frq, lag5_extra_lag)
            model_lag_array[5,] = corrfunc.lag_section_phase(lag6_distribtn,    lag6_start_frq, lag6_start_lag, lag6_close_frq, lag6_close_lag, lag6_extra_frq, lag6_extra_lag)

        else:
            ##    - Set the lag model using corrfunc.lag_section_time.
            ##    - NOTE: Time lags cannot go to or across zero. Please split distributions instead.
            ##                                                distribution        start_freq,    close_freq,    start_lag,    close_lag=0
            model_lag_array[0,] = corrfunc.lag_section_time(lag1_distribtn,    lag1_start_frq, lag1_start_lag, lag1_close_frq, lag1_close_lag, lag1_extra_frq, lag1_extra_lag)
            model_lag_array[1,] = corrfunc.lag_section_time(lag2_distribtn,    lag2_start_frq, lag2_start_lag, lag2_close_frq, lag2_close_lag, lag2_extra_frq, lag2_extra_lag)
            model_lag_array[2,] = corrfunc.lag_section_time(lag3_distribtn,    lag3_start_frq, lag3_start_lag, lag3_close_frq, lag3_close_lag, lag3_extra_frq, lag3_extra_lag)
            model_lag_array[3,] = corrfunc.lag_section_time(lag4_distribtn,    lag4_start_frq, lag4_start_lag, lag4_close_frq, lag4_close_lag, lag4_extra_frq, lag4_extra_lag)
            model_lag_array[4,] = corrfunc.lag_section_time(lag5_distribtn,    lag5_start_frq, lag5_start_lag, lag5_close_frq, lag5_close_lag, lag5_extra_frq, lag5_extra_lag)
            model_lag_array[5,] = corrfunc.lag_section_time(lag6_distribtn,    lag6_start_frq, lag6_start_lag, lag6_close_frq, lag6_close_lag, lag6_extra_frq, lag6_extra_lag)

        return time_or_phase, model_lag_array


    ## Execute CorrSim
    def execute_CorrSim(self, output_dir, fileprefix, obs_length, time_res, mean_counts_A, mean_counts_B, F_rms_A, F_rms_B, \
        apply_red_A, apply_red_B, F_rms_rednoise_A, F_rms_rednoise_B, rednoise_slope_A, rednoise_slope_B, \
        apply_poisson_A, apply_poisson_B, apply_readout_A, apply_readout_B, read_noise_A, read_noise_B, \
        apply_scintillation_A, apply_scintillation_B, \
        empirical_value_A, telescope_diameter_A, exposure_time_A, \
        target_altitude_A, telescope_altitude_A, turbulence_height_A, \
        empirical_value_B, telescope_diameter_B, exposure_time_B, \
        target_altitude_B, telescope_altitude_B, turbulence_height_B, \
        power_model_type, ps_power, break_freq, coh_constant, coh_power, coh_break_freq, \
        lorA_1_Nrm,    lorA_2_Nrm,    lorA_3_Nrm,    lorA_4_Nrm,    lorA_5_Nrm,    lorA_6_Nrm, \
        lorA_1_Wdt,    lorA_2_Wdt,    lorA_3_Wdt,    lorA_4_Wdt,    lorA_5_Wdt,    lorA_6_Wdt, \
        lorA_1_Mid,    lorA_2_Mid,    lorA_3_Mid,    lorA_4_Mid,    lorA_5_Mid,    lorA_6_Mid, \
        lorB_1_Nrm,    lorB_2_Nrm,    lorB_3_Nrm,    lorB_4_Nrm,    lorB_5_Nrm,    lorB_6_Nrm, \
        lorB_1_Wdt,    lorB_2_Wdt,    lorB_3_Wdt,    lorB_4_Wdt,    lorB_5_Wdt,    lorB_6_Wdt, \
        lorB_1_Mid,    lorB_2_Mid,    lorB_3_Mid,    lorB_4_Mid,    lorB_5_Mid,    lorB_6_Mid, \
        lorB_1_Coh,    lorB_2_Coh,    lorB_3_Coh,    lorB_4_Coh,    lorB_5_Coh,    lorB_6_Coh, \
        time_or_phase,    overall_lag, \
        lag1_distribtn,    lag2_distribtn,    lag3_distribtn,    lag4_distribtn,    lag5_distribtn,    lag6_distribtn, \
        lag1_start_frq,    lag2_start_frq,    lag3_start_frq,    lag4_start_frq,    lag5_start_frq,    lag6_start_frq, \
        lag1_start_lag,    lag2_start_lag,    lag3_start_lag,    lag4_start_lag,    lag5_start_lag,    lag6_start_lag, \
        lag1_close_frq,    lag2_close_frq,    lag3_close_frq,    lag4_close_frq,    lag5_close_frq,    lag6_close_frq, \
        lag1_close_lag,    lag2_close_lag,    lag3_close_lag,    lag4_close_lag,    lag5_close_lag,    lag6_close_lag, \
        lag1_extra_frq,    lag2_extra_frq,    lag3_extra_frq,    lag4_extra_frq,    lag5_extra_frq,    lag6_extra_frq, \
        lag1_extra_lag,    lag2_extra_lag,    lag3_extra_lag,    lag4_extra_lag,    lag5_extra_lag,    lag6_extra_lag, \
        plot_rednoise, plot_model_ccf, calculate_ccf, ccf_segment, ccf_binning, \
        plot_model_fourier, calculate_fourier, fourier_bins, fourier_rebins, reference_freq, remove_whitenoise):

        ## Check for errors:
        if obs_length < time_res:        raise Exception("Time Resolution must be less than Observation Length")
        if obs_length <= 0:              raise Exception("Observation Length must be greater than 0")
        if time_res <= 0:                raise Exception("Time Resolution must be greater than 0")
        if telescope_diameter_A <= 0:    raise Exception("Telescope Diameter must be greater than 0")
        if telescope_diameter_B <= 0:    raise Exception("Telescope Diameter must be greater than 0")
        if exposure_time_A > time_res and apply_scintillation_A == 1:    raise Exception("Exposure Time (in Scintillation Noise) must be less than Time Resolution")
        if exposure_time_B > time_res and apply_scintillation_B == 1:    raise Exception("Exposure Time (in Scintillation Noise) must be less than Time Resolution")
        if exposure_time_A < 0:          raise Exception("Exposure Time must be greater than 0")
        if exposure_time_B < 0:          raise Exception("Exposure Time must be greater than 0")
        if target_altitude_A < 0 or target_altitude_A > 90:    raise Exception("Target Altitude must be between 0 and 90")
        if target_altitude_B < 0 or target_altitude_B > 90:    raise Exception("Target Altitude must be between 0 and 90")
        if telescope_altitude_A < 0:     raise Exception("Telescope Altitude must be greater than or equal to 0")
        if telescope_altitude_B < 0:     raise Exception("Telescope Altitude must be greater than or equal to 0")
        if turbulence_height_A <= 0:     raise Exception("Turbulence Height must be greater than 0")
        if turbulence_height_B <= 0:     raise Exception("Turbulence Height must be greater than 0")
        if empirical_value_A < 0:        raise Exception("Empirical Value must be greater than 1")
        if empirical_value_B < 0:        raise Exception("Empirical Value must be greater than 1")

        self.handle_powspec_errors(power_model_type,\
            ps_power, break_freq, coh_constant, coh_power, coh_break_freq, \
            lorA_1_Nrm, lorA_2_Nrm, lorA_3_Nrm, lorA_4_Nrm, lorA_5_Nrm, lorA_6_Nrm, \
            lorA_1_Wdt, lorA_2_Wdt, lorA_3_Wdt, lorA_4_Wdt, lorA_5_Wdt, lorA_6_Wdt, \
            lorB_1_Nrm, lorB_2_Nrm, lorB_3_Nrm, lorB_4_Nrm, lorB_5_Nrm, lorB_6_Nrm, \
            lorB_1_Wdt, lorB_2_Wdt, lorB_3_Wdt, lorB_4_Wdt, lorB_5_Wdt, lorB_6_Wdt, \
            lorB_1_Coh, lorB_2_Coh, lorB_3_Coh, lorB_4_Coh, lorB_5_Coh, lorB_6_Coh)

        if  lag1_start_frq < 0 or lag1_close_frq < 0 or lag1_extra_frq < 0 or \
            lag2_start_frq < 0 or lag2_close_frq < 0 or lag2_extra_frq < 0 or \
            lag3_start_frq < 0 or lag3_close_frq < 0 or lag3_extra_frq < 0 or \
            lag4_start_frq < 0 or lag4_close_frq < 0 or lag4_extra_frq < 0 or \
            lag5_start_frq < 0 or lag5_close_frq < 0 or lag5_extra_frq < 0 or \
            lag6_start_frq < 0 or lag6_close_frq < 0 or lag6_extra_frq < 0:        raise Exception("All lag frequencies must be greater than or equal to 0")

        ## CCF:
        if calculate_ccf == 1:
            if ccf_segment < time_res:                  raise Exception("CCF Segment Size must be greater than Time Resolution")
            if ccf_segment > obs_length:                raise Exception("CCF Segment Size must be less than Observation Length")
            if ccf_binning > (ccf_segment/time_res):    raise Exception("CCF Binning must be less than number of bins per segment (CCF Segment Size / Time Resolution)")
        if calculate_fourier == 1:
            if fourier_bins > (obs_length/time_res):    raise Exception("Fourier Segment Size must be less than the total number of bins")
            if fourier_rebins < 0:                      raise Exception("Fourier Rebinning must be greater than or equal to 0")
            if reference_freq <= 0:                     raise Exception("Reference Frequency must be greater than or equal to 0")

        ## Red Noise
        if apply_red_A == 1 and apply_red_B == 1:   apply_red = "AB"
        elif apply_red_A == 1 and apply_red_B != 1: apply_red = "A"
        elif apply_red_A != 1 and apply_red_B == 1: apply_red = "B"
        else: apply_red = "n"

        ## Make sure that it doesn't try to plot any rednoise if there isn't any
        if apply_red_A == 0 and apply_red_B == 0:
            plot_rednoise = 0

        if apply_poisson_A == 1 and apply_poisson_B == 1:   apply_poisson = "AB"
        elif apply_poisson_A == 1 and apply_poisson_B != 1: apply_poisson = "A"
        elif apply_poisson_A != 1 and apply_poisson_B == 1: apply_poisson = "B"
        else: apply_poisson = "n"

        if apply_scintillation_A == 1 and apply_scintillation_B == 1:   apply_scintillation = "AB"
        elif apply_scintillation_A == 1 and apply_scintillation_B != 1: apply_scintillation = "A"
        elif apply_scintillation_A != 1 and apply_scintillation_B == 1: apply_scintillation = "B"
        else: apply_scintillation = "n"

        if apply_readout_A == 1 and apply_readout_B == 1:   apply_readout = "AB"
        elif apply_readout_A == 1 and apply_readout_B != 1: apply_readout = "A"
        elif apply_readout_A != 1 and apply_readout_B == 1: apply_readout = "B"
        else: apply_readout = "n"

        if apply_red_A == 1 and F_rms_A <= F_rms_rednoise_A:
            raise Exception("F_rms_A must be greater than F_rms_rednoise_A.")
        if apply_red_B == 1 and F_rms_B <= F_rms_rednoise_B:
            raise Exception("F_rms_B must be greater than F_rms_rednoise_B.")

        scin_noise_A = 0
        if apply_scintillation_A == 1:
            scin_noise_A = corrfunc.calculate_scintillation_noise(empirical_value_A, telescope_diameter_A,\
                exposure_time_A, target_altitude_A, telescope_altitude_A, turbulence_height_A)

        scin_noise_B = 0
        if apply_scintillation_B == 1:
            scin_noise_B = corrfunc.calculate_scintillation_noise(empirical_value_B, telescope_diameter_B,\
                exposure_time_B, target_altitude_B, telescope_altitude_B, turbulence_height_B)

        ## Power Spectra

        if power_model_type == 2: power_model_type = "lorentzians"
        else: power_model_type = "broken_powerlaw"

        ## Construct Lorentzian Array

        Lorentz_params_A    = [
            ##        Norm        Width        Midpoint
                     lorA_1_Nrm,    lorA_1_Wdt,    lorA_1_Mid
                ,    lorA_2_Nrm,    lorA_2_Wdt,    lorA_2_Mid
                ,    lorA_3_Nrm,    lorA_3_Wdt,    lorA_3_Mid
                ,    lorA_4_Nrm,    lorA_4_Wdt,    lorA_4_Mid
                ,    lorA_5_Nrm,    lorA_5_Wdt,    lorA_5_Mid
                ,    lorA_6_Nrm,    lorA_6_Wdt,    lorA_6_Mid
        ]

        Lorentz_params_B    = [
            ##        Norm        Width        Midpoint    Coherence
                     lorB_1_Nrm,    lorB_1_Wdt,    lorB_1_Mid, lorB_1_Coh
                ,    lorB_2_Nrm,    lorB_2_Wdt,    lorB_2_Mid, lorB_2_Coh
                ,    lorB_3_Nrm,    lorB_3_Wdt,    lorB_3_Mid, lorB_3_Coh
                ,    lorB_4_Nrm,    lorB_4_Wdt,    lorB_4_Mid, lorB_4_Coh
                ,    lorB_5_Nrm,    lorB_5_Wdt,    lorB_5_Mid, lorB_5_Coh
                ,    lorB_6_Nrm,    lorB_6_Wdt,    lorB_6_Mid, lorB_6_Coh
        ]

        model_lags_full = self.make_model_lags(time_or_phase, overall_lag,
            lag1_distribtn, lag2_distribtn, lag3_distribtn, lag4_distribtn, lag5_distribtn, lag6_distribtn, \
            lag1_start_frq, lag2_start_frq, lag3_start_frq, lag4_start_frq, lag5_start_frq, lag6_start_frq, \
            lag1_start_lag, lag2_start_lag, lag3_start_lag, lag4_start_lag, lag5_start_lag, lag6_start_lag, \
            lag1_close_frq, lag2_close_frq, lag3_close_frq, lag4_close_frq, lag5_close_frq, lag6_close_frq, \
            lag1_close_lag, lag2_close_lag, lag3_close_lag, lag4_close_lag, lag5_close_lag, lag6_close_lag, \
            lag1_extra_frq, lag2_extra_frq, lag3_extra_frq, lag4_extra_frq, lag5_extra_frq, lag6_extra_frq, \
            lag1_extra_lag, lag2_extra_lag, lag3_extra_lag, lag4_extra_lag, lag5_extra_lag, lag6_extra_lag)

        time_or_phase      = model_lags_full[0]
        model_lag_array    = model_lags_full[1]

        CorrSim.CorrSim(output_dir, fileprefix, obs_length, time_res, \
            mean_counts_A, mean_counts_B, F_rms_A, F_rms_B, \
            apply_red, F_rms_rednoise_A, F_rms_rednoise_B, rednoise_slope_A, rednoise_slope_B, \
            apply_poisson, apply_readout, read_noise_A, read_noise_B, \
            apply_scintillation, scin_noise_A, scin_noise_B, \
            plot_rednoise, plot_model_ccf, calculate_ccf, ccf_segment, ccf_binning, \
            plot_model_fourier, calculate_fourier, fourier_bins, fourier_rebins, reference_freq, remove_whitenoise, \
            power_model_type, ps_power, break_freq, coh_constant, coh_power, coh_break_freq, \
            Lorentz_params_A, Lorentz_params_B, time_or_phase, overall_lag, model_lag_array)

    ## Shut down the GUI
    def shutdown_ttk_repeat(self):
        # print("Woof!")
        self.quit()
        self.destroy()

    ### --------------- End of Functions ---------------
    ### ================================================

if __name__ == "__main__":
    app = App()
    app.mainloop()
