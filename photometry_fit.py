
import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
from scipy.optimize import curve_fit
from astropy.table import Table
from astropy.constants import c #c=speed of light 
from numpy import loadtxt
from scipy import interpolate
import emcee
import sys
import corner
from scipy.stats import norm
from numpy import where 
from numpy import asarray
from astropy.io import fits
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
import timeit
from astropy.cosmology import FlatLambdaCDM
import pdb
import glob
import os
import platform



### GLOBAL VARIABLES ###

global kappa_eff
global wl_eff
global z
global filt_array
global fname_out_long
global fname_out_long_red
global currentmode
global nparam
global whatprior
global beta
global temperature

wavelengths = 10**np.arange(1,3.5,0.0005) *(u.um)   # Wavelength array for model
frequencies = (c/wavelengths.to(u.m)).value         # Same, converted to frequency
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)  # Cosmological parameters
dist_0 = 100                                        # Distance to adopt (in Mpc) when z = 0

# Kernel Density Estimation of the M31 beta profile (Smith et al., 2012) for the prior
fitsFile = fits.open('Smith_etal_2012_M31_BetaMap.fits', memmap=True)
image = fitsFile[0].data 
image_numbers_only_array = image[~np.isnan(image)]
X = image_numbers_only_array[:,np.newaxis] #np.newaxis = Adds a new axis to the array.
# Bandwidth
n = len(X)  #No. of data points 
d = 1       #Measuring only beta therefore 1D
bandwidth_calc = (n * (d+2) / 4.) ** (-1. / (d+4)) #Silverman's Rule
# Gaussian Kernel Density Estimation
kde = KernelDensity(kernel='cosine', bandwidth=bandwidth_calc).fit(X) 



### SYNTHETIC PHOTOMETRY MODEL ###

def point_source2(frequencies, mass, temperature, z=0.):
        CMB_corrfact = np.ones(( len(frequencies) ))
        rest_freq = frequencies*(1+z)
        if z == 0:
            dist = dist_0 * u.Mpc
        else:
            dist = cosmo.luminosity_distance(z)
            CMB_corrfact = CMB_corrfact - blackbody_nu(rest_freq, cosmo.Tcmb(z))/blackbody_nu(rest_freq, temperature)
        kappa = (kappa_eff / (wl_eff.to(u.m).value ** (-1 * beta))) * ((c.value/rest_freq) ** (-1 * beta))
        try:
            mass.unit
        except:
            mass = mass * u.g
        S_nu = mass * kappa * blackbody_nu(rest_freq, temperature) * (1+z) * CMB_corrfact / dist**2 / 1 * (u.sr)
        return S_nu.to(u.Jy)

def point_source3(frequencies, mass, temperature, beta, z=0.):
        CMB_corrfact = np.ones(( len(frequencies) ))
        rest_freq = frequencies*(1+z)
        if z == 0:
            dist = dist_0 * u.Mpc
        else:
            dist = cosmo.luminosity_distance(z)
            CMB_corrfact = CMB_corrfact - blackbody_nu(rest_freq, cosmo.Tcmb(z))/blackbody_nu(rest_freq, temperature)
        kappa = (kappa_eff / (wl_eff.to(u.m).value ** (-1 * beta))) * ((c.value/rest_freq) ** (-1 * beta))
        try:
            mass.unit
        except:
            mass = mass * u.g
        S_nu = mass * kappa * blackbody_nu(rest_freq, temperature) * (1+z) * CMB_corrfact / dist**2 / 1 * (u.sr)
        return S_nu.to(u.Jy)

# Convolution of model by filters
def read_fileterresponse(frequencies, filters):
        filter_files_all = ("PACS_70mu.dat", 
                            "PACS_100mu.dat", 
                            "PACS_160mu.dat", 
                            "SPIRE_250mu.dat", 
                            "SPIRE_350mu.dat", 
                            "SPIRE_500mu.dat",
                            "SCUBA2_450mu.dat",
                            "SCUBA2_850mu.dat",
                            "ALMA_10_Cycle6.2.dat",
                            "ALMA_9_Cycle6.2.dat",
                            "ALMA_8_Cycle6.2.dat",
                            "ALMA_7_Cycle6.2.dat",
                            "ALMA_6_Cycle6.2.dat",
                            "ALMA_5_Cycle6.2.dat",
                            "ALMA_4_Cycle6.2.dat",
                            "ALMA_3_Cycle6.2.dat")
        filter_names_all = ("PACS70",
                            "PACS100",
                            "PACS160",
                            "SPIRE250",
                            "SPIRE350",
                            "SPIRE500",
                            "SCUBA2_450",
                            "SCUBA2_850",
                            "ALMA_10",
                            "ALMA_9",
                            "ALMA_8",
                            "ALMA_7",
                            "ALMA_6",
                            "ALMA_5",
                            "ALMA_4",
                            "ALMA_3")
        filt_used = []
        for i in range(len(filter_names_all)):
                if filter_names_all[i] in filters:
                        filt_used.append(filter_files_all[i])
        #Convolving modified-bb with the filter response curves 
        filter_array = np.zeros([len(filt_used),len(frequencies)])  #To regrid all filters to the same frequencies
        for i in range(len(filt_used)):
                filename = filt_used[i]
                instr = filename[:4]
                if (instr == "SCUB"): #SCUBA2 filter response curves are in GHz
                        x_col, y_col = loadtxt(filename, usecols = (0,-1), unpack=True) #usecols = (1stcol, lastcol)
                        x_col = (x_col*1e9)  #Unit conversion GHz --> Hz
                else:                        #All other files are in AA
                        x_col, y_col = loadtxt(filename, usecols = (0,-1), unpack=True)
                        x_col = (c/(x_col*1e-10)) # Unit conversion AA --> Hz
                interpolated_function = interpolate.interp1d(x_col, y_col, kind='cubic', fill_value=0., bounds_error=False)
                interpolated_filter = interpolated_function(frequencies) 
                filter_array[i,:] = interpolated_filter 
        return filter_array

#Generating Synthetic Photometry
def synthetic_photometry(blackbody, filter_array):
        syn_phot = np.sum(blackbody*filter_array, axis = 1) / np.sum(filter_array, axis = 1)
        return syn_phot


### BAYESIAN STATISTICS ###

#Priors - our data set. All in log values
def prior(theta): 
        # Parameters for temperature normal distribution (if needed)
        temperature_mean = 40.
        temperature_sigma = 40.
        if (nparam == 2): 
            temperature, masslog = theta
            if (whatprior == 'TB-flat') or (whatprior == 'T-flat'):
                    if cosmo.Tcmb(z).value < temperature < 300 and 35 < masslog < 45:
                            return 0.
            else:
                    if cosmo.Tcmb(z).value < temperature < 300 and 35 < masslog < 45:
                            return norm.logpdf(temperature, temperature_mean, temperature_sigma)
        else:
            temperature, masslog, beta = theta
            #data for Beta from Global KDE 
            probbeta = kde.score_samples(np.array(beta).reshape(1,-1))
            if (whatprior == 'TB-flat'):
                    if cosmo.Tcmb(z).value < temperature < 300 and 35 < masslog < 45 and .5 < beta < 4.:
                            return 0.
            elif (whatprior == 'T-flat'):
                    if cosmo.Tcmb(z).value < temperature < 300 and 35 < masslog < 45 and .5 < beta < 4.:
                            return probbeta
            elif (whatprior == 'B-flat'):
                    if cosmo.Tcmb(z).value < temperature < 300 and 35 < masslog < 45 and .5 < beta < 4.:
                            return norm.logpdf(temperature, temperature_mean, temperature_sigma)
            else:
                    if cosmo.Tcmb(z).value < temperature < 300 and 35 < masslog < 45 and .5 < beta < 4.:
                            return norm.logpdf(temperature, temperature_mean, temperature_sigma) + probbeta
        return -np.inf 

#Building Model using Synthetic Photometry and bb curve equation.   
def likelihood_function(theta, y, yerr): #theta - the input values, x = frequencies, y = our fluxes from the residual profile.
                                         #all in log values
        # Intermediate model: MBB spectrum
        if (nparam == 1):
            masslog = theta
            model = point_source1(frequencies, 10**masslog, z=z)
        elif (nparam ==2): 
            temperature, masslog = theta
            model = point_source2(frequencies, 10**masslog, temperature=temperature, z=z)
        else:
            temperature, masslog, beta = theta
            model = point_source3(frequencies, 10**masslog, temperature=temperature, beta=beta, z=z)
        # Final model (to fit): synthetic photometry
        model = synthetic_photometry(model.value, filter_array)
        #1. Gaussian likelihood function for real detections
        detection = np.where (y >= 3*yerr)
        #pdb.set_trace()
        likelihood = (-0.5) * np.sum ( ((y[detection] - model[detection])**2 / (yerr[detection] ** 2)) + \
					       (np.log(2*np.pi*(yerr[detection]**2))) )         
        #2. Cumulative distribution function for values where we only have an upper limit.
        non_detection = np.where (y < 3*yerr)	
        likelihood = likelihood + np.sum(norm.logcdf(3*yerr[non_detection], model[non_detection], yerr[non_detection]))
        return likelihood

# Posterior
def posterior(theta, y, yerr):
        prob = prior(theta)
        if not np.isfinite(prob):
                return -np.inf
        probablilty = prob + likelihood_function(theta, y, yerr)
        return probablilty



### MCMC FITTING ###

#Running MCMC and Saving the results
def run_emcee(flux_data, flux_unc_data, nwl, mcmcoutput_filename, nparam, nwalkers, steps, burn_in):
        # Set up the Sampler: starting positions for T, log(mass), beta (+ some random perturbation) walkers in the moddle grid.
        masslog = np.log10(10**8 * (u.Msun).to(u.g)) #Initial mass guess: 10^8 Msun
        if (nparam == 2):
            pos=[(28.,masslog) + np.random.rand(nparam) for i in range(nwalkers)] 
        else:
            pos=[(28.,masslog,1.) + np.random.rand(nparam) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers, nparam, posterior, args=(flux_data, flux_unc_data)) #Here we use the flux data as our sample set of data and unc.
        # Run the production chain.
        print("Running MCMC...") 
        sampler.run_mcmc(pos, steps) #steps given in main body.
        print("Done.")
        #Plotting Walker plots
        if (nparam == 3):
            fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 9))
            axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
            axes[0].set_ylabel("$temperature$")
            axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
            axes[1].set_ylabel("$log(mass)$")
            axes[2].plot((sampler.chain[:, :, 2]).T, color="k", alpha=0.4)
            axes[2].set_ylabel("$beta$")
            axes[2].set_xlabel("step number")
            if (currentmode == 'raw'):
                plt.savefig(mcmcplotfolder + 'walkerplot_' + fname_out_long + '.png')
            elif (currentmode == 'red'):
                plt.savefig(mcmcplotfolder + 'walkerplot_' + fname_out_long_red + '.png')
            else:
                pass
        else:
            fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 9))
            axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
            axes[0].set_ylabel("$temperature$")
            axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
            axes[1].set_ylabel("$log(mass)$")
            axes[1].set_xlabel("step number")
            if (currentmode == 'raw'):
                plt.savefig(mcmcplotfolder + 'walkerplot_' + fname_out_long + '.png')
            elif (currentmode == 'red'):
                plt.savefig(mcmcplotfolder + 'walkerplot_' + fname_out_long_red + '.png')
            else:
                pass
        plt.close()        
        #Saving the data for Plotting corner plots and analysis
        samples = sampler.chain[:, burn_in:, :].reshape((-1, nparam)) #sampler.chain[no.of.chains, burn-in-step, no.of.variables].
        #burnin - given in main body - point/step from which we are to consider the data. The point where the walkers all converge as seen in the walker plots.
        #Save the data after the burnin has taken place. 
        f = open(mcmcoutput_filename, "w") #creating and writing(w) to a file. mcmcoutput_filename - given in mainbody.
        f.close()
        f = open(mcmcoutput_filename, "a") #append(a) the created data file 
        for i in range(samples.shape[0]):
            if (nparam == 2):
                outstring = "\t" + str(samples[i,0]) + "\t" + str(samples[i,1]) + "\t" + "\n" 
            else:
                outstring = "\t" + str(samples[i,0]) + "\t" + str(samples[i,1]) + "\t" + str(samples[i,2]) + "\t" + "\n" 
            f.write(outstring)
            #"\t" = tab space. "\n" = start new line. str() = converting/casting a value into a string-Textual data in Python is handled with str objects. We're just printing the data as text format(although they are numbers) into the file.
        f.close()
        #return - In this case since there is already an output generated we do not need a "return".


#Analysing emcee results and plotting corner plots
def analyse_emcee(mcmcoutput_filename, nparam): 
        data = np.loadtxt(mcmcoutput_filename) #Loading the file contaning the T, mass, beta array data from mcmc
        samples = data.reshape([len(data),nparam]) #Reshaping the saved array to a '3xlength' array (in this case it's the same shape as the orginal array) that can be read and used to plot corner plots.
        #nparam - Given in main body. len(data) - length of the data array. 
        #Plotting Corner plots
        if (nparam == 2):
            fig = corner.corner(samples, labels=["$temperature$", "$log(mass)$"])
        else:
            fig = corner.corner(samples, labels=["$temperature$", "$log(mass)$", "$beta$"])
        #Saving plots
        if (currentmode == 'raw'):
                fig.savefig(mcmcplotfolder + 'cornerplot_' + fname_out_long + '.png')
        elif (currentmode == 'red'):
                fig.savefig(mcmcplotfolder + 'cornerplot_' + fname_out_long_red + '.png')
        else:
                print('ERROR (PHOTOMETRY_FIT): Keyword "Currentmode" not recognized. Please choose "raw" or "red".')
        plt.close()
        #Final values and uncertainties for temperature, density and beta (expressed as: value, unc+, unc-)
        #*np.percentile(samples, [16-v0, 50-v1, 84-v2] = symmetrical/nornal distribution arround the median assuming the distribution is gaussian. 50th percentile: median=our result. 84th percentile: +1sigma value. 16th percentile: -1sigma value.
        if (nparam == 2):
            temperature_mcmc, mass_mcmc = np.percentile(samples, [16, 50, 84], axis=0).T
            return temperature_mcmc, mass_mcmc
        else:
            temperature_mcmc, mass_mcmc, beta_mcmc = np.percentile(samples, [16, 50, 84], axis=0).T
            return temperature_mcmc, mass_mcmc, beta_mcmc



################### Main Body - Script ##############################################

if __name__=="__main__":

        currentmode = 'None'
        
        ### FILE I/O: UPDATE EACH TIME ###
	# Nota: If you are doing T-correction, update format of %Treal on lines 407, 494
        kappa_eff = .7 * (u.cm)**2 / (u.g)  # kappa_0 from James et al. 2002
        wl_eff = 850. * (u.um)              # lambda_0 from James et al. 2002

        nparam = 3     # Accepted values: 3 (free beta fit), 2 (fixed beta fit)
        if (nparam == 3):
                pass
        elif (nparam == 2):
            beta = 1.5 # Value of beta to adopt in fixed beta fit
        else:
            print('ERROR (PHOTOMETRY_FIT.PY): NPARAM = ', nparam, '. This quantity needs to be either 2 or 3.')

        sedtype = '1T' # Accepted values: '1T' (single-temperature model),
                       #                  '1T+corr' (single-temperature with CMB heating correction
                       #                  '2T' (two-temperatures model)
                       # NOTA: This refers to the SED to be read, not the fitting model (which is always single-temperature)
        if (sedtype == '1T'):
            wl_min = 100.
        else:
            wl_min = 50.

        whatprior = 'TB-flat' # Accepted values: 'TB-flat' = flat prior in T and beta
                              #                  'T-flat' = flat prior in T, beta prior from M31 (Smith et al. 2012)
                              #                  'B-flat' = flat prior in beta, Gaussian prior in T with 40 K mean and 40 K sigma
                              #                  'nonflat' = Gaussian prior in T with 40 K mean and 40 K sigma, beta prior from M31 (Smith et al. 2012)

        # Composition of SED to fit (comment/uncomment/add as needed):
        #comp = 'MBBtest'
        comp = 'E30R-70.0+BE-30.0'

        # Whether to use the fast or long run
        # 'fast': nwalkers = 30, steps = 1000, burn_in = 200; 'long': nwalkers = 100, steps = 10000, burn_in = 1000
        whatrun = 'fast'

        # Folder paths (syntax depends on operating system)
        plat = platform.system()
        if plat == 'Windows':
            SEDfolder = r"synthetic_SEDs\Photometry\\"         # Synthetic potometry files
            outfolder = r"fit_result_files\\"                  # Fit results files
            mcmcplotfolder = r"plots\MCMC\\"   # MCMC PDF plots
        else:
            SEDfolder = 'synthetic_SEDs/Photometry/'
            outfolder = 'fit_result_files/'
            mcmcplotfolder = 'plots/MCMC/'

        ##########################

        # Input & Output filenames

        fixstring = '_allbd'
        
        if (nparam ==2): 
            fixstring = fixstring + '-beta' + '{:.2f}'.format(beta)
        else:
            fixstring = fixstring + '-freeparams'

        if (whatprior == 'TB-flat'):
                fixstring = fixstring + '_flatpriors'
        elif (whatprior == 'T-flat'):
                fixstring = fixstring + '_betaM31-prior'
        elif (whatprior == 'B-flat'):
                fixstring = fixstring + '_T40K-prior'
        else:
             fixstring = fixstring + '_T40K-betaM31-prior'   

        if (sedtype == '1T'):
            allSED = [fn for fn in glob.glob(SEDfolder + '*Phot*' + comp + '*oneT-20.0K*.dat', recursive = False) if not (((os.path.basename(fn).__contains__('CMBcorr')) or (os.path.basename(fn).__contains__('beta'))))]
            if comp == 'MBBtest':
                    fname_out_raw = 'Fit_' + comp + '_oneT' + fixstring
            else:
                    fname_out_raw = 'Fit_' + comp + '-raw_oneT' + fixstring
            fname_out_red = 'Fit_' + comp + '-red_oneT' + fixstring
        elif (sedtype == '1T+corr'):
            allSED = glob.glob(SEDfolder+'*Phot*' + comp + '*oneT*CMBcorr*.dat', recursive = False)
            if comp == 'MBBfit':
                    fname_out_raw = 'Fit_' + comp + '_oneT-CMBcorr' + fixstring
            else:
                    fname_out_raw = 'Fit_' + comp + '-raw_oneT-CMBcorr' + fixstring
            fname_out_red = 'Fit_' + comp + '-red_oneT-CMBcorr' + fixstring
        elif (sedtype == '2T'):
            allSED =  glob.glob(SEDfolder+'Phot*' + comp + '*fw*.dat', recursive = False)
            if comp == 'MBBfit':
                    fname_out_raw = 'Fit_' + comp + '_twoT' + fixstring
            else:
                    fname_out_raw = 'Fit_' + comp + '-raw_twoT' + fixstring
            fname_out_red = 'Fit_' + comp + '-red_twoT' + fixstring
        else:
            print('SED type %s not recognized. Please use 1T, 1T+corr or 2T.' % {sedtype})

        allSED = np.array(allSED)
        

        #Find SED files
        nsedall = len(allSED)
        if nsedall == 0:
                print()
                print('ERROR (PHOTOMETRY_FIT.PY): No files found using the specified pattern')
                print()
        files_raw = []
        files_red = []
        for i in range(nsedall):
                if '-red' in allSED[i]:
                        if plat == 'Windows':
                                files_red.append( ((allSED[i].split('\\'))[-1])[0:-4] )  # Assuming .dat file extension
                        else:
                                files_red.append( ((allSED[i].split('/'))[-1])[0:-4] )
                else:
                        if plat == 'Windows':
                                files_raw.append( ((allSED[i].split('\\'))[-1])[0:-4] )
                        else:
                                files_raw.append( ((allSED[i].split('/'))[-1])[0:-4] )
        files_red.sort(key=str.lower)
        files_raw.sort(key=str.lower)
        nsed = len(files_raw)

        
        #Set up fit parameters
        #Fast Run: nwalkers=30, steps=1000, burn_in=200 // Long Run: nwalkers=100, steps=10000, burn_in=1000
        if whatrun == 'long':
                nwalkers = 100
                steps = 10000
                burn_in = 1000
                fname_out_raw = fname_out_raw + '_longrun'
                fname_out_red = fname_out_red + '_longrun'
        if whatrun == 'fast':
                nwalkers = 30
                steps = 1000
                burn_in = 200
        mcmcoutput_filename = 'Results_temp.dat'

        
        #"Raw MAC" results
        currentmode = 'raw'
        for i in range(nsed):
                z = float((files_raw[i].split('_')[-1])[1:])
                if (sedtype == '2T'):
                    Tstring = ((files_raw[i].split('_')[-2]).split('-')[-1])
                    Treal = float(Tstring[2:])
                else:
                    Tstring = ((files_raw[i].split('_')[-2]).split('-')[1])
                    Treal = float(Tstring[:-1])
                comp = files_raw[i].split('_')[1]
                        
                # Output filenames for pics #
                if (sedtype == '1T'):
                    fname_out_long = comp + '_oneT-' + Tstring + '_z' + '{:0.2f}'.format(z) + fixstring
                elif (sedtype == '1T+corr'):
                    fname_out_long = comp + '_oneT-' + Tstring + '-CMBcorr_z' + '{:0.2f}'.format(z) + fixstring
                elif (sedtype == '2T'):
                    fname_out_long = comp + '_fw' + str(Treal) + '_z' + '{:0.2f}'.format(z) + fixstring
                else:
                    print('SED type %s not recognized. Please use 1T, 1T+corr or 2T.' % {sedtype})
                filt_all = loadtxt(SEDfolder + files_raw[i] + '.dat', dtype = 'str', delimiter = '|', usecols = (0), unpack=True)
                for j in range(len(filt_all)):
                    filt_all[j] = filt_all[j].strip()
                x, y, dy = loadtxt(SEDfolder + files_raw[i] + '.dat', delimiter = '|', usecols = (1, 2, 3), unpack=True)

                # Filter selection
                wlmask = x >= (1+z) * wl_min
                filtmask00 = np.isin(filt_all, "PACS70")
                filtmask01 = np.isin(filt_all, "PACS100")
                filtmask02 = np.isin(filt_all, "PACS160")
                filtmask03 = np.isin(filt_all, "SPIRE250")
                filtmask04 = np.isin(filt_all, "SPIRE350")
                filtmask05 = np.isin(filt_all, "SPIRE500")
                filtmask06 = np.isin(filt_all, "ALMA_10")
                filtmask07 = np.isin(filt_all, "ALMA_9")
                filtmask08 = np.isin(filt_all, "ALMA_8")
                filtmask09 = np.isin(filt_all, "ALMA_7")
                filtmask10 = np.isin(filt_all, "ALMA_6")
                filtmask11 = np.isin(filt_all, "ALMA_5")
                filtmask12 = np.isin(filt_all, "ALMA_4")
                filtmask13 = np.isin(filt_all, "ALMA_3")
                mask_final = np.logical_and(wlmask,
                                            np.logical_or(filtmask00,
                                                          np.logical_or(filtmask01,
                                                                        np.logical_or(filtmask02,
                                                                                      np.logical_or(filtmask03,
                                                                                                    np.logical_or(filtmask04,
                                                                                                                  np.logical_or(filtmask05,
                                                                                                                                np.logical_or(filtmask06,
                                                                                                                                              np.logical_or(filtmask07,
                                                                                                                                                            np.logical_or(filtmask08,
                                                                                                                                                                          np.logical_or(filtmask09,
                                                                                                                                                                                        np.logical_or(filtmask10,
                                                                                                                                                                                                      np.logical_or(filtmask11,
                                                                                                                                                                                                                    np.logical_or(filtmask12, filtmask13)
                                                                                                                                                                                                      )
                                                                                                                                                                                        )
                                                                                                                                                                          )
                                                                                                                                                            )
                                                                                                                                              )
                                                                                                                                )
                                                                                                                  )
                                                                                                    )
                                                                                      )
                                                                        )
                                                          )
                                            )
                )
                wl_data = x[mask_final]
                flux_data = y[mask_final]
                flux_unc_data = dy[mask_final]
                filt_final = filt_all[mask_final]
                nwl = len(wl_data)
                filter_array = read_fileterresponse(frequencies, filt_final)

                # Actual fit
                print('Start time for ', comp, ', T = ', Treal, ', z = ', z, ' (raw MAC):')
                start_time = timeit.default_timer() # Start timer 
                print (start_time)
                run_emcee(flux_data, flux_unc_data, nwl, mcmcoutput_filename, nparam, nwalkers, steps, burn_in)
                if (nparam ==1):
                        mass_output = analyse_emcee(mcmcoutput_filename, nparam)
                elif (nparam == 2):
                        temperature_output, mass_output = analyse_emcee(mcmcoutput_filename, nparam)
                else:
                        temperature_output, mass_output, beta_output = analyse_emcee(mcmcoutput_filename, nparam)
                elapsed = timeit.default_timer() - start_time # End timer
                print ('end time:', elapsed)
                
                mass_rescaled = 10**(mass_output - np.log10((u.Msun).to(u.g)))
                temperature = temperature_output[1]
                dtemperatureplus = temperature_output[2] - temperature_output[1]
                dtemperatureminus = temperature_output[1] - temperature_output[0]
                mass = mass_rescaled[1]
                dmassplus = mass_rescaled[2] - mass_rescaled[1]
                dmassminus = mass_rescaled[1] - mass_rescaled[0]
                beta = beta_output[1]
                dbetaplus = beta_output[2] - beta_output[1]
                dbetaminus = beta_output[1] - beta_output[0]
                # Output fit results if you want:
                #print('T_fit:', temperature, ' +', dtemperatureplus, ' -', dtemperatureminus)
                #print('mass_fit:', mass, ' +', dmassplus, ' -', dmassminus)
                #print('beta_fit:', beta, ' +', dbetaplus, ' -', dbetaminus)
                        
                if i == 0:
                        fname_out_raw = fname_out_raw + '.dat'
                        f = open(fname_out_raw, "w")#open('test', "w")#
                        bdhdr = 'nbands'.rjust(7)
                        if 'fw' in Tstring:
                                Tstr = 'f_w'
                                f.write('\\     Two-temperature model: ' + Tstring.split('-fw')[0] + '\n')
                        else:
                                Tstr = 'T_real'
                        if (nparam == 1):
                                header = '\\' + 'Composition'.rjust(20)  + ',' +  Tstr.rjust(7) + ',' + 'z'.rjust(7) + ',' + bdhdr + ',' + \
                                    'M_fit'.rjust(10) + ',' + 'dM_fit+'.rjust(10) + ',' + 'dM_fit-'.rjust(10) + "\n"
                        elif (nparam == 2):
                                header = '\\' + 'Composition'.rjust(20)  + ',' +  Tstr.rjust(7) + ',' + 'z'.rjust(7) + ',' + bdhdr + ',' + \
                                    'M_fit'.rjust(10) + ',' + 'dM_fit+'.rjust(10) + ',' + 'dM_fit-'.rjust(10) + ',' + 'T_fit'.rjust(7) + ',' + \
	                            'dT_fit+'.rjust(8) + ',' + 'dT_fit-'.rjust(8) + "\n"
                        else:
                                header = '\\' + 'Composition'.rjust(20)  + ',' +  Tstr.rjust(7) + ',' + 'z'.rjust(7) + ',' + bdhdr + ',' + \
                                    'M_fit'.rjust(10) + ',' + 'dM_fit+'.rjust(10) + ',' + 'dM_fit-'.rjust(10) + ',' + 'T_fit'.rjust(7) + ',' + \
	                            'dT_fit+'.rjust(8) + ',' + 'dT_fit-'.rjust(8) + ',' + 'beta_fit'.rjust(9) + ',' + 'dbeta_fit+'.rjust(11) + ',' + \
                                    'dbeta_fit-'.rjust(11) + "\n"
                        f.write(header)
                        f.close()
                f = open(fname_out_raw, "a")#open('test', "a")#
                bdstring = '%7i' %nwl
        
                if 'fw' in Tstring:
                        if (nparam == 1):
                                outstring = comp.rjust(21) + ',' + '%7.4f' %Treal + ',' + '%7.2f' %z + ',' + bdstring + ',' + \
                                    "%0.3e".rjust(6) %mass + ',' + "%0.3e".rjust(6) %dmassplus + ',' + \
				    "%0.3e".rjust(6) %dmassminus + ',' +  "\n"
                        elif (nparam == 2):
                                outstring = comp.rjust(21) + ',' + '%7.4f' %Treal + ',' + '%7.2f' %z + ',' + bdstring + ',' + \
                                    "%0.3e".rjust(6) %mass + ',' + "%0.3e".rjust(6) %dmassplus + ',' + \
				    "%0.3e".rjust(6) %dmassminus + ',' + "%7.2f" %temperature + ',' + \
        			    "%8.2f" %dtemperatureplus + ',' + "%8.2f" %dtemperatureminus + ',' +  "\n"
                        else:
                                outstring = comp.rjust(21) + ',' + '%7.4f' %Treal + ',' + '%7.2f' %z + ',' + bdstring + ',' + \
                                    "%0.3e".rjust(6) %mass + ',' + "%0.3e".rjust(6) %dmassplus + ',' + \
				    "%0.3e".rjust(6) %dmassminus + ',' + "%7.2f" %temperature + ',' + \
        			    "%8.2f" %dtemperatureplus + ',' + "%8.2f" %dtemperatureminus + ',' + \
				    "%9.2f" %beta + ',' + "%11.2f" %dbetaplus + ',' + "%11.2f" %dbetaminus + "\n"
                else:
                        if (nparam == 1):
                                outstring = comp.rjust(21) + ',' + '%7.1f' %Treal + ',' + '%7.2f' %z + ',' + bdstring + ',' + \
        			    "%0.3e".rjust(6) %mass + ',' + "%0.3e".rjust(6) %dmassplus + ',' + \
				    "%0.3e".rjust(6) %dmassminus + "\n"
                        elif (nparam == 2):
                                outstring = comp.rjust(21) + ',' + '%7.1f' %Treal + ',' + '%7.2f' %z + ',' + bdstring + ',' + \
        			    "%0.3e".rjust(6) %mass + ',' + "%0.3e".rjust(6) %dmassplus + ',' + \
				    "%0.3e".rjust(6) %dmassminus + ',' + "%7.2f" %temperature + ',' + \
        			    "%8.2f" %dtemperatureplus + ',' + "%8.2f" %dtemperatureminus + "\n"
                        else:
                                outstring = comp.rjust(21) + ',' + '%7.1f' %Treal + ',' + '%7.2f' %z + ',' + bdstring + ',' + \
				    "%0.3e".rjust(6) %mass + ',' + "%0.3e".rjust(6) %dmassplus + ',' + \
				    "%0.3e".rjust(6) %dmassminus + ',' + "%7.2f" %temperature + ',' + \
				    "%8.2f" %dtemperatureplus + ',' + "%8.2f" %dtemperatureminus + ',' + \
				    "%9.2f" %beta + ',' + "%11.2f" %dbetaplus + ',' + "%11.2f" %dbetaminus + "\n"
                f.write(outstring)
                f.close()
        print('Raw MAC fits completed for ' + comp)
        print()


        if comp == 'MBBtest':   # No opacity correction needed for MBBtest
                pdb.set_trace()

        
        # Opacity-corrected results #
        currentmode = 'red'
        for i in range(nsed):
                z = float((files_red[i].split('_')[-1])[1:])
                if (sedtype == '2T'):
                        Tstring = ((files_red[i].split('_')[-2]).split('-')[-1])
                        Treal = float(Tstring[2:])
                else:
                        Tstring = ((files_red[i].split('_')[-2]).split('-')[1])
                        Treal = float(Tstring[:-1])
                comp = files_red[i].split('_')[1]

                # Output filenames for pics #
                if (sedtype == '1T'):
                    fname_out_long_red = comp + '_oneT-' + Tstring + '_z' + '{:0.2f}'.format(z) + fixstring
                elif (sedtype == '1T+corr'):
                    fname_out_long_red = comp + '_oneT-' + Tstring + '-CMBcorr_z' + '{:0.2f}'.format(z) + fixstring
                elif (sedtype == '2T'):
                    fname_out_long_red = comp + '_twoT-fw' + str(Treal) + '_z' + '{:0.2f}'.format(z) + fixstring
                else:
                    print('SED type %s not recognized. Please use 1T, 1T+corr or 2T.' % {sedtype})
                    
                filt_all = loadtxt(SEDfolder + files_red[i] + '.dat', dtype = 'str', delimiter = '|', usecols = (0), unpack=True)
                for j in range(len(filt_all)):
                    filt_all[j] = filt_all[j].strip()
                x, y, dy = loadtxt(SEDfolder + files_red[i] + '.dat', delimiter = '|', usecols = (1, 2, 3), unpack=True)

                wlmask = x >= (1+z) * wl_min
                filtmask00 = np.isin(filt_all, "PACS70")
                filtmask01 = np.isin(filt_all, "PACS100")
                filtmask02 = np.isin(filt_all, "PACS160")
                filtmask03 = np.isin(filt_all, "SPIRE250")
                filtmask04 = np.isin(filt_all, "SPIRE350")
                filtmask05 = np.isin(filt_all, "SPIRE500")
                filtmask06 = np.isin(filt_all, "ALMA_10")
                filtmask07 = np.isin(filt_all, "ALMA_9")
                filtmask08 = np.isin(filt_all, "ALMA_8")
                filtmask09 = np.isin(filt_all, "ALMA_7")
                filtmask10 = np.isin(filt_all, "ALMA_6")
                filtmask11 = np.isin(filt_all, "ALMA_5")
                filtmask12 = np.isin(filt_all, "ALMA_4")
                filtmask13 = np.isin(filt_all, "ALMA_3")
                mask_final = np.logical_and(wlmask, np.logical_or(filtmask00, np.logical_or(filtmask01, np.logical_or(filtmask02, np.logical_or(filtmask03, np.logical_or(filtmask04, np.logical_or(filtmask05, np.logical_or(filtmask06, np.logical_or(filtmask07, np.logical_or(filtmask08, np.logical_or(filtmask09, np.logical_or(filtmask10, np.logical_or(filtmask11, np.logical_or(filtmask12, filtmask13))))))))))))))
                wl_data = x[mask_final]
                flux_data = y[mask_final]
                flux_unc_data = dy[mask_final]
                filt_final = filt_all[mask_final]
                nwl = len(wl_data)
                filter_array = read_fileterresponse(frequencies, filt_final)

                # Actual fit
                print('Start time for ', comp, ', T = ', Treal, ', z = ', z, ' (reduced opacity):')
                start_time = timeit.default_timer() #Start timer 
                print(start_time)
                run_emcee(flux_data, flux_unc_data, nwl, mcmcoutput_filename, nparam, nwalkers, steps, burn_in)
                if (nparam == 1):
                        mass_output = analyse_emcee(mcmcoutput_filename, nparam)
                elif (nparam == 2):
                        temperature_output, mass_output = analyse_emcee(mcmcoutput_filename, nparam)
                else:
                        temperature_output, mass_output, beta_output = analyse_emcee(mcmcoutput_filename, nparam)
                elapsed = timeit.default_timer() - start_time #End timer
                print ('end time:', elapsed)

                mass_rescaled = 10**(mass_output - np.log10((u.Msun).to(u.g)))
                temperature = temperature_output[1]
                dtemperatureplus = temperature_output[2] - temperature_output[1]
                dtemperatureminus = temperature_output[1] - temperature_output[0]
                mass = mass_rescaled[1]
                dmassplus = mass_rescaled[2] - mass_rescaled[1]
                dmassminus = mass_rescaled[1] - mass_rescaled[0]
                beta = beta_output[1]
                dbetaplus = beta_output[2] - beta_output[1]
                dbetaminus = beta_output[1] - beta_output[0]
                # Output fit results if you want:
                #print('T_fit:', temperature, ' +', dtemperatureplus, ' -', dtemperatureminus)
                #print('mass_fit:', mass, ' +', dmassplus, ' -', dmassminus)
                #print('beta_fit:', beta, ' +', dbetaplus, ' -', dbetaminus)
                
                if i == 0:
                        fname_out_red = fname_out_red + '.dat'
                        f = open(fname_out_red, "w")#open('test', "a")#
                        bdhdr = 'nbands'.rjust(7)

                        if 'fw' in Tstring:
                                Tstr = 'f_w'
                                f.write('\\     Two-temperature model: ' + Tstring.split('-fw')[0] + '\n')
                        else:
                                Tstr = 'T_real'
                        if (nparam == 1):
                                header = '\\' + 'Composition'.rjust(20)  + ',' +  Tstr.rjust(7) + ',' + 'z'.rjust(7) + ',' + bdhdr + ',' + \
                                    'M_fit'.rjust(10) + ',' + 'dM_fit+'.rjust(10) + ',' + 'dM_fit-'.rjust(10) + "\n"
                        elif (nparam == 2):
                                header = '\\' + 'Composition'.rjust(20)  + ',' +  Tstr.rjust(7) + ',' + 'z'.rjust(7) + ',' + bdhdr + ',' + \
                                    'M_fit'.rjust(10) + ',' + 'dM_fit+'.rjust(10) + ',' + 'dM_fit-'.rjust(10) + ',' + 'T_fit'.rjust(7) + ',' + \
                                    'dT_fit+'.rjust(8) + ',' + 'dT_fit-'.rjust(8) + "\n"
                        else:
                                header = '\\' + 'Composition'.rjust(20)  + ',' +  Tstr.rjust(7) + ',' + 'z'.rjust(7) + ',' + bdhdr + ',' + \
                                    'M_fit'.rjust(10) + ',' + 'dM_fit+'.rjust(10) + ',' + 'dM_fit-'.rjust(10) + ',' + 'T_fit'.rjust(7) + ',' + \
                                    'dT_fit+'.rjust(8) + ',' + 'dT_fit-'.rjust(8) + ',' + 'beta_fit'.rjust(9) + ',' + 'dbeta_fit+'.rjust(11) + ',' + \
                                    'dbeta_fit-'.rjust(11) + "\n"
                        f.write(header)
                        f.close()
                f = open(fname_out_red, "a")#open('test', "a")#
                bdstring = '%7i' %nwl

                if 'fw' in Tstring:
                        if (nparam == 1):
                                outstring = comp.rjust(21) + ',' + '%7.4f' %Treal + ',' + '%7.2f' %z + ',' + bdstring + ',' + \
                                    "%0.3e".rjust(6) %mass + ',' + "%0.3e".rjust(6) %dmassplus + ',' + \
                                    "%0.3e".rjust(6) %dmassminus + ',' +  "\n"
                        elif (nparam == 2):
                                outstring = comp.rjust(21) + ',' + '%7.4f' %Treal + ',' + '%7.2f' %z + ',' + bdstring + ',' + \
                                    "%0.3e".rjust(6) %mass + ',' + "%0.3e".rjust(6) %dmassplus + ',' + \
                                    "%0.3e".rjust(6) %dmassminus + ',' + "%7.2f" %temperature + ',' + \
                                    "%8.2f" %dtemperatureplus + ',' + "%8.2f" %dtemperatureminus + ',' +  "\n"
                        else:
                                outstring = comp.rjust(21) + ',' + '%7.4f' %Treal + ',' + '%7.2f' %z + ',' + bdstring + ',' + \
                                    "%0.3e".rjust(6) %mass + ',' + "%0.3e".rjust(6) %dmassplus + ',' + \
                                    "%0.3e".rjust(6) %dmassminus + ',' + "%7.2f" %temperature + ',' + \
                                    "%8.2f" %dtemperatureplus + ',' + "%8.2f" %dtemperatureminus + ',' + \
                                    "%9.2f" %beta + ',' + "%11.2f" %dbetaplus + ',' + "%11.2f" %dbetaminus + "\n"
                else:
                        if (nparam == 1):
                                outstring = comp.rjust(21) + ',' + '%7.1f' %Treal + ',' + '%7.2f' %z + ',' + bdstring + ',' + \
                                    "%0.3e".rjust(6) %mass + ',' + "%0.3e".rjust(6) %dmassplus + ',' + \
                                    "%0.3e".rjust(6) %dmassminus + "\n"
                        elif (nparam == 2):
                                outstring = comp.rjust(21) + ',' + '%7.1f' %Treal + ',' + '%7.2f' %z + ',' + bdstring + ',' + \
                                    "%0.3e".rjust(6) %mass + ',' + "%0.3e".rjust(6) %dmassplus + ',' + \
                                    "%0.3e".rjust(6) %dmassminus + ',' + "%7.2f" %temperature + ',' + \
                                    "%8.2f" %dtemperatureplus + ',' + "%8.2f" %dtemperatureminus + "\n"
                        else:
                                outstring = comp.rjust(21) + ',' + '%7.1f' %Treal + ',' + '%7.2f' %z + ',' + bdstring + ',' + \
                                    "%0.3e".rjust(6) %mass + ',' + "%0.3e".rjust(6) %dmassplus + ',' + \
                                    "%0.3e".rjust(6) %dmassminus + ',' + "%7.2f" %temperature + ',' + \
                                    "%8.2f" %dtemperatureplus + ',' + "%8.2f" %dtemperatureminus + ',' + \
                                    "%9.2f" %beta + ',' + "%11.2f" %dbetaplus + ',' + "%11.2f" %dbetaminus + "\n"
                f.write(outstring)
                f.close()
        print('Reduced MAC fits completed for ' + comp)
        print()

        #pdb.set_trace()        
