
This Git repository contains the code and data needed to reporduce the results from Fanciullo et al. 2020, including figures. The code is mostly written in IDL 7.0, with one script in Python 3.6. It has not been tested for other distributions.

The scripts perform four different functions: rewriting the original opacity data files into a format that can be used by the rest of the code; use these opacities to create a grid of synthetic SEDs (spectra and photometry); fit the synthetic photometry thus obtained; plot the fit results. The rest of this document explains how to do this step-by-step.



0 --- Obtaining the original opacity files
---------------------------------------

The silicate opacity from Demyk et al. 2017a, b (D17A, D17B) can be downloaded from https://www.sshade.eu/db/stopcoda (one needs an SSHADE profile to download the data). Unpack the files in the MAC_files_original/ subfolder.

The carbon data from Mennella et al. 1998 (M98) can be requested to Vito Mennella at the following address: vito.mennella@inaf.it and unpacked in the subdirectory MAC_files_original/M98/



1 --- Reading and rewriting the opacity files
------------------------------------------

Open IDL in the main directory and run the following:

~~~IDL
    opacity_reprocess, 'M98', /savedata
    opacity_reprocess, 'D17A', /savedata
    opacity_reprocess, 'D17B', /savedata
~~~

The output files are saved in the *MAC_files_reprocessed/* subdirectory with the filename *{material name}_{temperature}K_interpol.xcat*.



2 --- Creating synthetic spectroscopy and photometry
-------------------------------------------------

(Work in progress)

The synthetic SEDs can be created in IDL from the main repository. Before calling the SED-making code [create_FIR_SED](create_FIR_SED.pro) you will need to specify:
* The dust composition 
* The dust temperature distribution
* The redshift range

**Dust composition:** Requires one string vector for the material names and one float vector for the materials' mass fraction. For instance:
~~~IDL
    comp = ['E30R', 'BE']
    compfrac = [.7, .3]
~~~
defines dust made of 70% E30R silicates (from D17B) and 30% BE carbon (from M98), which is the standard dust composition used in the article.

**Temperature distribution:** ...

**Redshift:** ...

~~~IDL
    .r grams_synthphot
    .r physconst
    create_FIR_SED, comp = comp, fcomp = compfrac, T_all = T_array, fT_all = T_frac_array, z_all = z_array, $
                /savesed, /plotsmoothing 
~~~


3a --- Fitting synthetic photometry in the general case (Python)
-------------------------------------------------------------

(Work in progress)

Call the following from shell in the main repository:

~~~Python
    python photometry_fit.py
~~~

When doing the consistency test ('MBBtest'), the script will stop halfway (no need for reduced opacity) with the following message:

~~~Python
    > [Path of the folder containing the Git repository]/Fanciullo_etal_dust_mass_systematics/photometry_fit.py(633)<module>()
    -> currentmode = 'red'
~~~

Just write 'exit' to close.



3b --- Fitting high-redshift, two-band synthetic photometry (IDL)
--------------------------------------------------------------

(Work in progress)

Open IDL in the main repository and write the following:

~~~IDL
    .r grams_synthphot
    .r physconst
    comp = 'E30R-70.0+BE-30.0'
    fit2bands, comp = comp, /saveres
    fit2bands, comp = comp, /red, /saveres
~~~

NOTA: For the standard dust composition in the article --- 'E30R-70.0+BE-30.0', i.e. 70% E30R silicates and 30% BE carbon --- it is not necessary to specify the 'comp' keyword, which is the code's default. One must specify 'comp' for any other composition they wish to fit, though.

The output files are saved in the *fit_result_files/* folder, under the names *2bdfit_{composition}-{raw or red}_oneT_{bands used}.dat*.


4 --- Making the plots
-------------------

Coming soon

