
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

This will populate the MAC_files_reprocessed/ subdirectory with the files needed to create the synthetic SEDs.



2 --- Creating synthetic spectroscopy and photometry
-------------------------------------------------

Coming soon



3a --- Fitting synthetic photometry in the general case (Python)
-------------------------------------------------------------

Coming soon



3b --- Fitting high-redshift, two-band synthetic photometry (IDL)
--------------------------------------------------------------

Coming soon



4 --- Making the plots
-------------------

Coming soon

