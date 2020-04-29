
This Git repository contains the code and data needed to reporduce the results from "Systematic errors in dust mass determinations: Insights fromlaboratory opacity measurements" (Fanciullo et al. 2020), including figures. 

The scripts perform four different functions: rewriting the original opacity data files into a format that can be used by the rest of the code; use these opacities to create a grid of synthetic SEDs (spectra and photometry); fit the synthetic photometry thus obtained; plot the fit results. The rest of this document explains how to do this step-by-step.

All scripts have to be run from the main directory unless otherwise indicated.


0 --- Software requirements
---------------------------

The code is mostly written in IDL 7.0, with one script in Python 3.6. It has not been tested for other distributions. The following code is needed for the scripts to work:
* The `GRAMS_synthphot` package for IDL, developed by S. Srinivasan, is included in this Git repository.
* `match2.pro`, which can be downloaded off the Internet. This procedure needs to be in your machine's IDL path (`~/idl-libs/`) for [GRAMS_synthphot](grams_synthphot.pro) to work.
* `writecol.pro`, which can be downloaded off the Internet.
* `medabsdev.pro`, a statistics procedure that can be downloaded from [Chris Beaumont's IDL Library](https://www.ifa.hawaii.edu/users/beaumont/code/medabsdev.html). 
* `percentiles.pro`, a statistics procedure that can be downloaded from [www.heliodocs.com](http://www.heliodocs.com/xdoc/xdoc_print.php?file=$SSW/packages/s3drs/idl/util/percentiles.pro).


1 --- Preparation: recovering original opacity and filter profile files
-----------------------------------------------------------------------

**Silicate opacity files:** The silicate opacity from [Demyk et al. 2017a](https://ui.adsabs.harvard.edu/abs/2017A%26A...600A.123D/abstract), [2017b](https://ui.adsabs.harvard.edu/abs/2017A%26A...606A..50D/abstract) (D17A, D17B) can be downloaded from https://www.sshade.eu/db/stopcoda (one needs an SSHADE profile to download the data). Unpack the files in the `MAC_files_original/` directory.

**Carbon opacity files:** The carbon data from [Mennella et al. 1998](https://ui.adsabs.harvard.edu/abs/1998ApJ...496.1058M/abstract) (M98) can be requested to Vito Mennella at the following address: `vito.mennella [at] inaf.it`. They should be unpacked in the `MAC_files_original/M98/` subdirectory.

**Filter profiles:** The `band_profiles/` folder contains files for ALMA bands (in .dat format) and `(GRAMS_synhtphot)[band_profiles/GRAMS_synhtphot.fits]`, which includes PACS, SPIRE and SCUBA2 files that can be used by IDL. However, `GRAMS_synhtphot` cannot be used by the Python script: it is necessary to download additional profiles from the Internet.

* *Herschel bands:* Open the [Herschel page of the SVO Filter Profile Service](http://svo2.cab.inta-csic.es/svo/theory/fps/index.php?mode=browse&gname=Herschel) and download the following files (in ASCII table format) to the `band_profiles/` folder:
    * Herschel_Pacs.blue.dat (PACS 70 um)
    * Herschel_Pacs.green.dat (PACS 100 um)
    * Herschel_Pacs.red.dat (PACS 160 um)
    * Herschel_SPIRE.PSW.dat (SPIRE 250 um)
    * Herschel_SPIRE.PMW.dat (SPIRE 350 um)
    * Herschel_SPIRE.PLW.dat (SPIRE 500 um)
* *SCUBA2 bands:* The SCUBA2 bandpass profiles (model prediction) can be downloaded from [their JCMT page](https://www.eaobservatory.org/jcmt/instrumentation/continuum/scuba-2/filters/) (section "Estimated filter profiles") as `model450.pdf` and `model850.pdf`. Despite their extension, these are text files.
    * Only the first and last column of each file are necessary: save the first and last column of each file in new files named `SCUBA2_450mu.dat` and `SCUBA2_850m.dat`. This procedure can be done automatically by running [rewrite_scuba2_bands](rewrite_scuba2_bands.idl) from within the `band_profiles/` subfolder:
    ~~~IDL
        @rewrite_scuba2_bands.idl
    ~~~
    The files `model850.pdf` and `model450.pdf` can be erased afterwards. 


2 --- Reading and rewriting the opacity files
---------------------------------------------

This section corresponds to section 2.1 of the main article. Open IDL in the main directory and run the following:

~~~IDL
    opacity_reprocess, 'M98', /savedata
    opacity_reprocess, 'D17A', /savedata
    opacity_reprocess, 'D17B', /savedata
~~~

The output files are saved in the `MAC_files_reprocessed/` subdirectory with the filename convention `{material name}_{temperature}K_interpol.xcat`.



3 --- Creating synthetic spectroscopy and photometry
----------------------------------------------------

This corresponds to section 3 of the main article.

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

**Temperature distribution:** Two float arrays, one for the temperatures and one for the (mass) fraction of dust at each temperature. One-dimensional arrays of length *N* are interpreted as a single temperature distribution with *N* different temperature values. It can be more convenient to use two-dimensional arrays of size $N_{T} \times N_{\dist}$, which define $N_{\dist}$ separate distributions of $N_{T}$ values each. For instance:
~~~IDL
    T_array = [[20.], [25.], [30.], [35.], [40.], [45.], [50.], [60.], [80.], [100.]]
    T_frac_array = [[1.], [1.], [1.], [1.], [1.], [1.], [1.], [1.], [1.], [1.]]
~~~
for 10 single-temperature distribution, while:
~~~IDL
    T_array = [[30., 100.], [30., 100.], [30., 100.], [30., 100.], [30., 100.], [30., 100.], [30., 100.], [30., 100.]]
    T_frac_array = [[.9999, .0001], [.9997, .0003], [.999, .001], [.997, .003], [.99, .01], [.97, .03], [.9, .1], [.7, .3]]
~~~
for 8 two-temperature distributions.

**Redshift:** Just a simple array
~~~IDL
    z_array = [0., 1., 2., 3., 4., 5., 6., 7.]
~~~

Other parameters have a default value, but the script also accepts user-defined values. You can see more in `[samplesession.pro](samplesession.pro)`.

Finally, call:
~~~IDL
    .r grams_synthphot
    .r physconst
    create_FIR_SED, comp = comp, fcomp = compfrac, T_all = T_array, fT_all = T_frac_array, z_all = z_array, $
                /savesed, /plotsmoothing 
~~~

The output files are saved in the `synthetic_SEDs/` subdirectory, further subdivided into `synthetic_SEDs/Spectra/` and `synthetic_SEDs/Photometry/`.


4a --- Fitting synthetic photometry in the general case (Python)
----------------------------------------------------------------

(Work in progress)

This corresponds to sections 4.1--4.3 in the main article.

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



4b --- Fitting high-redshift, two-band synthetic photometry (IDL)
-----------------------------------------------------------------

(Work in progress)

This corresponds to section 4.4 in the main article.

Open IDL in the main repository and write the following:

~~~IDL
    .r grams_synthphot
    .r physconst
    comp = 'E30R-70.0+BE-30.0'
    fit2bands, comp = comp, /saveres
    fit2bands, comp = comp, /red, /saveres
~~~

NOTA: For the standard dust composition in the article --- 'E30R-70.0+BE-30.0', i.e. 70% E30R silicates and 30% BE carbon --- it is not necessary to specify the 'comp' keyword, which is the code's default. You must specify 'comp' for any other composition they wish to fit, though.

The output files are saved in the `fit_result_files/` folder, with the naming convention `2bdfit_{composition}-{raw or red}_oneT_{bands used}.dat`.


4 --- Making the plots
----------------------

All plots from the article can be remade using `[makeplots.pro](makeplots.pro)`. Open IDL in the main directory and run the following:
~~~IDL
    makeplots, 1, /saveplot
~~~
after changing the first argument to the number of the figure you want to reproduce. Accepted values are the number from 1 to 15 included (either as integers or as strings, i.e. '1', '2', ...) and the strings 'A1', 'A2', ... to 'A8' (appendix figures).

The plots are saved in the `plots/` folder in EPS format, under filenames starting by `File{number}_`.

A few of the figures are peculiar:
* Fig. 2 and Fig. 9 are composed of two pictures each, so the script saves two `Fig2_...` files and two `Fig9_...` files with different names.
* Since `[makeplots.pro](makeplots.pro)` calls `[create_FIR_SED.pro](create_FIR_SED.pro)` to make Fig. 2, you need to make sure that `grams_synthphot.pro` and `physconst.pro` are compiled before creating the figure:
~~~IDL
    .r grams_synthphot
    .r physconst
    makeplots, 2, /saveplot
~~~


Contact
-------

For any questions, please contact Lapo Fanciullo at `lfanciullo [at] asiaa.sinica.edu.tw`.