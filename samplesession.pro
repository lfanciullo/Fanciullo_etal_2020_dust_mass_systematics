
PRO samplesession


;;; REPROCESSING OPACITY ;;;

opacity_reprocess, 'M98', /savedata
opacity_reprocess, 'D17A', /savedata
opacity_reprocess, 'D17B', /savedata



;;; CREATING SYNTHETIC SPECTROSCOPY/PHOTOMETRY ;;;

;; Settings
savesed = 1        ; Save the synthetic spectra + photometry?
plot_smooth = 0   ; Plot the effect of opacity smoothing?
plotSED = 0         ; Plot the SED (spectrum + photometry)?
saveplot = 0       ; Save the plots?

;; Parameters
folddata = 'MAC_files_reprocessed/'  ; Folder containing the (reprocessed) opacity
foldSEDs = 'synthetic_SEDs/'         ; Folder for the synthetic spectra + photometry
foldpics = 'plots/'                  ; Plots folder
Tcmb0 = 2.725  ; K
Md = 1d8       ; Solar masses
beta = 1.5      ; Beta value for MBBtest 'material'
wlsed = 10.^(findgen(450)/150 + 1.) ; Wavelength range: 10 um - 1 cm, 450 values (150 per order of magnitude)
filters = ['PACS70', 'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500', 'SCUBA2_450', 'SCUBA2_850', $
           'ALMA_10', 'ALMA_9', 'ALMA_8', 'ALMA_7', 'ALMA_6', 'ALMA_5', 'ALMA_4', 'ALMA_3']

;; Dust model characteristics (comment/uncomment/modify as needed):
comp = $
   ;['E30R', 'BE']
   ['MBBtest']
compfrac = $
   ;[.7, .3]
   [1.]
T_array = $
   ;[[20.], [25.], [30.], [35.], [40.], [45.], [50.], [60.], [80.], [100.]]
   [[30., 100.], [30., 100.], [30., 100.], [30., 100.], [30., 100.], [30., 100.], [30., 100.], [30., 100.]]
T_frac_array = $
   ;[1.] & Tdist = 'singleT'  ; One single-temperature model
   ;[[1.], [1.], [1.], [1.], [1.], [1.], [1.], [1.], [1.], [1.]] & Tdist = 'singleT'  ; Multiple single-temperature models
   ;[.9, .1] & Tdist = 'customT'  ; One two-temperatures model
   [[.9999, .0001], [.9997, .0003], [.999, .001], [.997, .003], [.99, .01], [.97, .03], [.9, .1], [.7, .3]] ;& Tdist = 'customT'
z_array = [0., 1., 2., 3., 4., 5., 6., 7.]

;; Launching the script
.r grams_synthphot
.r physconst
create_FIR_SED, fold_data = folddata, fold_sed = foldseds, fold_pics = foldpics, $
                save_sed = savesed, plotsmoothing = plot_smooth, plot_SED = plotSED, save_plot = saveplot, $
                comp = comp, fcomp = compfrac, T_all = T_array, fT_all = T_frac_array, z_all = z_array, $
                filt = filters, T_cmb0 = Tcmb0, M_d = Md, beta = beta, wl_sed = wlsed



;;; 2-BAND PHOTOMETRY FIT ;;;
;; Nota: all-band fit is done with Python

beta_2bd = [1.5, 1.75, 2.]

fit2bands, ...
fit2bands, ..., /red


;;; PLOTS ;;;

;; Work in progress




END
