
;; This sets up everything you need to run the IDL scripts
;; WIP

;; TO INTEGRATE:
;nwl_sed = 450
;wl_fineness = 150. ;Points per order of magnitude in wl
;wl_sed = 10.^(findgen(nwl_sed)/wl_fineness + 1.) ; 10 um - 1 cm

;; (The following in a case/of loop)
;material = ...
;beta = ...
;material_frac = ...
;nmat = n_elements(material)
;carbon = where(material EQ 'AC' OR material EQ 'BE')
;sil = where(material NE 'AC' AND material NE 'BE')
;T_all = ...
;if (size(T_all))[0] EQ 1 then T_all = reform(T_all, (size(T_all))[1], 1) ;Script needs 2D arrays
;T_frac_all = ... & Tdist = ...
;if (size(T_frac_all))[0] EQ 1 then T_frac_all = reform(T_frac_all, (size(T_frac_all))[1], 1)
;n_Tall = (size(T_all))[2] ; Number of different T distributions in current run
;if Tdist EQ 'singleT' then wlmin = 100. else wlmin = 50.


   
;; Compile scripts needed for SED synthesis and fitting
.r grams_synthphot
.r physconst

;; Set up directories
folder_master = $   ; Path to the main folder on your machine. Modify as needed.
     '~/Documents/Postdoc/Science/IDL_stuff/Fanciullo_etal_dust_mass_systematics/'
folder_kappafiles_orig = folder_master + 'MAC_files_original/'
folder_kappa = folder_master + 'MAC_files_regrid/'
folder_SED = folder_master + 'synthetic_SEDs/'
folder_fits = folder_master + 'fit_result_files/'
folder_plots = folder_master + 'plots/'

;; Parameters:
;; Instrumental
filt = [$
       'PACS70', 'PACS100', 'PACS160', $
       'SPIRE250', 'SPIRE350', 'SPIRE500', $
       'SCUBA2_450', 'SCUBA2_850', $
       'ALMA_10', 'ALMA_9', 'ALMA_8', 'ALMA_7', 'ALMA_6', 'ALMA_5', 'ALMA_4', 'ALMA_3'$
       ]
;; Apparently wrong by up to 50%
;fwhm = [$                      ; In arcsec
;       5.67, 7.04, 11.18, $    ; Aniano et al. 2011
;       18.15, 24.88, 36.09, $  ; Ditto
;       7.9, 13., $             ; Dempsey et al. 2013
;       0., 0., 0., 0., 0., 0., 0., 0. $  ; ALMA is a special case
;       ]
;beam_surf = fwhm^2 * !pi/(4. * alog(2.))
beam_surf = [$                   ; In arcsec^2
            [5.67, 7.04, 11.18]^2 * !pi/(4.*alog(2.)), $ ; PACS (Aniano+11 & Dempsey+13)
            469., 831., 1804., $ ; SPIRE from http://herschel.esac.esa.int/twiki/pub/Public/SpireDocsEditableTable/SPIRE_at_a_glance_v1.pdf
            104., 228., $        ; Dempsey+13
            1., 1., 1., 1., 1., 1., 1., 1. $ ; ALMA (arbitrary)
            ]
;; Wavelength
nwl_sed = 450
wl_fineness = 150. ;Points per order of magnitude in wl
wl_sed = 10.^(findgen(nwl_sed)/wl_fineness + 1.) ; 10 um - 1 cm
;; Other
Tmodel = [20., 25., 30., 35., 40., 45., 50., 60., 80., 100.]
nTmodel = n_elements(Tmodel)
zmodel = [0., 1., 3., 5., 7.]
nzmodel = n_elements(zmodel)

;; Set up plots
; ... ?








			
