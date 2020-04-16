
FUNCTION uncphot, flux, filters
  
nbands = (size(flux))[1] ;n_elements(filters)
nsed = (size(flux))[2]
if (size(flux))[0] EQ 2 then nmat = 1 else nmat = (size(flux))[3]
nsample = 10000
err_samples = fltarr(nbands, nsed, nmat, nsample)


;; Opacity uncertainty
unc_oprop = .10 

;; Photon noise (ALMA only)
unc_gamma =  [      .10,      .10,      .10,      .10,      .10,      .10,      .10,      .10]
filt_gamma = ['ALMA_10', 'ALMA_9', 'ALMA_8', 'ALMA_7', 'ALMA_6', 'ALMA_5', 'ALMA_4', 'ALMA_3']
n_unc_gamma = n_elements(unc_gamma)

;; Calibration uncertainty
unc_calib =   [   .05,    .055,          .12,          .08]
filt_calibs = ['PACS', 'SPIRE', 'SCUBA2_450', 'SCUBA2_850']
n_unc_calib = n_elements(unc_calib)
;; PACS: Poglitsch et al. '10
;; SPIRE: SPIRE handbook
;; SCUBA2: Dempsey et al. '13

;; Confusion noise (not for ALMA)
unc_conf =  [      .1,        .3,        1.,         6.,         6.,         7.,           .5,           .7]*1e-3 ;In Jy/beam
filt_conf = ['PACS70', 'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500', 'SCUBA2_450', 'SCUBA2_850']
n_unc_conf = n_elements(unc_conf)

;; Roll of the dice
dice_oprop = randomn(seed, nsed, nmat, nsample)
dice_gamma = randomn(seed, n_unc_gamma, nsed, nsample)
dice_calib = randomn(seed, n_unc_calib, nsed, nsample)
dice_conf = randomn(seed, n_unc_conf, nsed, nsample)
samples = fltarr(nbands, nsed, nmat, nsample)

;; Optical properties systematics (same effect on all bands)
for i = 0, nbands-1 do begin
   for j = 0, nsed-1 do begin
      for k = 0, nmat-1 do err_samples[i, j, k, *] = flux[i, j, k] * (1 + unc_oprop * dice_oprop[j, k, *])
   endfor
endfor
;; Sum on materials
if nmat GT 1 then err_samples = total(err_samples, 3) else err_samples = reform(err_samples)
;; ; Test example
;; stop
;; print, flux[*, 0, 0], flux[*, 0, 1]
;; print, reform(1 + unc_oprop * dice_oprop[0, 0, 0:5]), reform(1 + unc_oprop * dice_oprop[0, 1, 0:5])
;; print, reform(err_samples[*, 0, 0:2])

;; Photon noise
for i = 0, n_unc_gamma-1 do begin
   indtemp_all = where(strmatch(filters, '*' + filt_gamma[i] + '*') EQ 1, nband_temp) ; Select relevant bands
   if nband_temp GT 0 then begin
      for j = 0, nband_temp-1 do begin
         indtemp = indtemp_all[j] ; Single band
         for k = 0, nsed-1 do err_samples[indtemp, k, *] *= reform(1 + unc_gamma[i] * dice_gamma[i, k, *])
      endfor
   endif
endfor

;; Calibration systematics
for i = 0, n_unc_calib-1 do begin
   ;; NOTA: This is the calibration uncertainty and therefore the same
   ;; (for each "roll") on all bands of the instrument. 
   indtemp_all = where(strmatch(filters, '*'+filt_calibs[i]+'*') EQ 1, nband_temp)
   if nband_temp GT 0 then begin
      for j = 0, nband_temp-1 do begin
         indtemp = indtemp_all[j] ; Single band
         for k = 0, nsed-1 do err_samples[indtemp, k, *] *= reform(1 + unc_calib[i] * dice_calib[i, k, *])
      endfor
   endif
endfor

;; Confusion noise for Herschel + SCUBA2 bands
for i = 0, n_unc_conf-1 do begin
   indtemp_all = where(strmatch(filters, filt_conf[i]) EQ 1, nband_temp)
   if nband_temp GT 0 then begin
      for j = 0, nband_temp-1 do begin
         indtemp = indtemp_all[j] ; Single band
         for k = 0, nsed-1 do err_samples[indtemp, k, *] += unc_conf[i] * dice_conf[i, k, *]
      endfor
   endif
endfor

;; Extracting median + error bars from distribution
;res = fltarr(nbands, 2)
;res[*, 0] = median(err_samples, dimension = 2)
;for i = 0, nbands-1 do res[i, 1] = stddev(err_samples[i,*])
statarr = fltarr(nbands, nsed, 3)
for i = 0, nbands-1 do begin
   for k = 0, nsed-1 do begin
      statarr[i, k, *] = percentiles(err_samples[i, k, *], value = [.16, .5, .84])
   endfor
endfor

phot_median = statarr[*, *, 1]
phot_err_plus = statarr[*, *, 2] - statarr[*, *, 1]  ; Positive error bars
phot_err_minus = statarr[*, *, 1] - statarr[*, *, 0] ; Negative error bars
err_asymm = (phot_err_plus - phot_err_minus) / (phot_err_plus + phot_err_minus)
if max(abs(err_asymm)) GT .1 then print, 'UNCPHOT.PRO WARNING: Some error bars are asymmetric (max asymm = ' + $
                                         strtrim(string(max(abs(err_asymm)), format = 'F0.2'), 1)+ '). We suggest looking into it.'
phot_err = (phot_err_plus + phot_err_minus)/2

;stop

return, [[[phot_median]], [[phot_err]]]

END
