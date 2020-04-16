FUNCTION makesynthphot, wl_in, sed_in, filters_in = filters_in

; Input:
;  WL_IN: wavelenght (1-D array)
;  SED_IN: Input SEDs (1- or 2_D array). If 2D array (multiple SEDs)
;          they all share the same wavelength (wl_in)
;  FILTERS_IN: (optional) filters array of length N
; Output:
;  2-D photometry of size Nx3, where:
;   N = # of filters
;   First row is central filter wavelengths
;   2nd and 3rd row are (CMB-corrected) photometry and uncertainty, respectively


;; If filters are undefined: default filters
if not keyword_set(filters_in) then begin
   filters_in = ['PACS100', 'SPIRE250', 'SPIRE350', 'SPIRE500', 'SCUBA2_850']
   print, 'MAKESYNTHPHOT: Filters not defined. Defaulting to:', filters_in
endif
filters_in = strupcase(filters_in)

;; Find value of relevant parameters
n_filt = n_elements(filters_in)
n_wl = n_elements(wl_in)
;n_sed = (size(sed_in))[2]
CASE (size(sed_in))[0] OF
   1: begin
      n_sed = 1
      n_mat = 1
   end
   2: begin
      n_sed = (size(sed_in))[2]
      n_mat = 1
   end
   3: begin
      n_sed = (size(sed_in))[2]
      n_mat = (size(sed_in))[3]
   end
   else: print, 'ERROR (MAKESYNTHPHOT): Array SED_IN is ', strtrim(string((size(sed_in))[0]),1), '-D. Only 1- to 3-D arrays accepted.'
   ENDCASE
;if (size(sed_in))[0] EQ 2 then n_mat = 1 else n_mat = (size(sed_in))[1]
;stop
sed2use = reform(sed_in, n_wl, n_sed, n_mat)

;; Separate photometry into ALMA ('red') & non-ALMA ('blue') sections
;;  non-ALMA to be calculated with synthphot (Sundar's code);
;;  ALMA with my own code
filters_blue = WHERE(STRMATCH(filters_in, 'ALMA_*') EQ 0) ; Non-ALMA bands
filters_red = WHERE(STRMATCH(filters_in, 'ALMA_*') EQ 1)  ; ALMA bands
wl_out = fltarr(n_filt)
fphot_out = fltarr(n_filt, n_sed, n_mat)
fphot_out = reform(fphot_out, n_filt, n_sed, n_mat) ; Don't allow IDL to erase (eventual) shallow dimensions
wl2use = reform(wl_in, n_wl, 1) ; Synthphot needs 2-D arrays
;sed2use = reform(sed_in)

;stop

;; Get non-ALMA photometry with synthphot
if filters_blue NE [-1] then begin
   ;; Select only the filters you're interested in
   filters_blue_all = mrdfits('GRAMS_filters.fits', 1)
   match2, strtrim(filters_in, 2), strtrim(filters_blue_all.filter_name, 2), a, b, na, nb
   ;; a, b = matches; na, nb = non-matches. The filters wanted are filters_blue_all[b].filter_name or filters_in[a]
   if (size(a))[0] NE 0 then wl_out[filters_blue] = filters_blue_all[b].lamref
   for i = 0, n_mat-1 do begin
      synthphot, wl_in, transpose(sed2use[*, *, i]), fphot_blue, fitspath = './'
      fphot_blue = transpose(fphot_blue)
      if (size(a))[0] NE 0 then fphot_out[filters_blue, *, i] = fphot_blue[b, *]
   endfor
endif

;stop

;; Get ALMA photometry
if filters_red NE [-1] then begin
   for i = 0, n_mat-1 do begin
      synthphot_alma, wl_in, sed2use[*, *, i], filters_in, wl_red, photom_red, /silent
      wl_out[filters_red] = wl_red  ; A bit redundant, but it works
      fphot_out[filters_red, *, i] = photom_red
   endfor
endif

;stop
;;st1 = {filt_name:filters_in}
;st2 = {filt_wl:0.}
;if n_mat EQ 1 then fphot_1filt = fltarr(n_sed) else fphot_1filt = fltarr(n_sed, n_mat)
;st3 = {filt_phot:fphot_1filt}

res = {filt_name:filters_in, filt_wl:wl_out, n_mat:n_mat, fluxes:fphot_out}

;res = [[wl_out], [fphot_out]]

return, res

end
