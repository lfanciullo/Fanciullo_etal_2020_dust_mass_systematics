
PRO makesed, wl_sed, $
             sed_out, $
             T_in, $
             fT_in, $
             comp, $
             f_comp, $
             photometry = photometry, $
             beta = beta, $
             T_0 = T_0, $
             fold_data = fold_data, $
             red = red, $
             errbars = errbars, $
             z_in = z_in, $
             M_d = M_d, $
             filters_in = filters_in, $
             plot = plot, $
             saveplot = saveplot, $
             nocmbcorr = nocmbcorr, $
             T_split = T_split, $
             comp_split = comp_split, $
             silent = silent
                
;; ADD DESCRIPTION OF PROCEDURE
;; sed_out
;; T_in
;; fT_in
;; Ttype
;; comp
;; f_comp
;; photometry
;; beta
;; T_0
;; fold_data
;; red
;; errbars
;; z_in
;; M_d
;; filters_in
;; plot
;; saveplot
;; nocmbcorr
;; T_split
;; comp_split
;; silent
  
@makesed_params.idl
  
nwl_sed = n_elements(wl_sed)
nsed = n_elements(z_in)
nfT = (size(fT_in))[1]
nmat = n_elements(comp)
op_all = dblarr(nwl_sed, nsed, nfT, nmat)
sed_split = dblarr(nwl_sed, nsed, nfT, nmat)
if not keyword_set(T_0) then begin
   T_0 = 2.725
   if not keyword_set(silent) then print, 'MAKESED: Keyword T_0 not specified. Defaulting to ', $
                                          strtrim(string(T_0, format = '(F0.3)'), 1), 'K'
endif
if not keyword_set(M_d) then begin
   M_d = 1d8                    ; In Msun
   if not keyword_set(silent) then print, 'MAKESED: Keyword M_D not defined. Defaulting to ', $
                                           strtrim(string(M_d, format = '(E7.1)'), 1), ' M_sun'
endif
M_d_cgs = M_d * 1.9885d33       ; Conversion M_sun --> g


for i = 0, nsed-1 do begin
   z = z_in[i]
   T = T_in[*, i]
   fT = fT_in[*, i]
   if z NE 0 then wl_eff = wl_sed/(1 + z) else wl_eff = wl_sed

   ;; Multiplicative factor, assuming MAC opacity and point source
   ;; (Emission *4 pi to get total emission (assume isotropy), /4 pi d^2 to get flux)
   if z NE 0 then dist = lumdist(z, /silent) else dist = d_0  ; Distance in Mpc
   dist_cgs = dist * 3.0857d24                                ; Conversion to cm
   em_factor = M_d_cgs * (1+z) / dist_cgs^2
   
   ;; Material opacity
   if (comp NE ['MBBtest']) then begin
      ;; Opacity folder path
      if not keyword_set(fold_data) then begin
         fold_data = 'MAC_files_reprocessed/'
         if not keyword_set(silent) then print, 'MAKESED: FIR opacity directory not defined. Defaulting to ', fold_data
      endif

      ;; Comp keywords
      op_all[*, i, *, *] = make_opac(fold_data, comp, T, wl_eff, plot = plot, saveplot = saveplot)
      nmat = n_elements(comp)
      if keyword_set(beta) and not keyword_set(silent) then print, $
         'MAKESED: Specified keyword BETA is not needed when using experimental MAC measurements'
   endif else begin
      f_comp = dblarr(nmat)+1.
      if keyword_set(beta) then sig_data[2] = -beta else begin
         if not keyword_set(silent) then print, 'MAKESED: No value specified for BETA. Defaulting to BETA = ', $
                                                string(-sig_data[2], format = '(F0.2)')
      endelse
      op_all_1T = sig_data[0] * (wl_eff/sig_data[1])^(sig_data[2])
      for j = 0, nfT-1 do op_all[*, i, j] = op_all_1T
   endelse

   ;; Reduced dust opacity case
   if keyword_set(red) then begin
      carbon = where(comp EQ 'AC' OR comp EQ 'BE')
      sil = where(comp NE 'AC' AND comp NE 'BE')
      if carbon NE [-1] then op_all[*, i, *, carbon] /= redfact_c
      if sil NE [-1] then op_all[*, i, *, sil] /= redfact_sil
   endif

   ;; SED creation
   ;  ss_bbfunct.pro outputs in W m^-2 Hz^-1 (default) or W m^-2 Hz^-1 sr^-1 (/INTENSITY)
   ;  1 Jy = 10^-23 erg s^-1 cm^-2 Hz^-1 = 10^-26 W m^-2 Hz^-1
   bb = dblarr(nwl_sed, nfT)
   for j = 0, nfT-1 do begin
      bb[*, j] = ss_bbfunc(wl_eff, T[j], /intensity) * 1e3 ; mks --> cgs conversion
      if not keyword_set(nocmbcorr) then bb[*, j] -= ss_bbfunc(wl_eff, T_0 * (1+z), /intensity) * 1e3
      ;; CMB background subtraction (da Cunha et al. 2013)
      for k = 0, nmat-1 do sed_split[*, i, j, k] = em_factor * op_all[*, i, j, k] * fT[j] * bb[*, j] * f_comp[k]
   endfor
   sed_split[*, i, *, *] *= 1e23  ; Conversion --> Jy (point source) or Jy/sr (extended)
endfor

;; SEDs split by temperature
sed_split = reform(sed_split, nwl_sed, nsed, nfT, nmat) ; Don't allow IDL to erase (eventual) shallow dimensions
sed_multi_T = total(sed_split, 4) 

;; SEDs split by comp
sed_split = reform(sed_split, nwl_sed, nsed, nfT, nmat) ; Don't allow IDL to erase (eventual) shallow dimensions
sed_multi_mat =  total(sed_split, 3)

;; Overall SEDs
sed_multi_T = reform(sed_multi_T, nwl_sed, nsed, nfT) ; Don't allow IDL to erase (eventual) shallow dimensions
sed =  total(sed_multi_T, 3)

;; Photometry
if keyword_set(photometry) then begin
   if keyword_set(errbars) then begin
      str_phot_multi_mat = makesynthphot(wl_sed, sed_multi_mat, filters_in = filters_in) ; Need multi-material version for error bars
      phot_and_err = uncphot(str_phot_multi_mat.fluxes, filters_in)   
      photometry = {filt_name:filters_in, filt_wl:str_phot_multi_mat.filt_wl, flux:phot_and_err[*, *, 0], dflux:phot_and_err[*, *, 1]}
   endif else begin
      str_phot = makesynthphot(wl_sed, sed, filters_in = filters_in)
      photometry = {filt_name:filters_in, filt_wl:str_phot.filt_wl, flux:str_phot.fluxes}
   endelse
   sed_out = photometry
endif else begin
   sed_out = sed
endelse


END
