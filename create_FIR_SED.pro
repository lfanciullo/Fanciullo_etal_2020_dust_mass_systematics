
FUNCTION make_plaw_tdist, Tparams, gamma

  logTmin = alog10(Tparams[0])
  logTmax = alog10(Tparams[1])
  nT = fix(Tparams[2])
  nfT = n_elements(gamma)
  T = fltarr(nT, nfT)
  f_T = fltarr(nT, nfT)
  T_out = dblarr(nT, nfT, 2)
  for i = 0, nfT-1 do T[*, i] = 10.^(logTmin + (logTmax-logTmin) * findgen(nT)/(nT-1))
  for i = 0, nfT-1 do begin
     f_T[*, i] = (T[*, i]/Tparams[0])^(-gamma[i])
     f_T[*, i] /= total(f_T[*, i])
  endfor
  T_out[*, *, 0] = T & T_out[*, *, 1] = f_T
  
  return, T_out
END


PRO create_FIR_SED, fold_data = fold_data, fold_sed = fold_SEDs, fold_pics = fold_pics, $
                save_sed = save_sed, plotsmoothing = plotsmoothing, plot_SED = plot_SED, save_plot = save_plot, $
                comp = comp, fcomp = fcomp, T_all = T_all, fT_all = fT_all, z_all = z_all, plaw = plaw, filt = filt, $
                T_cmb0 = T_cmb0, M_d = M_d, beta = beta, wl_sed = wl_sed, silent = silent
  
;; NEED TO COMPILE GRAMS_SYNTHPHOT BEFORE USING THIS

;; TCORR: ...
;; SILENT: ...


;;; Initial settings ;;;

if not keyword_set(save_sed) then save_sed = 0
if not keyword_set(plotsmoothing) then plotsmoothing = 0
if not keyword_set(plot_SED) then plot_SED = 0
if not keyword_set(save_plot) then save_plot = 0

if not keyword_set(fold_data) then begin
   fold_data = 'MAC_files_reprocessed/'
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword FOLD_DATA not defined. Defaulting to ', fold_data
endif
if not keyword_set(fold_SEDs) then begin
   fold_SEDs = 'synthetic_SEDs/'
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword FOLD_SEDS not defined. Defaulting to ', fold_SEDs
endif
if not keyword_set(fold_pics) then begin
   fold_pics = 'plots/'
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword FOLD_PICS not defined. Defaulting to ', fold_pics
endif
if not keyword_set(T_cmb0) then begin
   T_cmb0 = 2.725
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword T_CMB0 not defined. Defaulting to ', strtrim(string(T_cmb0), 1)
endif
if not keyword_set(M_d) then begin
   M_d = 1d8  ; Solar masses
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword M_D not defined. Defaulting to ', $
                                           strtrim(string(M_d, format = '(E7.1)'), 1), ' M_sun'
endif
if not keyword_set(comp) then begin
   comp = ['E30R', 'BE']
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword COMP not defined. Defaulting to ', comp
endif
if not keyword_set(fcomp) then begin
   fcomp = [.7, .3]
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword FCOMP not defined. Defaulting to ', $
                                           strtrim(string(fcomp, format = '(F0.2)'), 1)
endif else begin
   if n_elements(fcomp) NE n_elements(comp) then print, 'CREATE_FIR_SED: COMP and FCOMP must have the same number of elements'
endelse
if not keyword_set(T_all) then begin
   T_all = [30.]
   if keyword_set(plaw) then print, 'ERROR (CREATE_FIR_SED): Keyword PLAW cannot be used if T_ALL is not defined. Ignoring it.'
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword T_ALL not defined. Defaulting to ', $
                                           strtrim(string(T_all, format = '(F0.1)'), 1), ' K'
endif
if not keyword_set(fT_all) then begin
   fT_all = [1.]
   if keyword_set(plaw) then print, 'ERROR (CREATE_FIR_SED): Keyword PLAW cannot be used if T_FRAC_ALL is not defined. Ignoring it.'
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword T_FRAC_ALL not defined. Defaulting to ', $
                                           strtrim(string(fT_all, format = '(F0.3)'), 1)
endif
if keyword_set(plaw) and (keyword_set(T_all) and keyword_set(fT_all)) then begin
   T_all_orig = T_all    ; Save original T_all & fT_all because we are going to reuse the variable names
   fT_all_orig = fT_all
   Tdist = make_plaw_tdist(T_all, fT_all)
   T_all = Tdist[*, *, 0]
   fT_all = Tdist[*, *, 1]
endif
if n_elements(fcomp) NE n_elements(comp) then print, 'ERROR (CREATE_FIR_SED): T_ALL and T_FRAC_ALL must have the same number of elements'
if n_elements(z_all) EQ 0 then begin 
   z_all = [1.]
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword Z_ALL not defined. Defaulting to ', $
                                           strtrim(string(z_all, format = '(F0.2)'), 1)
endif
if comp EQ ['MBBtest'] and not keyword_set(beta) then begin
   beta = 1.5
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword BETA not defined. Defaulting to ', $
                                           strtrim(string(beta, format = '(F0.2)'), 1)
endif
if not keyword_set(wl_sed) then begin
   nwl_sed = 450
   wl_sed = 10.^(findgen(nwl_sed)/150 + 1.) ; 10 um - 1 cm, 150 points per order of magnitude
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword WL_SED not defined. Using default value'
endif else begin
   nwl_sed = n_elements(wl_sed)
endelse

if not keyword_set(filt) then begin
   filt = $
      [$
      'PACS70', 'PACS100', 'PACS160', $
      'SPIRE250', 'SPIRE350', 'SPIRE500', $
      'SCUBA2_450', 'SCUBA2_850', $
      'ALMA_10', 'ALMA_9', 'ALMA_8', 'ALMA_7', 'ALMA_6', 'ALMA_5', 'ALMA_4', 'ALMA_3'$
      ]
   if not keyword_set(silent) then print, 'CREATE_FIR_SED: Keyword FILT not defined. Using default value'
endif

nmat = n_elements(comp)
carbon = where(comp EQ 'AC' OR comp EQ 'BE')
sil = where(comp NE 'AC' AND comp NE 'BE')
if (size(T_all))[0] EQ 1 then T_all = reform(T_all, (size(T_all))[1], 1) ;Script needs 2D arrays
if (size(fT_all))[0] EQ 1 then fT_all = reform(fT_all, (size(fT_all))[1], 1)
nT_all = (size(T_all))[2]       ; Number of different T distributions in current run
nT_dist = (size(T_all))[1]      ; Number of individual Ts for each distribution
if (size(T_all))[1] EQ 1 then wlmin = 100. else wlmin = 50.  ; min wl = 100 um for single-T dust, 50 um for multi-T
nz = n_elements(z_all)
nfilt = n_elements(filt)


;; Setting up plots (UPDATE)
loadct, 39
!p.color = 0
tvlct, 125, 125, 125, 1   ; Dark grey
tvlct, 200, 200, 200, 2   ; Light grey
tvlct, 225,   0, 225, 3   ; Fuchsia
tvlct, 255, 182, 193, 4   ; Light pink
tvlct, 150,  75,   0, 5   ; Brown
tvlct, 135, 206, 250, 6   ; Light sky blue
tvlct,   0, 100,   0, 7   ; Dark green
tvlct, 160,   42, 85, 8   ; Sienna
col_all = $
   [0, 250, 75, 1, 7, 3, 0, 250, 75, 1, 7, 3]
if save_plot then begin
   !x.charsize = 1.2
   !x.thick = 6
   !y.charsize = 1.05
   !y.thick = 6
   !p.charthick = 3
   !p.thick = 5
   !p.charsize = 1.1
endif


;;; Synthetic SED & photometry ;;;

;; Construct arrays containing the parameters for each individual SED
nsed = nz * nT_all
T_list = fltarr(nT_dist, nsed)
Tfrac_list = fltarr(nT_dist, nsed)
z_list = fltarr(nsed)
for i = 0, nsed-1 do begin
   T_list[*, i] = T_all[*, [i/nz]]
   Tfrac_list[*, i] = fT_all[*, [i/nz]]
   z_list[i] = z_all[[i MOD nz]]
endfor

;; Create spectra
makesed, wl_sed, sed, T_list, Tfrac_list, comp, fcomp, z_in = z_list, M_d = M_d, filters_in = filt, $
            T_0 = T_cmb0, fold_data = fold_data, plot = plotsmoothing, saveplot = save_plot, silent = silent
;; Create a reduced-opacity version as well (unless comp = MBBtest)
if comp NE ['MBBtest'] then begin
   makesed, wl_sed, sed_reduc, T_list, Tfrac_list, comp, fcomp, z_in = z_list, M_d = M_d, filters_in = filt, /red, $
               T_0 = T_cmb0, fold_data = fold_data, plot = plotsmoothing, saveplot = save_plot, silent = silent
endif

;; Save spectra
if save_sed then begin
   for i = 0, nt_all-1 do begin
      ;; Create output file name
      if comp EQ ['MBBtest'] then begin  ; Composition string
         compstring = 'MBBtest' 
         if (beta GT 0. AND beta NE 1.5) then compstring += '_beta' + $
            strtrim(string(beta, format = '(F0.2)'), 1)  ; Filename includes beta value if it's not the standard one
         compstring += '_'
      endif else begin
         compstring = strjoin(comp + '-' + strtrim(string(100 * fcomp, format = '(F4.1)'), 1), '+')
         compstring_reduc = compstring + '-red_'
         compstring += '-raw_'
      endelse

      CASE nT_dist OF           ; Temperature distribution string
         1: Tstring = 'oneT-' + strtrim(string(T_all[i], format = '(F0.1)'),1) + 'K'
         2: Tstring = 'twoT-' + strjoin(strtrim(string(T_all[*, i], format = '(F0.1)'),1), 'K+') + 'K-fw' + $
                      strtrim(string(fT_all[1, i], format = '(f0.'+strtrim(string(ceil(abs(alog10(fT_all[1, i])))),1) + ')'), 1)
         else: begin
            if keyword_set(plaw) then begin
               Tstring = 'plawT' + strtrim(string(fix(T_all_orig[2])), 1) + '-' + $
                         strjoin(strtrim(string(T_all_orig[0:1], format = '(F0.1)'),1), 'K-') + $
                         'K-ind' + strtrim(string(ft_all_orig[i], format = '(F0.2)')) 
            endif else begin
               nsigfig = ceil(abs(alog10(min(fT_all[*, i]))) -1)
               if nsigfig LT 0 OR not finite(nsigfig) then print, $
                  'ERROR (FIR_SCRIPT_COMPREHENSIVE): degenerate temperature fraction distribution (contains 1 or 0)'
               Tfracformat = '(E0.' + strtrim(nsigfig, 1) + ')'
               Tstring = 'multiT-' + strjoin(strtrim(string(T_all[*, i], format = '(F0.1)'),1), 'K+') + 'K-' + $
                         strjoin(strtrim(string(fT_all[*, i], format = Tfracformat),1), '+')
            endelse
         end
      ENDCASE
      fname = 'Spec_' + compstring + Tstring + '.dat'
;stop
      ;; Writing raw opacity file
      ised = i * nz
      goodspec = where(finite(sed[*, ised]))
      writecol, fold_SEDs + 'Spectra/' + fname, wl_sed[goodspec], sed[goodspec, ised], FMT='(F8.2, "|", E10.3)'
      ;; Writing reduced opacity file
      if comp NE ['MBBtest'] then begin
         fname = 'Spec_' + compstring_reduc + Tstring + '.dat'
         goodspec = where(finite(sed_reduc[*, ised]))
         writecol, fold_SEDs + 'Spectra/' + fname, wl_sed[goodspec], sed_reduc[goodspec, ised], FMT='(F8.2, "|", E10.3)'
      endif
   endfor
endif

;; Create photometry
makesed, wl_sed, photstruct, T_list, Tfrac_list, comp, fcomp, /photometry, /errbars, z_in = z_list, M_d = M_d, $
         filters_in = filt, T_0 = T_cmb0, fold_data = fold_data, plot = plotsmoothing, saveplot = save_plot, silent = silent
wl_phot = photstruct.filt_wl
fl_phot = photstruct.flux
dfl_phot = photstruct.dflux
;; Create a reduced-opacity version as well (unless comp = MBBtest)
if comp NE ['MBBtest'] then begin 
   makesed, wl_sed, photstruct_reduc, T_list, Tfrac_list, comp, fcomp, /photometry, /errbars, z_in = z_list, M_d = M_d, /red, $
            filters_in = filt, T_0 = T_cmb0, fold_data = fold_data, plot = plotsmoothing, saveplot = save_plot, silent = silent
   fl_phot_reduc = photstruct_reduc.flux
   dfl_phot_reduc = photstruct_reduc.dflux
endif

;; Save photometry
if save_sed then begin
   for i = 0, nsed-1 do begin
      ;; Create output file name
      if comp EQ ['MBBtest'] then begin  ; Composition string
         compstring = 'MBBtest' 
         if (beta GT 0. AND beta NE 1.5) then compstring += '_beta' + strtrim(string(beta, format = '(F0.2)'), 1)
         compstring += '_'
      endif else begin
         compstring = strjoin(comp + '-' + strtrim(string(100 * fcomp, format = '(F4.1)'), 1), '+')
         compstring_reduc = compstring + '-red_'
         compstring += '-raw_'
      endelse
      CASE nT_dist OF           ; Temperature distribution string
         1: Tstring = 'oneT-' + strtrim(string(T_list[i], format = '(F0.1)'),1) + 'K'
         2: Tstring = 'twoT-' + strjoin(strtrim(string(T_list[*, i], format = '(F0.1)'),1), 'K+') + 'K-fw' + $
                      strtrim(string(fT_all[1, i/nz], format = '(f0.'+strtrim(string(ceil(abs(alog10(fT_all[1, i/nz])))),1) + ')'), 1)
         else: begin
            if keyword_set(plaw) then begin
               iT = i mod nT_all
               Tstring = 'plawT' + strtrim(string(fix(T_all_orig[2])), 1) + '-' + $
                         strjoin(strtrim(string(T_all_orig[0:1], format = '(F0.1)'),1), 'K-') + $
                         'K-ind' + strtrim(string(ft_all_orig[iT], format = '(F0.2)'))
            endif else begin
               nsigfig = ceil(abs(alog10(min(fT_all[*, i/nz]))) -1)
               if nsigfig LT 0 OR not finite(nsigfig) then print, $
                  'ERROR (FIR_SCRIPT_COMPREHENSIVE): degenerate fT distribution (contains 1 or 0)'
               Tfracformat = '(E0.' + strtrim(nsigfig, 1) + ')'
               Tstring = 'multiT-' + strjoin(strtrim(string(T_all[*, i], format = '(F0.1)'),1), 'K+') + 'K-' + $
                         strjoin(strtrim(string(fT_all[*, i/nz], format = Tfracformat),1), '+')
            endelse
         end
      ENDCASE
      zstring = '_z' + strtrim(string(z_list[i], format = '(F4.2)'), 1) ; Redshift string
      fname = 'Phot_' + compstring + Tstring + zstring

      ;stop
      
      ;; RAW OPACITY FILE
      ;; Band selection
      goodbd_bool_wl = wl_phot GE wlmin * (1 + z_list[i])  ; Min wavelength cut
      goodbd_bool_finite = finite(fl_phot[*, i])           ; Finite values only
      goodbd_bool_snr = (fl_phot/dfl_phot)[*, i] GE 3.     ; Require decent S/N.
      goodbd_bool_bdchoice = 1 - strmatch(filt, 'SCUBA2*') ; Final band selection to avoid redundancy
      if (fl_phot/dfl_phot)[where(filt EQ 'SPIRE350'), i] GE 3. then goodbd_bool_bdchoice *= 1 - strmatch(filt, 'ALMA_10')
      if (fl_phot/dfl_phot)[where(filt EQ 'SPIRE500'), i] GE 3. then goodbd_bool_bdchoice *= 1 - strmatch(filt, 'ALMA_9')

      goodbd_strict = where(goodbd_bool_finite * goodbd_bool_wl * $
                            goodbd_bool_snr * goodbd_bool_bdchoice)

      ;; Write output file (raw opacity)
      if goodbd_strict EQ [-1] then begin
         print, 'WARNING: No valid data for T = ' + strtrim(string(T_list[*, i], format = '(F0.2)'), 1) + ', z = ' + $
                strtrim(string(z_list[i], format = '(F0.2)'), 1) + '. No file recorded.'
      endif else begin
         nlines_temp = n_elements(goodbd_strict)
         openw, 1, fold_SEDs + 'Photometry/' + fname + '.dat'  
         printf, 1, '#   Filter |   wl (um) | Flux (Jy) |   dFlux   |'
         printf, 1, '# ----------------------------------------------'
         printform = '(A11, "|", F11.2, "|", E11.3, "|", E11.3)'
         for j = 0, nlines_temp-1 do printf, 1, filt[goodbd_strict[j]], wl_phot[goodbd_strict[j]], fl_phot[goodbd_strict[j], i], $
                                        dfl_phot[goodbd_strict[j], i], format = printform
         close, 1
      endelse

      ;; REDUCED OPACITY FILE
      if comp NE ['MBBtest'] then begin
         fname_reduc = 'Phot_' + compstring_reduc + Tstring + zstring
         ;; Band selection
         goodbd_bool_wl = wl_phot GE wlmin * (1 + z_list[i])          ; Min wavelength cut
         goodbd_bool_finite = finite(fl_phot_reduc[*, i])             ; Finite values only
         goodbd_bool_snr = (fl_phot_reduc/dfl_phot_reduc)[*, i] GE 3. ; Require decent S/N.
         goodbd_bool_bdchoice = 1 - strmatch(filt, 'SCUBA2*')         ; Final band selection to avoid redundancy
         if (fl_phot_reduc/dfl_phot_reduc)[where(filt EQ 'SPIRE350'), i] GE 3. then goodbd_bool_bdchoice *= 1 - strmatch(filt, 'ALMA_10')
         if (fl_phot_reduc/dfl_phot_reduc)[where(filt EQ 'SPIRE500'), i] GE 3. then goodbd_bool_bdchoice *= 1 - strmatch(filt, 'ALMA_9')
         
         goodbd_strict = where(goodbd_bool_finite * goodbd_bool_wl * $
                               goodbd_bool_snr * goodbd_bool_bdchoice)

         ;; Write output file
         if n_elements(goodbd_strict) LT 4 then begin
            print, 'WARNING: Not enough valid bands for T = ' + strtrim(string(T_list[*, i], format = '(F0.2)'), 1) + ', z = ' + $
                   strtrim(string(z_list[i], format = '(F0.2)'), 1) + '. No file recorded.'
         endif else begin
            nlines_temp = n_elements(goodbd_strict)
            openw, 1, fold_SEDs + 'Photometry/' + fname_reduc + '.dat'  
            printf, 1, '#   Filter |   wl (um) | Flux (Jy) |   dFlux   |'
            printf, 1, '# ----------------------------------------------'
            printform = '(A11, "|", F11.2, "|", E11.3, "|", E11.3)'
            for j = 0, nlines_temp-1 do printf, 1, filt[goodbd_strict[j]], wl_phot[goodbd_strict[j]], fl_phot_reduc[goodbd_strict[j], i], $
                                                dfl_phot_reduc[goodbd_strict[j], i], format = printform
            close, 1
         endelse
      
      endif
        
   endfor

endif


END
