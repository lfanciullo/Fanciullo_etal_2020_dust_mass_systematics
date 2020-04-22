
FUNCTION parfitter, x, y, x2plot = x2plot, y2plot = y2plot

;; Fits a parabolic function to the input data, returns x of the
;; minumum and half-width at y = min(y) + 1.
;; Optional outputs: parabola to plot over data
  
x = reform(x) & y = reform(y)
coeff = poly_fit(x, y, 2)
a = coeff[2]  ; Parabola coefficients: y = a * x^2 + b * x + c
b = coeff[1]
c = coeff[0]
fit_res_best = - b/(2*a)
fit_res_unc = 1./sqrt(a)
res = [fit_res_best, fit_res_unc]
return, res

END



PRO fit2bands, bands = bands, comp = comp, Treal = Treal, zmod = zmod, beta = beta, red = red, silent = silent, saveres = saveres
  
  
;; Fits a 2-band photometry using the same technique as Knudsen et
;; al. 2017, Hashimoto et al. 2019
;;
;; Needs to compile grams_synthphot, physconst
;; Usage example: fit2bands, beta = [1.5, 2., 2.5], bands = ['ALMA_6', 'ALMA_3'], comp = 'FAYA-70.0+AC-30.0', /red, /saveres
;;
;; BETA: ...
;; BANDS: ...
;; COMP: ...
;; SAVERES: ...
;; RED: ...
;; SILENT: ...
  
;; Presets
if not keyword_set(beta) then begin
   beta = [1.5, 1.75, 2.]
   if not keyword_set(silent) then print, 'FIT2BANDS: BETA not set. Defaulting to ' + strjoin(strsplit(strjoin(string(beta)), /extract), ', ')
endif
if not keyword_set(bands) then begin
   bands = $
      [['ALMA_7', 'ALMA_6'], $
       ['ALMA_6', 'ALMA_3'], $
       ['ALMA_7', 'ALMA_3']]
   if not keyword_set(silent) then begin
      nsets = (size(bands))[2]
      print, 'BANDS not set. Defulting to the following sets:'
      for i = 0, nsets-1 do print, strjoin(strsplit(strjoin(bands[*, i]), /extract), ', ')
      ;; band:   ALMA_10,  _9,  _8,  _7,   _6,   _5,   _4,   _3
      ;; wl (um):    345, 442, 744, 873, 1287, 1477, 2068, 3075
   endif
endif
if (size(bands))[1] EQ 1 then bands = reform(bands, (size(bands))[1], (size(bands))[0])
nbands = (size(bands))[1]
nsets = (size(bands))[2]
if nbands NE 2 then print, 'ERROR (FIT2BANDS.PRO): EACH BAND SET MUST HAVE 2 BANDS'
if not keyword_set(comp) then begin
   comp = $
      'E30R-70.0+BE-30.0'
      ;;'MBBtest'
   if not keyword_set(silent) then print, 'FIT2BANDS: Variable COMP not defined. Defaulting to ', comp
endif
if not keyword_set(Treal) then begin
   Treal = [20., 25., 30., 35., 40., 45., 50., 60., 80., 100.]
   if not keyword_set(silent) then print, 'FIT2BANDS: Variable TREAL_ALL not defined. Defaulting to: ', Treal
endif
if not keyword_set(zmod) then begin
   zmod = [5., 6., 7.]
   if not keyword_set(silent) then print, 'FIT2BANDS: Variable ZMOD not defined. Defaulting to ', zmod
endif
if not keyword_set(T_CMB) then begin
   T_CMB = 2.725
   if not keyword_set(silent) then print, 'FIT2BANDS: Variable T_CMB not defined. Defaulting to ', strtrim(string(T_CMB, format = '(F0.3)'), 1)
endif

nwl = 450 & wl = 10.^(dindgen(nwl)/150. + 1.) ; Wavelength 10 um - 1 cm
nTreal = n_elements(Treal)
nzmod = n_elements(zmod)
nbeta = n_elements(beta)
;zmod = [5., 6., 7.]
nzmod = n_elements(zmod)
conv_factor = 1d26  ; Conversion factor to Jy for ss_bbfunc output
kappa_0 = .7        ; Standard opacity (James et al. 2002)
wl_0 = 850.         ; Standard opacity (James et al. 2002)
Mgrid = 10.^(7 + findgen(201)/50) & nMgrid = n_elements(Mgrid)


;; Defining folders used
folder_phot = 'synthetic_SEDs/Photometry/'
folder_fit = 'fit_result_files/'

;; Input file name (partial)
fname_start = 'Phot_' + comp
if comp NE 'MBBtest' then begin
   if keyword_set(red) then fname_start += '-red' else fname_start += '-raw'
endif
fname_start += '_oneT-'


;; Let's start!
for iset = 0, nsets-1 do begin
   bands_temp = bands[*, iset]
   if keyword_set(saveres) then begin
      fname_out = '2bdfit_' + comp
      if comp NE 'MBBtest' then begin
         if keyword_set(red) then fname_out += '-red_oneT_' else fname_out += '-raw_oneT_'
      endif else begin
         fname_out += '_oneT_'
      endelse
      fname_out += bands_temp[0] + '+' + bands_temp[1] + '.dat'
      openw, 1, folder_fit + fname_out
      printf, 1, '# Two-band fit (min. chi^2) for beta = ' + strjoin(strsplit(strjoin(string(beta)), /extract), ', ')
      printf, 1, '# Bands used:', bands_temp
      printf, 1, '#     z| beta| Treal| wl_rest_1| wl_rest_2|   Tfit|     Mfit|    dMfit|'
      printf, 1
      close, 1
   endif
   
   for i = 0, nzmod-1 do begin
      ;; Setting z-related quantities
      z = zmod[i]
      if (z EQ 0) then dist_factor = (1+z)/(100.*3.0857d24)^2 else dist_factor = (1+z)/(lumdist(z, /silent)*3.0857d24)^2
      wl_rest = wl/(1+z)
      T_lowbound = T_CMB * (1+z)
      Tgrid = T_lowbound + (findgen(1000)+1)/5. & nTgrid = n_elements(Tgrid)
      
      ;; Cycle on beta
      for ibeta = 0, nbeta-1 do begin
         betatemp = beta[ibeta]
         ;; Creating SED array
         opacity = kappa_0 * (wl_rest/wl_0) ^(-betatemp)
         SED_unit_array = dblarr(nwl, nTgrid)
         phot_unit_array = fltarr(nTgrid, nbands)
         if not keyword_set(silent) then print, 'Computing grid for z = ' + strtrim(string(z, format = '(F0.2)'), 1) + ', beta = ' + $
                strtrim(string(betatemp, format = '(F0.2)'), 1)
         for j = 0, nTgrid-1 do begin
            SED_unit_array[*, j] = (ss_bbfunc(wl_rest, Tgrid[j], /intensity) - ss_bbfunc(wl_rest, T_CMB*(1+z), /intensity)) $
                                   * opacity * 1.9885d33 * dist_factor * conv_factor  ;; SEDs for 'unitary' mass (1 Msun)
            ;if not keyword_set(silent) then print, 'T # ' + strtrim(string(j)+1, 1) + ' of ' + strtrim(string(nTgrid, format = '(I0)'), 1) + $
            ;                                       ' (z = ' + strtrim(string(z, format = '(F0.2)'), 1) + $
            ;                                       ', beta = ' + strtrim(string(betatemp, format = '(F0.2)'), 1) + ') concluded'
         endfor
         str_photunitarray = (makesynthphot(wl, SED_unit_array, filters_in = bands_temp))
         phot_unit_array = str_photunitarray.fluxes
         
         ;; Fits (cycling on Treal)
         for j = 0, nTreal-1 do begin
            Ttemp = Treal[j]
            fname = fname_start + strtrim(string(Ttemp, format = '(F0.1)'), 1) + 'K_z' + $
                    strtrim(string(z, format = '(F0.2)'), 1) + '.dat'
            checkfile = file_test(folder_phot + fname)
            if checkfile then begin
               ;; Read photometry to fit
               readcol, folder_phot + fname_start + strtrim(string(Ttemp, format = '(F0.1)'), 1) + 'K_z' + $
                        strtrim(string(z, format = '(F0.2)'), 1) + '.dat', bands2fit, wl2fit, phot2fit, dphot2fit, $
                        format = '(A, D, D, D)', delimiter = '|', /silent
               ;; Band selection
               bdsel = where(bands2fit EQ bands_temp[0])
               for k = 1, nbands-1 do bdsel = [bdsel, where(bands2fit EQ bands_temp[k])]
               missingband = where(bdsel EQ -1, count)
               if count then begin
                  ;; Skip SED if one or more needed bands are missing
                  print, bands_temp[missingband] + ' band(s) missing for T_real = ' + strtrim(string(Ttemp, format = '(F0.1)'), 1) + $
                         ' K, z = ' + strtrim(string(z, format = '(F0.2)'), 1) + ', beta = ' + strtrim(string(betatemp, format = '(F0.2)'), 1) + $
                         '. Skipping to next SED.'
               endif else begin
                  ;; Determine Tfit
                  ratio = phot2fit[bdsel[0]] / phot2fit[bdsel[1]] ; Band ratio to reproduce
                  ratio_tcurve = phot_unit_array[0, *] / phot_unit_array[1, *]
                  Tfit = interpol(Tgrid, ratio_tcurve, ratio) & Tfit = Tfit[0]
                  if Tfit LE T_lowbound then begin
                     print, 'FIT ERROR: T_fit <= T_CMB'
                     print, 'The problem is probably caused by T_real being too close to T_CMB'
                     print, 'Results discarded for T_real = ', strtrim((string(Ttemp, format = '(F0.1)')),1), ', z = ', $
                            strtrim((string(z, format = '(F0.2)')),1)
                  endif else begin
                     ;; Reproduce synthetic photometry using Tfit
                     SED_unit = (ss_bbfunc(wl_rest, Tfit, /intensity) - ss_bbfunc(wl_rest, T_CMB*(1+z), /intensity)) $
                                * opacity * 1.9885d33 * dist_factor * conv_factor  ;; SED for T = Tfit, 'unitary' mass
                     phot_unit = (makesynthphot(wl, reform(SED_unit), filters_in = bands_temp)).fluxes
                     phot = phot_unit#Mgrid
                     ;; Chi^2 minimization
                     chi2 = fltarr(nMgrid)
                     for k = 0, nMgrid-1 do chi2[k] += total( ( (phot[*, k] - phot2fit[bdsel]) / dphot2fit[bdsel])^2 )
                     Mgrid_ind = where(chi2 LT min(chi2) + 9.)
                     Mfit_res = parfitter(Mgrid[Mgrid_ind], chi2[Mgrid_ind]) ;, x2plot = M2plot, y2plot = chi2plot)
                     Mfit = Mfit_res[0] & dMfit = Mfit_res[1]
                     if not keyword_set(silent) then begin
                        print, 'Results for band set #', strtrim(string(iset+1), 1), ', Treal = ', $
                               strtrim(string(Ttemp, format = '(F0.1)'), 1), ', z = ',  strtrim(string(z, format = '(F0.2)'), 1), ':'
                        print, 'Tfit = ', strtrim(string(Tfit, format = '(F0.2)'), 1)
                        print, 'Mfit = ' + strtrim(string(Mfit, format = '(E0.2)'), 1), ' +- ', strtrim(string(dMfit, format = '(E0.2)'), 1)
                     endif
                     if keyword_set(saveres) then begin
                        openw, 1, folder_fit + fname_out, /append
                        printf, 1, z, betatemp, Ttemp, wl2fit[bdsel]/(1+z), Tfit, Mfit, dMfit, format = $
                                '(F8.2, F6.2, F7.1, F11.2, F11.2, F8.2, E10.2, E10.2)'
                        close, 1
                     endif   
                  endelse
               endelse
            endif else begin
               print, 'No valid files found for T = ', strtrim(string(Ttemp),1), ', z = ', strtrim(string(z),1)
            endelse
         endfor
         
         if keyword_set(saveres) then begin
            openw, 1, folder_fit + fname_out, /append
            printf, 1
            close, 1
         endif
         
      endfor
      
      if keyword_set(saveres) then begin
         openw, 1, folder_fit + fname_out, /append
         printf, 1
         close, 1
      endif
      
   endfor
endfor


;stop


END
