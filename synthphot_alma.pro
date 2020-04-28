
PRO synthphot_alma, wl_in, sed_in, filters_in, wl_alma, photom_alma, writefilt = writefilt, silent = silent

; Current band parameters are from the ALMA Observing Tool Cycle 6, Phase 2
; WRITEFILT option overrides everything else just to output the filter profiles used
;           Need to run ".r physconst.pro" to use it

if (size(sed_in))[0] EQ 1 then sed_in = reform(sed_in, (size(sed_in))[1], 1)
n_sed = (size(sed_in))[2]

nu = 1e-3*speedoflight()/wl_in  ;nu in GHz for c in m/s, wl_in in um
alma_bands_all = ['ALMA_3', 'ALMA_4', 'ALMA_5', 'ALMA_6', 'ALMA_7', 'ALMA_8', 'ALMA_9', 'ALMA_10']
nu_central_all = [97.5, 145., 203., 233., 343.5, 403., 679., 869.25]
bandwidth = 1.875
resol = .03125
nu_relative = resol*(indgen(641)-320) ;20 GHz width for resol = 31.25 MHz

if keyword_set(writefilt) then begin
   filters_alma = alma_bands_all
   nu_alma = nu_central_all
endif else begin
   match2, strtrim(filters_in, 2), alma_bands_all, a, b, na, nb
   filters_alma = filters_in[a]
   nu_alma = nu_central_all[b]
endelse
wl_alma = 1e-3 * speedoflight()/nu_alma ;nu in GHz, wl in um
nfiltout = n_elements(filters_alma)
photom_alma = dblarr(nfiltout, n_sed)


for i = 0, nfiltout-1 do begin
   nu_temp = nu_relative + nu_alma[i]
   filtshape = nu_temp*0.
   CASE filters_alma[i] OF
      'ALMA_6': IF_central = [-9., -7., 7., 9.] 
      'ALMA_9': IF_central = [-3., -1., 1., 3.] 
      'ALMA_10': IF_central = [-3., -1., 1., 3.] 
      ELSE: IF_central = [-7., -5., 5., 7.]
   ENDCASE
   nu_IF_temp = nu_alma[i] + IF_central
   nif = n_elements(IF_central)
   for k = 0, nif-1 do filtshape[where(abs(nu_temp - nu_IF_temp[k]) LT bandwidth/2)] = 1.
   if keyword_set(writefilt) then begin
      folder_temp = 'band_profiles/'
      fname_temp = alma_bands_all[i] + '.dat'
      wl_temp = 1e1*speedoflight()/nu_temp ;In Angstrom
      writecol, folder_temp+fname_temp, wl_temp, filtshape
   endif else begin
      for j = 0, n_sed-1 do begin
         goodsed = where(finite(sed_in[*, j]))
         if (min(nu_temp) LT min(nu[goodsed])) or (max(nu_temp) GT max(nu[goodsed])) then begin
            photom_alma[i, j] = !values.d_nan
            if not keyword_set(silent) then begin
               print
               print, 'SYNTHPHOT_ALMA: Band ', filters_alma[i], ', SED #',  strtrim(string(j), 1), $
                      ' is out of wavelength bounds. Flux for this band defaulting to NaN.' 
               print
            endif
         endif else begin
            sed_temp = interpol(sed_in[goodsed, j], nu[goodsed], nu_temp)
            bandwidth_total = nif * bandwidth
            photom_alma[i, j] = resol * total(sed_temp * filtshape) / bandwidth_total
         endelse
      endfor
   endelse
   ;;plot, /free, nu_temp, filtshape, yrange = [-.1, 1.1], col = 0
   ;;stop
endfor

;stop

end
