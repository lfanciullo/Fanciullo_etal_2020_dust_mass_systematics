FUNCTION make_opac, fold_data, material, T, wl_sed, plot = plot, saveplot = saveplot

;; NOTA BENE: The code assumes that all opacity files use the same
;;            wavelength array (which they do if they were created
;;            with opacity_reprocess.pro)

nwl_sed = n_elements(wl_sed)
nT = n_elements(T)
nmat = n_elements(material)
mac_complete = make_array(nwl_sed, nT, nmat, /float, VALUE=!values.f_nan)
if not keyword_set(plot) then plot = 0

;; Cycle on material
for j = 0, nmat-1 do begin
   ;; Find & read opacity files
   filepath_all = file_search(fold_data, material[j]+'*.xcat', expand_environment = 0)
   if (strmid(material[j], strlen(material[j])-1, 1) NE 'R') then $
      filepath_all = filepath_all(where(1-(strmatch(filepath_all, '*'+material[j]+'R*')))) ; Exclude E**R materials if looking for E**
   nT_fromfile = n_elements(filepath_all)
   T_fromfile = dblarr(nT_fromfile)
   for i = 0, nT_fromfile-1 do begin
      readcol, filepath_all[i], wl_temp, mac_temp, skipline = 4, /silent
      if (i EQ 0) then begin
         nwl_fromfile = n_elements(wl_temp)  ; NOTA: This assumes that all wavelength arrays have the same size (which should be the case if opacity_reprocess worked well) or that the first array is the shortest (!!!)
         wl_fromfile = wl_temp
         ;wl_fromfile = dblarr(nwl_fromfile, nT_fromfile) 
         mac_fromfile = dblarr(nwl_fromfile, nT_fromfile)
      endif
      ;wl_fromfile[*, i] = wl_temp
      mac_fromfile[*, i] = mac_temp
      junk = strsplit(filepath_all[i], '_/K.', /extract)
      T_fromfile[i] = junk[n_elements(junk)-3] ;(junk[11])
   endfor

   ;; Need to reorder the arrays
   realTorder = sort(T_fromfile)
   T_fromfile = T_fromfile(realTorder)
   mac_fromfile = mac_fromfile[*, realTorder]
   ;wl_fromfile = wl_fromfile[*, realTorder]
   
   ;; Interpolation on T
   mac_temp = make_array(nwl_fromfile, nT, /FLOAT, VALUE=!values.f_nan)
   mac_temp_smooth = make_array(nwl_fromfile, nT, /FLOAT, VALUE=!values.f_nan)
   lowT = where(T LT min(T_fromfile), count_lowT) ; Temperature array extremes
   hiT = where(T GT max(T_fromfile), count_hiT)
   for iwl = 0, nwl_fromfile-1 do begin
      Texists = where(finite(mac_fromfile[iwl, *]), count)
      mac_temp[iwl, *] = interpol(mac_fromfile[iwl, Texists], T_fromfile[Texists], T)
      ;; Avoid extrapolation beyond available temperature ranges:
      if count_lowT GT 0 then mac_temp[iwl, lowT] = mac_fromfile[iwl, 0]
      if count_hiT GT 0 then mac_temp[iwl, hiT] = mac_fromfile[iwl, nT_fromfile-1]
   endfor

   ;; Smoothing
   for iT = 0, nT-1 do mac_temp_smooth[*, iT] = smooth(mac_temp[*, iT], 10)  ; Boxcar smoothing, width = 10

   
   ;; PLOTS OF SMOOTHING SYSTEMATICS
   ;; Preparation
   if plot then begin
      loadct, 39
      !p.color = 0
      tvlct, 125, 125, 125, 1   ; Dark grey
      tvlct, 200, 200, 200, 2   ; Light grey
      tvlct, 225,   0, 225, 3   ; Pink
      tvlct, 150,  75,   0, 4   ; Brown
      ;allcol = [0, 250, 2, 75, 35, 3]
      ;allcol_BWfrly = [0, 0, 250, 250]
      ;allstyle_BWfrly = [0, 5, 0, 5]
      plotsym, 0, /fill
      Tstring = strtrim(string(fix(T)), 1) + 'K' & Tstring[0] = 'T = ' + Tstring[0]
      longwl = where(wl_fromfile GT 50.)
      ;vlongwl = where(wl_fromfile GT 500.)
      lstyles = [0, 5, 1, 2, 3]
      if keyword_set(saveplot) then begin
         ;caldat, systime(/julian), curr_month, curr_day, curr_year
         ;curr_date = strmid(curr_year, 10, 2) + '-' + string(curr_month, format='(I2.2)') + '-' + string(curr_day, format='(I2.2)')
         fold_pics = 'plots/'
         !x.charsize = 1.3
         !x.thick = 6
         !y.charsize = 1.4
         !y.thick = 6
         !p.charthick = 5
         !p.thick = 4
         !p.charsize = 1.05
      endif

      ;; Plot: percentage smoothing bias
      if keyword_set(saveplot) then begin
         set_plot, 'PS'
         device, filename = fold_pics + 'Fig2_smoothing_systematics_' +  material[j] + '.eps', /color, /encapsulated
      endif else begin
         window, /free
      endelse
      xlims = [40., 2000.]

      plot, wl_fromfile, 100 * (mac_temp_smooth[0, *]/mac_temp[*, 0]-1), col = 0, yrange = [-5., 5.], /nodata, $
            /xlog, xrange = xlims, /xstyle, /ystyle
      ;; Mark +/- 2% variation range
      polyfill, [xlims, reverse(xlims)], [2., 2., -2., -2.], color = 2
      ;; Overplot cage
      plot, wl_fromfile, 100 * (mac_temp_smooth[0, *]/mac_temp[*, 0]-1), col = 0, yrange = [-5., 5.], /nodata, /noerase, $
            xtitle = textoidl('\lambda (\mum)'), ytitle = 'Post-smoothing change (%)', /xlog, xrange = xlims, /xstyle, /ystyle
      ;; Main plots
      for i = 0, nT-1 do oplot, wl_fromfile[longwl], 100 * (mac_temp_smooth[longwl, i]/mac_temp[longwl, i] - 1.), $
                                               col = 0, linestyle = lstyles[i], thick = 1.5*!p.thick
      ;; Legend and labels
      legend, Tstring, psym = -0, linestyle = lstyles[0:(nT-1)], /right, /bottom, charsize = 1.4, thick = 1.5*!p.thick
      xyouts, sqrt(xlims[0] * xlims[1]), 3.75, material[j], alignment = .5, charsize = 2.5, charthick = 1.5*!p.charthick

      if keyword_set(saveplot) then begin
         device, /close
         set_plot, 'x'
      endif
   endif

   ; Interpolation on wl
   wlrange = where(wl_sed GE min(wl_fromfile) AND wl_sed LE max(wl_fromfile))
   for i = 0, nT-1 do mac_complete[wlrange, i, j] = interpol(mac_temp[*, i], wl_fromfile, wl_sed[wlrange])
endfor

return, mac_complete

END
