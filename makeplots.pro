
FUNCTION makembb, wl, params
  ;; This function creates a kappa_lambda vector from a wavelength and
  ;; a set of 3 (MBB) or 5 (broken power-law MBB) parameters 
  ;; * 3 parameters: kappa_0, lambda_0, beta
  ;; * 5 parameters: kappa_0, lambda_0, beta for lambda < lambda_break,
  ;;   lambda_break, excess at 500 um (see Gordon et al. 2014)

  parsel = where(finite(params)) & npar = n_elements(parsel)
  pargood = params[parsel]

  CASE npar OF
     3: begin
        opac = pargood[0] * (wl/pargood[1])^(-pargood[2])
     endcase
     5: begin
        opac = pargood[0] * (wl/pargood[1])^(-pargood[2])
        longwl = where(wl GE pargood[3])
        opac[longwl] *= (wl[longwl]/pargood[3])^(-alog(1 + pargood[4])/alog(pargood[3]/500.))
     endcase
     else: print, 'MAKEMBB error: # of parameters = ', npar, '; only 3 and 5 accepted.'
  ENDCASE

  return, opac
  
END


PRO plotfitres, whatplot, par2plot, comp, compfrac, Tdist, red = red, saveplot = saveplot

  fit_folder = 'fit_result_files/'
  pic_folder = 'plots/'
  col_all = [0, 30, 75, 250, 3]
  z2use = [0., 1., 3., 5., 7.]
  
  if comp EQ ['MBBtest'] then begin
     compstring = 'MBBtest'
     compstringadd = ''
  endif else begin
     compstring = strjoin(comp + '-' + strtrim(string(100 * compfrac, format = '(F4.1)'), 1), '+')
     if keyword_set(red) then compstringadd = '-red' else compstringadd = '-raw'
     compstring += compstringadd
  endelse
  if Tdist EQ '1T' OR Tdist EQ 'oneT' then begin
     Tdstring = '_oneT_'
     T2use = [20., 25., 30., 35., 40., 45., 50., 60., 80., 100.]
     xlims = [10., 110.]
     x_log = 0
     ymultfact = .85
     xtit = textoidl('T_{real} (K)')
  endif else begin
     ;;if Tdist = '2T' OR Tdist = 'twoT' then Tdstring = '_twoT_'
     Tdstring = '_twoT_'
     T2use = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1]
     xlims = [3e-5, 1.]
     x_log = 1
     ymultfact = .75
     xtit = textoidl('f_{W}')
  endelse
  
  fname = 'Fit_' + compstring + Tdstring + 'allbd-freeparams_flatpriors.dat'
  readcol, fit_folder + fname, Treal, z, Mfit, dMfit_hi, dMfit_lo, Tfit, dTfit_hi, $
           dTfit_lo, betafit, dbetafit_hi, dbetafit_lo, format = '(X,F,F,X,F,F,F,F,F,F,F,F,F,F)'

  x = Treal
  par2plot = strupcase(par2plot)
  CASE par2plot OF
     'B': begin
        y1 = betafit
        dy1_lo = dbetafit_lo
        dy1_hi = dbetafit_hi
        ylims = [.8, 2.4]
        xextra = xlims
        if (comp EQ ['MBBtest']) then yextra = [1.5, 1.5] else yextra = [0., 0.]
        ytit = textoidl('\beta_{fit}')
     end
     'M': begin
        y1 = Mfit / 1e8
        dy1_lo = dMfit_lo / 1e8
        dy1_hi = dMfit_hi / 1e8
        if (comp EQ ['MBBtest']) then begin
           ylims = [0., 2.4]
        endif else begin
           if keyword_set(red) then ylims = [0., 19.] else ylims = [0., 49.]
        endelse
        xextra = xlims
        yextra = [1., 1.]
        ytit = textoidl('M_{fit}/M_{real}')
     end
     'T': begin
        y1 = Tfit
        dy1_lo = dTfit_lo
        dy1_hi = dTfit_hi
        if (Tdist EQ '1T' OR Tdist EQ 'oneT') then begin
           ylims = [10., 190.]
           xextra = ylims
           yextra = ylims
        endif else begin
           ylims = [10., 110.]
           xextra = xlims
           yextra = [30., 30.]
           yextra2 = [100., 100.]
        endelse
        ytit = textoidl('T_{fit} (K)')
     end
     else: print, 'ERROR (MAKEPLOTS -> PLOTFITRES): ', par2plot, ' is not among the accepted values of PAR2PLOT'
  endcase

  if keyword_set(saveplot) then begin
     set_plot, 'PS'
     fname = 'Fig' + whatplot + '_' + par2plot + 'fit' + compstringadd + Tdstring + 'results.eps'
     device, filename = pic_folder + fname, /color, /encapsulated
  endif else begin
     window, /free
  endelse

  !p.multi = [0, 2, 3]
  multiplot, mxtitle = xtit, mytitle = ytit, mxTitSize = 1.3, myTitSize = 1.3
  if (Tdist EQ '2T' OR Tdist EQ 'twoT') then xtempname = textoidl(['10^{-4}', '0.001', '0.01', '0.1', '1'])
  
  for i = 0, n_elements(z2use)-1 do begin

     seltemp1 = where(z EQ z2use[i])
     match2, Treal[seltemp1], T2use, suba, subb
     seltemp2 = subb(where(subb GE 0))
     sel = seltemp1[seltemp2]
     
     if (i EQ n_elements(z2use)-1) then $
        plot, xextra, yextra, col = 0, yrange = ylims, xlog = x_log, xrange = xlims, /xstyle, /ystyle, xtickname = xtempname, /nodata else $
           plot, xextra, yextra, col = 0, yrange = ylims, xlog = x_log, xrange = xlims, /xstyle, /ystyle, /nodata
     oplot, xextra, yextra, col = 0, linestyle = 2
     if keyword_set(yextra2) then oplot, xlims, yextra2, col = 0, linestyle = 2
     ;oploterror, x[sel], y1[sel], dy1_hi[sel], ps = -8, col = 0, errcol = 0, /hibar
     ;oploterror, x[sel], y1[sel], dy1_lo[sel], ps = -8, col = 0, errcol = 0, /lobar
     oploterror, x[sel], y1[sel], dy1_hi[sel], ps = -8, col = col_all[i], errcol = col_all[i], /hibar ;, hatlength = !d.x_vsize * .02 
     oploterror, x[sel], y1[sel], dy1_lo[sel], ps = -8, col = col_all[i], errcol = col_all[i], /lobar ;, hatlength = !d.x_vsize * .02 

     ;; Printing z value with XYOUTS
     if x_log then xout = sqrt(max(xlims) * min(xlims)) else xout = (max(xlims) + min(xlims))/2
     xyouts, xout, max(ylims * ymultfact), 'z = ' + strtrim(string(z2use[i], format = '(F0.2)'), 1), charsize = 1.3, alignment = 0.5
     multiplot
     
  endfor

  ;; Final box
  plot, xextra, yextra, col = 0, xrange = xlims, yrange = ylims, xlog = x_log, /xstyle, /ystyle, xtickname = xtempname, /nodata
  oplot, xextra, yextra, col = 0, linestyle = 2
  if keyword_set(yextra2) then oplot, xlims, yextra2, col = 0, linestyle = 2
  for i = 0, n_elements(z2use)-1 do begin
     seltemp1 = where(z EQ z2use[i])
     match2, Treal[seltemp1], T2use, suba, subb
     seltemp2 = subb(where(subb GE 0))
     sel = seltemp1[seltemp2]
     oplot, x[sel], y1[sel], ps = -8, col = col_all[i]
  endfor
  
  multiplot, /reset

  if keyword_set(saveplot) then begin
     device, /close
     set_plot, 'x'
  endif
  
END



PRO plot2bands, whatplot, par2plot, bands, comp, optype = optype, z2plot = z2plot, beta_all = beta_all, saveplot = saveplot
  
  fit_folder = 'fit_result_files/'
  pic_folder = 'plots/'
  nbdcomb = (size(bands))[2]    ; # of band combinations
  nz2plot = n_elements(z2plot)
  nbeta = n_elements(beta_all)
  
  if comp EQ 'MBBtest' then begin
     opstring = ''
  endif else begin
     if not keyword_set(optype) then begin
        optype = 'red'
        print, 'KEYWORD ''OPTYPE'' NOT SET. DEFAULTING TO ''', strupcase(optype), ''''
     endif
     opstring = '-' + optype
  endelse
  
  ;; Plot quantities setup
  xtit = textoidl('T_{real} (K)')
  xlims = [10., 110.]
  if par2plot EQ 'T' then begin
     fname = 'Fig' + whatplot + '_2bdfit_Tfit_' + comp + opstring + '.eps'
     ytit = textoidl('T_{fit} (K)')
     ylims = [10., 110.]
     yextra = xlims
  endif else begin
     fname = 'Fig' + whatplot + '_2bdfit_Mfit_' + comp + opstring + '.eps'
     ytit = textoidl('M_{fit}/M_{real}')
     if (comp EQ ['MBBtest']) then begin
        ylims_all = [[0., 3.], [0., 9.], [0., 4.5]]
     endif else begin
        if (optype EQ 'red') then ylims_all = [[0., 15.], [0., 34.], [0., 24.]] else ylims = [[0., 55.], [0., 55.], [0., 55.]]
     endelse
     yextra = [1., 1.]
  endelse
  col_all = [250, 0, 75]
  
  ;; Multiplot preparation
  if keyword_set(saveplot) then begin
     set_plot, 'PS'
     device, filename = pic_folder + fname, /color, /encapsulated
  endif else begin
     window, /free
  endelse
  !p.multi = [0, nz2plot, nbdcomb]
  multiplot, mxtitle = xtit, mytitle = ytit, mxTitSize = 1.3, myTitSize = 1.3, xgap = .01

  ;; Cycle on band combinations
  for k = 0, nbdcomb-1 do begin
     bdtemp = bands[*, k]
     bdstring = bdtemp[0] + '+' + bdtemp[1]
     fname_temp = '2bdfit_' + comp + opstring + '_oneT_' + bdstring + '.dat'
     readcol, fit_folder + fname_temp, z, breal, Treal, Tfit, Mfit, dMfit, format = '(F,F,F,X,X,F,F,F)'
     if par2plot EQ 'M' then ylims = ylims_all[*, k]
     
     ;; Cycles on z
     for i = 0, nz2plot-1 do begin
        plot, xlims, yextra, xrange = xlims, /xstyle, yrange = ylims, /ystyle, col = 0, linestyle = 2

        ;; Cycles on beta
        for j = 0, nbeta-1 do begin
           selection = where(z EQ Z2plot[i] AND breal EQ beta_all[j])
           x = Treal[selection]
           if par2plot EQ 'T' then begin
              y = Tfit[selection]
              oplot, xlims, [2.725*(1+z2plot[i]), 2.725*(1+z2plot[i])], col = col_all[j], linestyle = 1 ; T_CMB
              oplot, x, y, psym = -8, col = col_all[j]
           endif else begin
              y = Mfit[selection] / 1e8
              dy = dMfit[selection] / 1e8
              oplot, x, y, psym = -8, col = col_all[j]
           endelse
        endfor

        ;; Legends and labels
        if k EQ 0 then begin
           xyouts, (xlims[0] + xlims[1]) / 2, ylims[1]*.85, 'z = ' + strtrim(string(z2plot[i], format = '(F0.2)'), 1), alignment = 0.5
           xyouts, (xlims[0] + xlims[1]) / 2, ylims[1]*.7, bdtemp[0] + ' + ' + bdtemp[1], alignment = 0.5
        endif else begin
           xyouts, (xlims[0] + xlims[1]) / 2, ylims[1]*.85, bdtemp[0] + ' + ' + bdtemp[1], alignment = 0.5
        endelse
        if ((i EQ nz2plot-1 AND k EQ 0) AND comp NE 'MBBtest')  then legend, $
           [textoidl( '\beta = ' + strtrim(string(beta_all[0], format = '(F0.2)'),1) ), [strtrim(string(beta_all[1:2], format = '(F0.2)'),1)]]$
           , charsize = .9, psym = 8, col = col_all, /bottom, /right, /clear
        
        multiplot
        
     endfor
     
  endfor
  
  multiplot, /reset
    
  if keyword_set(saveplot) then begin
     device, /close
     set_plot, 'x'
  endif
  
END




PRO makeplots, whatplot, Tdist = Tdist, saveplot = saveplot
  ;; Possible keys to add: comp, compfrac, Ts
  ;; Make opacity for arbitrary T #s rather than just 2 


;; Plot parameters setup
MAC_folder = 'MAC_files_reprocessed/'
SED_folder = 'synthetic_SEDs/'
fit_folder = 'fit_result_files/'
pic_folder = 'plots/'
loadct, 39
!p.color = 0
plotsym, 0, /fill
;; Defining additional colors
tvlct, 200, 200, 200, 1         ; Shades of grey
tvlct, 175, 175, 175, 2
tvlct, 150, 150, 150, 3
tvlct, 125, 125, 125, 4
tvlct, 100, 100, 100, 5
tvlct,  75,  75,  75, 6
tvlct, 225,   0, 225, 7         ; Fuchsia
tvlct, 255, 182, 193, 8         ; Light pink
tvlct, 150,  75,   0, 9         ; Brown
tvlct, 135, 206, 250, 10        ; Light sky blue
tvlct,   0, 100,   0, 11        ; Dark green
tvlct, 160,   42, 85, 12        ; Sienna
;col_all = [0, 30, 75, 250, 1]
linestyle_all = [3, 5, 1, 2, 0]
;pstyle_all = [0, 8, 4, 5, 3]
;psize_all = [1., 1., 1.3, 1.3, 1.2]
if keyword_set(saveplot) then begin
   !x.charsize = 1.2
   !x.thick = 6
   !y.charsize = 1.1
   !y.thick = 6
   !p.charthick = 5
   !p.thick = 4
   !p.charsize = 1.05
endif

if (size(whatplot))[1] NE 7 then whatplot = strtrim(string(whatplot), 1)



if whatplot EQ '1' then begin

   ;; Subset of opacities to plot
   matsample = ['E30R', 'X035', 'BE', 'AC', 'FAYA', 'FAY']
   nmat = n_elements(matsample)

   if keyword_set(saveplot) then begin
      set_plot, 'PS'
      device, filename = pic_folder + 'Fig' + whatplot + '_lab_opacity_selection.eps', /color, /encapsulated
   endif else begin
      window, /free
   endelse
   
   !p.multi=[0, 2, 3]
   xlims = [25., 2e3] & ylims = [.25, 5e3]
   multiplot, mxtitle = textoidl('\lambda (\mum)'), mytitle = textoidl('\kappa (cm^2 g^{-1})'), $
              mxTitSize = 1.3, myTitSize = 1.3
   
   ;; Cycle on materials
   for i = 0, nmat-1 do begin
      plot, xlims, ylims, col = 0, /xlog, /xstyle, /ylog, /ystyle, /nodata

      ;; Read all files for the material in question
      filepath_all_temp = file_search(MAC_folder, matsample[i]+'_*.xcat', expand_environment = 0)
      if (strmid(matsample[i], strlen(matsample[i])-1, 1) NE 'R') then filepath_all_temp = $
         filepath_all_temp(where(1-(strmatch(filepath_all_temp, '*'+matsample[i]+'R*'))))  ; Exclude E**R D17 files if looking for E**
      nT_fromfile = n_elements(filepath_all_temp) & T_fromfile = fltarr(nT_fromfile)

      ;; Create and sort list of available temperatures
      for j = 0, nT_fromfile-1 do begin
         junk = strsplit(filepath_all_temp[j], '_/K.', /extract)
         T_fromfile[j] = junk[(n_elements(junk) - 3)]
      endfor
      realTorder = sort(T_fromfile)  ; Sort T: alphabetical --> numerical order
      T_fromfile = T_fromfile(realTorder)
      filepath_all_temp = filepath_all_temp[realTorder]

      ;; Plot multi-T MAC
      for j = 0, nT_fromfile-1 do begin
         readcol, filepath_all_temp[j], wl_temp, mac_temp, skipline = 4, /silent ; Read wl + MAC from file
         oplot, wl_temp, mac_temp, col = 0, linestyle = linestyle_all[j]
      endfor

      ;; Labels and legends
      xyouts, sqrt(xlims[0]*xlims[1]), ylims[1]*.2, matsample[i], alignment = .5, charsize = 1.5, charthick = !p.charthick
      Tff_string  = strtrim(string(T_fromfile, format = '(I0)'), 1) + ' K'
      if (i EQ 0) then Tff_string[0] = 'T = ' + Tff_string[0]
      legend, Tff_string, linestyle = linestyle_all, charsize = .75, /bottom
      multiplot   
   endfor

   multiplot, /reset
   
   if keyword_set(saveplot) then begin
      device, /close
      set_plot, 'x'
   endif
   
endif



if whatplot EQ '2' then begin

   print
   print, 'Calling create_FIR_SED.pro to make Fig. 2. NOTA: You need to have compiled GRAMS_SYNTHPHOT.PRO and PHYSCONST.PRO'
   print

   ;; E30R silicates
   comp = ['E30R']
   compfrac = [1.]
   T_array = [10., 100., 300.]
   T_frac_array = [1., 0., 0.]  ; Content is not important as long as it has the same # of elements as T_array
   ;.r grams_synthphot
   ;.r physconst
   create_FIR_SED, save_sed = 0, plotsmoothing = 1, plot_SED = 0, save_plot = 1, $
                   comp = comp, fcomp = compfrac, T_all = T_array, fT_all = T_frac_array, z_all = z_array, /silent
   ;; BE carbon
   comp = ['BE']
   compfrac = [1.]
   T_array = [24., 100., 295.]
   T_frac_array = [1., 0., 0.]
   ;.r grams_synthphot
   ;.r physconst
   create_FIR_SED, save_sed = 0, plotsmoothing = 1, plot_SED = 0, save_plot = 1, $
                   comp = comp, fcomp = compfrac, T_all = T_array, fT_all = T_frac_array, z_all = z_array, /silent
   
endif



if whatplot EQ '3' then begin

   ;; Parameters
   col_all = [1, 2, 3, 4, 5, 6, 0]
   youtfactor = [.9, .8, .9, .9, 1.05, .9, .9]
   nwl = 201
   wls = 30.* 10.^(findgen(nwl)/100)

   ;; MBB opacity from the scientific literature: 7 models with 3 - 5 parameters each
   nmbb = 7
   opac_mbb_all = make_array(5, nmbb, /float, value = !values.f_nan)
   opac_mbb_all[0:2, 0] = [ 7.5,  230., 1.5]       ;Bertoldi et al. 2003
   opac_mbb_all[0:2, 1] = [30.0,  125., 2.0]       ;Robson et al. 2004
   opac_mbb_all[0:2, 2] = [  .4, 1200., 1.6]       ;Beelen et al. 2006
   opac_mbb_all[0:2, 3] = [34.7,  100., 2.2]       ;Weingartner & Draine 2001 (LMC)
   opac_mbb_all[0:2, 4] = [40.0,  100., 1.4]       ;Bianchi & Schneider 2007
   opac_mbb_all[*, 5] = [11.6, 160., 2.27, 294., .48] ;Gordon et al. 2014 (MW)
   opac_mbb_all[0:2, 6] = [  .7,  850., 2.0]          ;James et al. 2002
   opac_mbb_names = [$
                    'Bertoldi+03', $ ; 'Bertoldi et al. ''03', $ ;
                    'Robson+04', $   ; 'Robson et al. ''04', $ ;
                    'Beelen+06', $   ; 'Beelen et al. ''06', $ ;
                    'WD01', $        ; 'Weingartner & Draine ''01 (fit)', $ ;
                    'B&S07', $       ; 'Bianchi & Schneider ''07', $ ;
                    'Gordon+14', $   ; 'Gordon et al. ''14' $ ;
                    'James+02' $     ; 'James et al. ''02' $ ;
                    ]

   ;; Experimental opacity at different temperatures
   comp = ['E30R', 'BE']
   compfrac = [.7, .3]
   temp = [30., 100.]
   @makesed_params.idl
   op_lab_temp = make_opac(MAC_folder, comp, temp, wls, plot = plot)
   op_lab_cold = op_lab_temp[*, 0, 0] * compfrac[0] + op_lab_temp[*, 0, 1] * compfrac[1]  ; Raw
   op_lab_warm = op_lab_temp[*, 1, 0] * compfrac[0] + op_lab_temp[*, 1, 1] * compfrac[1]
   op_lab_cold_op = op_lab_temp[*, 0, 0] * compfrac[0] / redfact_sil + op_lab_temp[*, 0, 1] * compfrac[1] / redfact_c  ; Reduced
   op_lab_warm_op = op_lab_temp[*, 1, 0] * compfrac[0] / redfact_sil + op_lab_temp[*, 1, 1] * compfrac[1] / redfact_c
   
   if keyword_set(saveplot) then begin
      set_plot, 'PS'
      device, filename = pic_folder + 'Fig' + whatplot + '_opacity_lab_vs_lit_comparison.eps', /color, /encapsulated
   endif else begin
      window, /free
   endelse
   plot, [0., 0.], [0., 0.], col = 0, xrange = [50., 2000.], xstyle = 1, /xlog, yrange = [.1, 1e3], ystyle = 1, /ylog, /nodata   

   for i = 0, nmbb-1 do begin  ; Plotting modified blackbodies from the scientific literature
      ;; White margin for visibility
      oplot, wls, makembb(wls, opac_mbb_all[*, i]), thick = 1.5 * !p.thick, col = 255
      oplot, [opac_mbb_all[1, i]], [opac_mbb_all[0, i]], ps = 8, symsize = 1.5, col = 255
      ;; Grey or black actual line
      oplot, wls, makembb(wls, opac_mbb_all[*, i]), col = col_all[i]
      oplot, [opac_mbb_all[1, i]], [opac_mbb_all[0, i]], ps = 8, col = col_all[i]
   endfor
   for i = 0, nmbb-1 do xyouts, [opac_mbb_all[1, i]] * .95, [opac_mbb_all[0, i]] * youtfactor[i], opac_mbb_names[i], alignment = 1.  ; Carefully calibrated tag placement

   oplot, wls, op_lab_cold, col = 75  ; Lab dust opacity
   oplot, wls, op_lab_warm, col = 250
   
   cx = total(!x.window)/2.  ; X axis label
   cy = !y.window[0] - 0.08
   xyouts, cx, cy, textoidl('\lambda (\mum)'), /NORM, charsize = 1.5, align = 0.5
   cx = !x.window[0] - 0.08  ; Y axis label
   cy = total(!y.window)/2.
   xyouts, cx, cy, textoidl('\kappa (cm^2 g^{-1})'), /NORM, charsize = 1.5, align = 0.5, orientation = 90
   
   legend, ['T = 100 K', 'T = 30 K'], linestyle = [0, 0], col = [250, 75], /right
   
   if keyword_set(saveplot) then begin
      device, /close
      set_plot, 'x'
   endif
   
endif



if whatplot EQ '4' then begin

   ;; Experimental opacity at different temperatures
   nwl = 201
   wls = 30.* 10.^(findgen(nwl)/100)
   comp = ['E30R', 'BE']
   compfrac = [.7, .3]
   temp = [30., 100.]
   @makesed_params.idl
   op_lab_temp = make_opac(MAC_folder, comp, temp, wls, plot = plot)
   op_lab_cold = op_lab_temp[*, 0, 0] * compfrac[0] + op_lab_temp[*, 0, 1] * compfrac[1]  ; Raw
   op_lab_warm = op_lab_temp[*, 1, 0] * compfrac[0] + op_lab_temp[*, 1, 1] * compfrac[1]
   op_lab_cold_op = op_lab_temp[*, 0, 0] * compfrac[0] / redfact_sil + op_lab_temp[*, 0, 1] * compfrac[1] / redfact_c  ; Reduced
   op_lab_warm_op = op_lab_temp[*, 1, 0] * compfrac[0] / redfact_sil + op_lab_temp[*, 1, 1] * compfrac[1] / redfact_c
   
   if keyword_set(saveplot) then begin
      set_plot, 'PS'
      device, filename = pic_folder + 'Fig' + whatplot + '_opacity_correction.eps', /color, /encapsulated
   endif else begin
      window, /free
   endelse
   
   !p.multi=[0, 1, 2]   
   !p.position=[0.14, 0.3, 0.9, 0.9] ; Respectively: left, bottom, right, top axes

   ;; Opacity before/after reduction
   plot, [0., 0.], [0., 0.], col = 0, xrange = [50., 1000.], xstyle = 1, /xlog, xtickformat = '(A1)', $
         yrange = [.5, 1e3], ystyle = 1, /ylog, /nodata
   oplot, wls, op_lab_cold, col = 75
   oplot, wls, op_lab_warm, col = 250
   oplot, wls, op_lab_cold_op, col = 75, linestyle = 5
   oplot, wls, op_lab_warm_op, col = 250, linestyle = 5
   oplot, [850.], [.7], ps = 8, col = 0
   oplot, wls, makembb(wls, [.7,  850., 2.]), col = 0
   xyouts, [850.] * .95, [.7] * .9, 'James+02', alignment = 1.
   legend, ['T = 100 K', 'T = 30 K'], linestyle = [0, 0], col = [250, 75], position = [53., 9.]
   legend, ['Standard opacity', 'Reduced opacity'], linestyle = [0, 5], /bottom
   xyouts, .05, .6, textoidl('\kappa (cm^2 g^{-1})'), /NORM, charsize = 1.5, align = 0.5, orientation = 90 ; Y axis label

   ;; Ratio of before/after opacities
   !p.position = [.14, .1, .9, .3]
   plot, [0., 0.], [0., 0.], col = 0, xrange = [50., 1000.], xstyle = 1, /xlog, yrange = [.21, .49], /ystyle, /nodata
   oplot, wls, op_lab_cold_op/op_lab_cold, col = 75
   oplot, wls, op_lab_warm_op/op_lab_warm, col = 250
   xyouts, .52, .03, textoidl('\lambda (\mum)'), /NORM, charsize = 1.5, align = 0.5 ; X, Y axes labels
   xyouts, .05, .2, 'Ratio', /NORM, charsize = 1.5, align = 0.5, orientation = 90   ; Y axis label
   
   if keyword_set(saveplot) then begin
      device, /close
      set_plot, 'x'
   endif
   
endif



if whatplot EQ '5' then begin

   plotfitres, whatplot, 'M', ['E30R', 'BE'], [.7, .3], '1T', saveplot = saveplot
   
endif



if whatplot EQ '6' then begin

   plotfitres, whatplot, 'T', ['E30R', 'BE'], [.7, .3], '1T', saveplot = saveplot
   
endif



if whatplot EQ '7' then begin

   plotfitres, whatplot, 'B', ['E30R', 'BE'], [.7, .3], '1T', saveplot = saveplot   

endif



if whatplot EQ '8' then begin

   plotfitres, whatplot, 'M', ['E30R', 'BE'], [.7, .3], '1T', /red, saveplot = saveplot

endif



if whatplot EQ '9' then begin
   ;; T & beta (before/after reduction comparison)

   print
   print, 'Work in progress'
   print

   col_all = [0, 30, 75, 250, 3]
   comp = ['E30R', 'BE']
   compfrac = [.7, .3]
   compstring = strjoin(comp + '-' + strtrim(string(100 * compfrac, format = '(F4.1)'), 1), '+')
   fname_raw = 'Fit_' + compstring + '-raw_oneT_allbd-freeparams_flatpriors.dat'
   fname_red = 'Fit_' + compstring + '-red_oneT_allbd-freeparams_flatpriors.dat'
   readcol, fit_folder + fname_raw, Treal, z, Mfit_raw, dMfit_raw_hi, dMfit_raw_lo, Tfit_raw, dTfit_raw_hi, dTfit_raw_lo, betafit_raw, $
            dbetafit_raw_hi, dbetafit_raw_lo, format = '(X,F,F,X,F,F,F,F,F,F,F,F,F,F)'
   readcol, fit_folder + fname_red, Mfit_red, dMfit_red_hi, dMfit_red_lo, Tfit_red, dTfit_red_hi, dTfit_red_lo, betafit_red, dbetafit_red_hi, $
            dbetafit_red_lo, format = '(X,X,X,X,F,F,F,F,F,F,F,F,F,F)'
   z2plot = [0., 1., 3., 5., 7.]
   nz2plot = n_elements(z2plot)

   ;; Temperature comparison plot
   if keyword_set(saveplot) then begin
      set_plot, 'PS'
      device, filename = pic_folder + 'Fig' + whatplot + '_Tfit_opred_compare.eps', /color, /encapsulated
   endif else begin
      window, /free
   endelse
   ;; Setup
   xlims = [0., 200.]
   plot, xlims, xlims, xr = xlims, /xstyle, yr = xlims, linestyle = 2, col = 0, $
         xtit = textoidl('T_{fit} (K) for raw opacity'), ytit = textoidl('T_{fit} (K) for reduced opacity')
   ;; Cycle on z
   for i = 0, nz2plot-1 do begin
      z_temp = z2plot[i]
      indx = where(z EQ z_temp, count)
      ;npts = n_elements(indx)
      if count GT 0 then begin
         oploterror, Tfit_raw[indx], Tfit_red[indx], dTfit_raw_hi[indx], dTfit_red_hi[indx], ps = 8, /hibar, $
                     col = col_all[i], errcol = col_all[i]
         oploterror, Tfit_raw[indx], Tfit_red[indx], dTfit_raw_lo[indx], dTfit_red_lo[indx], ps = 8, /lobar, $
                     col = col_all[i], errcol = col_all[i]
      endif
   endfor
   if keyword_set(saveplot) then begin
      device, /close
      set_plot, 'x'
   endif

   ;; Beta comparison plot
   if keyword_set(saveplot) then begin
      set_plot, 'PS'
      device, filename = pic_folder + 'Fig' + whatplot + '_beta_opred_compare.eps', /color, /encapsulated
   endif else begin
      window, /free
   endelse
   ;; Setup
   xlims = [.8, 2.6]
   plot, xlims, xlims, xr = xlims, /xstyle, /ystyle, yr = xlims, linestyle = 2, col = 0, $
         xtit = textoidl('\beta for raw opacity'), ytit = textoidl('\beta for reduced opacity')
   ;; Cycle on z
   for i = 0, nz2plot-1 do begin
      z_temp = z2plot[i]
      indx = where(z EQ z_temp, count)
      ;npts = n_elements(indx)
      if count GT 0 then begin
         oploterror, betafit_raw[indx], betafit_red[indx], dbetafit_raw_hi[indx], dbetafit_red_hi[indx], ps = 8, /hibar, $
                     col = col_all[i], errcol = col_all[i]
         oploterror, betafit_raw[indx], betafit_red[indx], dbetafit_raw_lo[indx], dbetafit_red_lo[indx], ps = 8, /lobar, $
                     col = col_all[i], errcol = col_all[i]
      endif
   endfor
   if keyword_set(saveplot) then begin
      device, /close
      set_plot, 'x'
   endif
   
endif



if whatplot EQ '10' then begin
   ;; Effect of f_w on SED shape

   comp = ['E30R', 'BE']
   compfrac = [.7, .3]
   if not keyword_set(red) then red = 1
   compstring = strjoin(comp + '-' + strtrim(string(100 * compfrac, format = '(F4.1)'), 1), '+')
   if comp NE ['MBBtest'] then begin
      if red then compstring += '-red' else compstring += '-raw'
   endif
   fw_all = [.003, .03, .3]
   nfw = n_elements(fw_all)
   xlims = [30., 1000.]
   ylims = [.5, 5e3]

   if keyword_set(saveplot) then begin
      set_plot, 'PS'
      device, filename = pic_folder + 'Fig' + whatplot + '_twoT_SED_shapes.eps', /color, /encapsulated, xsize = 10., ysize = 15.
   endif else begin
      window, /free
   endelse

   !p.multi=[0, 1, nfw]
   multiplot, mxtitle = textoidl('\lambda (\mum)'), mytitle = textoidl('Flux (Jy)'), $
              mxTitSize = 1.2, myTitSize = 1.3
   
   for i = 0, nfw-1 do begin
      ;; Read SED files
      format2use = '(F0.' + strtrim(string(ceil(alog10(1/fw_all[i]))), 1) + ')'
      fwstring = strtrim(string(fw_all[i], format = format2use), 1)
      fname_spec = 'Spec_' + compstring + '_twoT-30.0K+100.0K-fw' + fwstring + '.dat'
      fname_phot = 'Phot_' + compstring + '_twoT-30.0K+100.0K-fw' + fwstring + '_z0.00.dat'
      readcol, SED_folder + 'Spectra/' + fname_spec, wl_spec, fl_spec
      readcol, SED_folder + 'Photometry/' + fname_phot, wl_phot, fl_phot, dfl_phot, format = '(X,F,F,F)'
      ;; Read fit file, plot fitted SED
      fname_fit = 'Fit_' + compstring + '_twoT_allbd-freeparams_flatpriors.dat'
      readcol, fit_folder + fname_fit, Tmod, zmod, Mfit, Tfit, betafit, format = '(X, F, F, X, F, X, X, F, X, X, F, X, X)'
      z2use = [0.]
      whatfit = where(Tmod EQ fw_all[i] AND zmod EQ z2use[0])
      Tfit2use = (Tfit[whatfit])[0]
      betafit2use = (betafit[whatfit])[0]
      Mfit2use = (Mfit[whatfit])[0]
      wl_fit = wl_spec * (1+z2use[0])
      makesed, wl_fit, fl_fit, Tfit2use, [1.], ['MBBtest'], [1.], M_d = Mfit2use, beta = betafit2use, z_in = z2use, /silent
      peak = where(fl_fit EQ max(fl_fit))
      ;; Plot
      plot, [0., 0.], [0., 0.], col = 0, xrange = xlims, xstyle = 1, /xlog, yrange = ylims, ystyle = 1, /ylog, /nodata
      oplot, [50., 50.], ylims, col = 3
      oplot, wl_fit, fl_fit, col = 250, thick = !p.thick * 1.5
      oplot, wl_fit[peak], fl_fit[peak], ps = 8, col = 250, symsize = 1.25
      oplot, wl_spec, fl_spec, col = 0
      oploterror, wl_phot, fl_phot, dfl_phot, ps = 4, col = 0, errcol = 0
      ;; Labels
      xyouts, (min(xlims)^.1 * max(xlims)^.9), (min(ylims)^.2 * max(ylims)^.8), textoidl('f_w = ') + fwstring, charsize = 1.3, alignment = 1.
      multiplot      
   endfor

   multiplot, /reset

   if keyword_set(saveplot) then begin
      device, /close
      set_plot, 'x'
   endif
   
endif



if whatplot EQ '11' then begin

   plotfitres, whatplot, 'M', ['E30R', 'BE'], [.7, .3], '2T', /red, saveplot = saveplot

endif



if whatplot EQ '12' then begin

   plotfitres, whatplot, 'T', ['E30R', 'BE'], [.7, .3], '2T', /red, saveplot = saveplot

endif




if whatplot EQ '13' then begin
   ;; Warm dust emission fraction

   comp = ['E30R', 'BE']
   compfrac = [.7, .3]
   if not keyword_set(red) then red = 1
   compstring = strjoin(comp + '-' + strtrim(string(100 * compfrac, format = '(F4.1)'), 1), '+')
   if comp NE ['MBBtest'] then begin
      if red then compstring += '-red' else compstring += '-raw'
   endif
   fw_all = [.3, .1, .03, .01]
   nfw = n_elements(fw_all)
   xlims = [30., 1000.]
   col_all = [0, 250, 75, 3]
   
   readcol, SED_folder + 'Spectra/Spec_' + compstring + '_oneT-30.0K.dat', wl, fl_30K
   readcol, SED_folder + 'Spectra/Spec_' + compstring + '_oneT-100.0K.dat', fl_100K, format = '(X,F)'
   
   if keyword_set(saveplot) then begin
      set_plot, 'PS'
      device, filename = pic_folder + 'Fig' + whatplot + '_warm_dust_contribution.eps', /color, /encapsulated
   endif else begin
      window, /free
   endelse

   plot, [0., 0.], [0., 0.], xrange = xlims, yrange = [0., 100.], /xstyle, /xlog, xtitle = textoidl('\lambda (\mum)'), $
         ytitle = 'Warm dust contribution (%)', xcharsize = 1.3, ycharsize = 1.5
   
   for i = 0, nfw-1 do begin
      format2use = '(F0.' + strtrim(string(ceil(alog10(1/fw_all[i]))), 1) + ')'
      fwstring = strtrim(string(fw_all[i], format = format2use), 1)
      fl_tot = fw_all[i] * fl_100K + (1 - fw_all[i]) * fl_30K
      oplot, wl, 100. * fw_all[i] * fl_100K/fl_tot, col = col_all[i], thick = !p.thick * 1.5
      xyouts, (min(xlims)^.05 * max(xlims)^.95), (min(100. * fw_all[i] * fl_100K/fl_tot) + 2.), textoidl('f_W = ') + fwstring, $
              col = col_all[i], charsize = 1.5, charthick = !p.charthick * 1.25, alignment = 1.
   endfor

   if keyword_set(saveplot) then begin
      device, /close
      set_plot, 'x'
   endif
   
endif



if whatplot EQ '14' then begin
   ;; 2-band fit (M)
   
   par2plot = 'M'
   comp = ['E30R', 'BE'] & compfrac = [.7, .3]
   compstring = strjoin(comp + '-' + strtrim(string(100 * compfrac, format = '(F4.1)'), 1), '+')
   z2plot = [6., 7.]
   beta_all = [1.5, 1.75, 2.]
   bands = $
      [['ALMA_7', 'ALMA_6'], $
       ['ALMA_6', 'ALMA_3'], $
       ['ALMA_7', 'ALMA_3']]
   
   plot2bands, whatplot, par2plot, bands, compstring, optype = 'red', z2plot = z2plot, beta_all = beta_all, saveplot = saveplot
   
endif



if whatplot EQ '15' then begin
   ;; 2-band fit (T)
   
   par2plot = 'T'
   comp = ['E30R', 'BE'] & compfrac = [.7, .3]
   compstring = strjoin(comp + '-' + strtrim(string(100 * compfrac, format = '(F4.1)'), 1), '+')
   z2plot = [6., 7.]
   beta_all = [1.5, 1.75, 2.]
   bands = $
      [['ALMA_7', 'ALMA_6'], $
       ['ALMA_6', 'ALMA_3'], $
       ['ALMA_7', 'ALMA_3']]
   
   plot2bands, whatplot, par2plot, bands, compstring, optype = 'red', z2plot = z2plot, beta_all = beta_all, saveplot = saveplot

endif




if whatplot EQ 'A1' then begin

   plotfitres, whatplot, 'M', ['MBBtest'], [1.], '1T', saveplot = saveplot

endif




if whatplot EQ 'A2' then begin

   plotfitres, whatplot, 'T', ['MBBtest'], [1.], '1T', saveplot = saveplot

endif




if whatplot EQ 'A3' then begin

   plotfitres, whatplot, 'B', ['MBBtest'], [1.], '1T', saveplot = saveplot

endif




if whatplot EQ 'A4' then begin

   plotfitres, whatplot, 'M', ['MBBtest'], [1.], '2T', saveplot = saveplot

endif




if whatplot EQ 'A5' then begin

   plotfitres, whatplot, 'T', ['MBBtest'], [1.], '2T', saveplot = saveplot

endif




if whatplot EQ 'A6' then begin

   plotfitres, whatplot, 'B', ['MBBtest'], [1.], '2T', saveplot = saveplot

endif




if whatplot EQ 'A7' then begin

   par2plot = 'M'
   comp = ['MBBtest'] ;& compfrac = [1.]
   compstring = comp[0] ;strjoin(comp + '-' + strtrim(string(100 * compfrac, format = '(F4.1)'), 1), '+')
   z2plot = [6., 7.]
   beta_all = [1.5, 1.75, 2.]
   bands = $
      [['ALMA_7', 'ALMA_6'], $
       ['ALMA_6', 'ALMA_3'], $
       ['ALMA_7', 'ALMA_3']]
   
   plot2bands, whatplot, par2plot, bands, compstring, z2plot = z2plot, beta_all = beta_all, saveplot = saveplot

endif



if whatplot EQ 'A8' then begin

   par2plot = 'T'
   comp = ['MBBtest'] ;& compfrac = [1.]
   compstring = comp[0] ;strjoin(comp + '-' + strtrim(string(100 * compfrac, format = '(F4.1)'), 1), '+')
   z2plot = [6., 7.]
   beta_all = [1.5, 1.75, 2.]
   bands = $
      [['ALMA_7', 'ALMA_6'], $
       ['ALMA_6', 'ALMA_3'], $
       ['ALMA_7', 'ALMA_3']]
   
   plot2bands, whatplot, par2plot, bands, compstring, z2plot = z2plot, beta_all = beta_all, saveplot = saveplot

endif




END
