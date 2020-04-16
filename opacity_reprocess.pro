
PRO opacity_reprocess, ref, savedata = savedata, logspace = logspace, plot = plot, saveplot = saveplot;, date = date

;; MODIFICATIONS:
;; * D17B is spread across 2 folders; code that
;; * skipline seems to be superfluous now
  
;; Reprocesses the original opacity (MAC) files for use by other scripts
;;  (add individual things it does)
;; 
;; KEYWORDS:
;; REF: What data to reprocess (indicated by article reference)
;; SAVEDATA: If set, saves the processed opacity (default: don't save)
;; LOGSPACE: ...
;; PLOT: ...
;; SAVEPLOT: ...
;; (DATE: ...)
  
;; Example call: opacity_reprocess, 'D17B', /savedata, /plot


if not keyword_set(savedata) then savedata = 0
if not keyword_set(plot) then plot = 0
if not keyword_set(saveplot) then saveplot = 0

fold_in = 'MAC_files_original/'
fold_out = 'MAC_files_reprocessed/'
wl_meta = 10.^(findgen(1000)/(1000/4)) ; "Overarching" wl array: 1000 pt from 1 um to 10 mm

ref = strupcase(ref)
CASE ref OF
   'M98':begin
      materials = ['AC', 'BE', 'FOR', 'FAY', 'FAYA']                   ; Available materials
      T_all = [24., 100., 160., 200., 295.]                            ; Available temperatures
      wl_intpol = wl_meta(where(wl_meta GT 25. AND wl_meta LT 2000.))  ; Available wavelength range
      ref_folder_all = ['Mennella_98/']                                ; Folders over which the data is saved
      x_type = 'wavelength'                                            ; Opacity is as a function of wavelength
      extension = '.FIR'                                               ; File extension
      material_fullnames = materials
      Tstring = strtrim(string(fix(T_all)), 1) ;+ 'K'
   end
   'D17A':begin
      materials = ['x035', 'x040', 'x050A', 'x050B']
      T_all = [10., 30., 100., 200., 300.]
      wl_intpol = wl_meta(where(wl_meta GT 25. AND wl_meta LT 1000.))
      ref_folder_all = ['EXPERIMENT_KD_20170711/']
      x_type = 'wavenumber'
      extension = '.data.csv'
      material_fullnames = 'MAC_MgO(1-x)-xSiO2_' + materials + '_5_1000mu' ; + Tstring
      Tstring = '_' + strtrim(string(fix(T_all)), 1) + 'K'
   end
   'D17B':begin
      materials = [ 'E10', 'E10R', 'E20', 'E20R', 'E30', 'E30R', 'E40', 'E40R']
      T_all = [10., 30., 100., 200., 300.]
      wl_intpol = wl_meta(where(wl_meta GT 5. AND wl_meta LT 1000.))
      ref_folder_all = ['EXPERIMENT_KD_20170822/', 'EXPERIMENT_KD_20170823/']
      x_type = 'wavenumber'
      extension = '.txt.data.csv'
      material_fullnames = 'MAC_Mg(1-x)FexSiO3_' + materials + '_5_1000mu'
      material_fullnames[0] = 'MAC_Mg(1-x)FexSiO3_E10_5_800mu'             ; E10 data only covers up to 800 um instead of 1000
      Tstring = '_' + strtrim(string(fix(T_all)), 1) + 'K'
   end
   else: begin
      print
      print, 'WARNING (OPACITY_REPROCESS): Unrecognized reference. Please use one of the following:'
      print, 'M98,   D17A,   D17B'
      print
   end
ENDCASE
nwl_intpol = n_elements(wl_intpol)
nreffold = n_elements(ref_folder_all)
nmat = n_elements(materials)
nT = n_elements(T_all)
if ref EQ 'M98' then Tstring[nT-1] = ''


for i = 0, nreffold-1 do begin  ; Cycle on folders
   ref_folder = ref_folder_all[i]
   ;; Search only for materials in local folder (not important for M98
   ;; and D17A, but D17B the E** and E**R materials are in different folders)
   if ref NE 'D17B' then begin
      materials_local = materials
      material_fullnames_local = material_fullnames
   endif else begin
      if strmatch(ref_folder, '*20170822*') then begin
         materials_local = materials[where(strmatch(materials, '*0'))]
         material_fullnames_local = material_fullnames[where(strmatch(material_fullnames, '*0_*'))]
      endif
      if strmatch(ref_folder, '*20170823*') then begin
         materials_local = materials[where(strmatch(materials, '*0R'))]
         material_fullnames_local = material_fullnames[where(strmatch(material_fullnames, '*0R*'))]
      endif
   endelse
   nmatloc = n_elements(materials_local)
   
   for j = 0, nmatloc-1 do begin  ; Cycle on material (in local folder)
      opac_1mat_allT = fltarr(nT, nwl_intpol)
      opac_temp_allT = make_array(nT, nwl_intpol, /FLOAT, VALUE=!values.f_nan)

      for k = 0, nT-1 do begin  ; Cycle on temperature
         if ref EQ 'D17A' and T_all[k] EQ 30. then filename = material_fullnames_local[j] + '_10K' else $
            filename = material_fullnames_local[j] + Tstring[k]
         ref_folder_temp = ref_folder
         if ref EQ 'D17A' then ref_folder_temp += filename + '/'
         if ref EQ 'D17B' then ref_folder_temp += filename + '.txt/'
         readcol, fold_in + ref_folder_temp + filename + extension, wl, opac, /silent
         ;; Reprocess data if needed:
         ;;  wavenumber (cm^-1) --> wl (um)
         if x_type EQ 'wavenumber' then begin 
            wl = 1e4/reverse(wl) & opac = reverse(opac) 
         endif
         ;;  M98 files (except three) are in logspace
         if ref EQ 'M98' and Tstring[K] NE '' then begin
            wl = 10.^wl
            opac = 10.^opac
         endif
         ;;  Interpolation on wavelength (taking care not to extrapolate)
         outsiders = where(wl_intpol LT min(wl) OR wl_intpol GT max(wl), count)
         opac_1mat_allT[k, *] = 10.^interpol(alog10(opac), alog10(wl), alog10(wl_intpol))
         if count GT 0 then opac_1mat_allT[k, outsiders] = !values.f_nan
      endfor
      
      ;; Interpolation on T
      keepthis = [-1]
      for k = 0, nwl_intpol-1 do begin
         opac_temp = opac_1mat_allT[*, k]
         gooddata = where(finite(opac_temp), count)
         ;; To do the interpolation, we need that at least some data
         ;; (including the two interpolation extremes) be finite:
         dothis = count GT 0 AND (finite(opac_temp[0]) AND finite(opac_temp[nT-1]))
         if dothis then begin
            opac_temp_allT[*,k] = interpol(opac_temp[gooddata], T_all[gooddata], T_all)
            if (keepthis EQ [-1]) then keepthis = [k] else keepthis = [keepthis, [k]]
         endif
      endfor
      if keepthis NE [-1] then begin
         opac_1mat_allT_intpol = opac_temp_allT[*,keepthis]
                                ;opac_1mat_allT_intpol = opac_1mat_allT[*,keepthis]
      endif else begin
         print, 'ERROR (OPACITY_REPROCESS): Interpolation on temperature of ', materials_local[j], ' failed'
         opac_1mat_allT_intpol = opac_1mat_allT
      endelse

      ;; Saving
      if keyword_set(savedata) then begin
         if ref EQ 'M98' then Tstring_out = '_' + strtrim(string(fix(T_all)), 1) + 'K' else Tstring_out = Tstring
         for k = 0, nT-1 do begin
            if keepthis NE [-1] then nwl_final = n_elements(keepthis) else nwl_final = nwl_intpol
            output_structure = {wl:0., opac:0.} & output_structure = replicate(output_structure, nwl_final)  ; Output structure
            if keepthis NE [-1] then output_structure.wl = wl_intpol[keepthis] else output_structure.wl = wl_intpol
            opac_temp = opac_1mat_allT_intpol[k, *]
            output_structure.opac = reform(opac_temp)
            fname_interpol = strupcase(materials_local[j]) + Tstring_out[k] + '_interpol.xcat'  ; Output file name
            write_xcat, output_structure, fold_out + fname_interpol, /silent
         endfor
      endif
      
      ;;; PLOTS ;;;
      if plot then begin  ;; (This can be removed once I have the initializer)
         loadct, 39
         !p.color = 0
         tvlct, 125, 125, 125, 1 ; Dark grey
         tvlct, 200, 200, 200, 2 ; Light grey
         tvlct, 225,   0, 225, 3 ; Pink
         tvlct, 150,  75,   0, 4 ; Brown
         linestyle_all = [3, 5, 1, 2, 0]
         allcol = [0, 250, 75, 2, 35, 3]
         Tstringlgd = Tstring & Tstringlgd[0] = 'T = ' + Tstringlgd[0] 
         plotsym, 0, /fill
         if keyword_set(saveplot) then begin
            caldat, systime(/julian), curr_month, curr_day, curr_year
            curr_date = strmid(curr_year, 10, 2) + '-' + string(curr_month, format='(I2.2)') + '-' + string(curr_day, format='(I2.2)')
            fold_pics = '~/Documents/Postdoc/Science/FIR_data/Plots_FIR_data/'
            !x.charsize = 1.2
            !x.thick = 6
            !y.charsize = 1.05
            !y.thick = 6
            !p.charthick = 3
            !p.thick = 5
            !p.charsize = 1.1   ;1.4;
         endif
         
         ;; Opacity plot for sanity check
         if keyword_set(saveplot) then begin
            set_plot, 'PS'
            device, filename = fold_pics +  materials_local[j] + '_Opacity_Interpolated.eps', /color, /encapsulated
         endif else begin
            window, /free
         endelse
         plot, wl_intpol, opac_1mat_allT[0,*], col = 0, /xlog, /ylog, /nodata, title = materials_local[j], $
               xtitle = textoidl('\lambda (\mum)'), ytitle = textoidl('\kappa (cm^2g^{-1})')
         for k = 0, nT-1 do oplot, wl_intpol, opac_1mat_allT[k, *], col = 0, linestyle = linestyle_all[k]
         legend, strmid(Tstring, 1), psym = -0, linestyle = indgen(nT), /bottom
         if keyword_set(saveplot) then begin
            device, /close
            set_plot, 'x'
         endif
         print, 'Plot made for material = ' + strtrim(materials_local[j], 1) + ' (#' + strtrim(string(j+1), 1) + ' of ' + strtrim(nmat, 1) $
                + '). Type ''.c'' to proceed to next material.'
         stop
      endif

   endfor
   
endfor      

;stop


END

