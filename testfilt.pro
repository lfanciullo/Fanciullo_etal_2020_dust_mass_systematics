
;readcol, 'band_profiles/SPIRE_250mu.dat', wl1, fl1  ; Filter I use
;readcol, 'Herschel_SPIRE.PSW.dat', wl2, fl2         ; SVO filter (point source)
;readcol, 'Herschel_SPIRE.PSW_ext.dat', wl3, fl3     ; SVO filter (extended)

;readcol, '../Fanciullo_etal_backup/band_profiles/SPIRE_350mu.dat', wl1, fl1
;readcol, 'band_profiles/Herschel_SPIRE.PMW.dat', wl2, fl2

readcol, '../Fanciullo_etal_backup/band_profiles/SPIRE_500mu.dat', wl1, fl1
readcol, 'band_profiles/Herschel_SPIRE.PLW.dat', wl2, fl2

plot, wl1, fl1, col = 0, xtitle = textoidl('\lambda (Angstrom)'), ytitle = 'Transmission'
oplot, wl2, fl2, col = 250
;oplot, wl3, fl3, col = 75
oplot, wl2, fl2/max(fl2), col = 250, linestyle = 2 ; Maybe if I normalize?
;oplot, wl3, fl3/max(fl3), col = 75, linestyle = 2
