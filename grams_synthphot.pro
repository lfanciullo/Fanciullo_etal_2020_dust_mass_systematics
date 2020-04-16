@~/idl-libs/match2
FUNCTION makestr,x
 return,strtrim(string(x),2)
END

PRO loadfits,fitspath=fitspath,fitsfile=fitsfile
;Updated 2013-06-18 added FITSPATH keyword
;Written July 14 2011
;Find the location of the GRAMS_filters.fits file
 common GRAMS_filters,s,path
 
 filename='GRAMS_filters.fits' & if keyword_set(fitsfile) then filename=makestr(fitsfile)+'.fits' 
 cd,current=curr
 if n_elements(path) eq 0 then path=curr & print,'Searching for '+makestr(filename)+' in current directory...'
 if 1-file_test(path+filename) then begin
  print,"GRAMS_filters access error: either the path hasn't been provided, or the file does not exist in the folder specified"
  junk=''
  read,junk,prompt="Enter path to file or folder in which to search (e.g., '~/Midlife_Chrysalis/DrKillinger/'): "
  if file_test(junk+filename) then path=junk $
  else path=file_dirname(file_search(junk,filename))+'/'
 endif
 print,'Reading in '+filename
 s=mrdfits(path+filename,1)
END

FUNCTION INTTRAP,x,y
;Updated July 16 2011
;	If there are only two wavelength points in X, output is 0.5*DeltaX*DeltaY (simplest trapezoid)
;Updated July 14 2011 to handle 2-D Y arrays
;	The format for Y is FIXED: NS by NX, where NX = n_elements(X)
;	The output will be a 1 by NS array of integrals over NX
;Integrate using trapezoid rule
 str='' & str1='' & str2='' & str3='' & sy=size(y)
 if sy[0] gt 1 then junk=execute("str+='*,' & str1+='transpose(cmreplicate(' & str2+=',sy[1]))' & str3+=',2'")
 n=n_elements(x) & k=bsort(x) & x=x[k] & junk=execute('y=y['+str+'k]')
 xu=x[1:n-1] & xl=x[0:n-2] & junk=execute('yu=y['+str+'1:n-1] & yl=y['+str+'0:n-2]')
 junk=execute('integrand='+str1+'(xu-xl)'+str2+'*(yu+yl)')
 if (size(integrand))[0] eq 1 then integral=0.5*integrand else integral=0.5*total(integrand,2)
 return,integral
END

FUNCTION FILTERCOVERAGE,LIN,LFIL,RFIL,KCOM
;Updated Dec 12 2012 - Restrict the filter wavelength range to the region with response > 5% of maximum
;Updated May 14 2012 - Return the location of the common wavelength points instead of the wavelength values
;Updated April 15 2011
;Written May 25 2009
;- Measure the overlap between the wavelength grids of the spectrum and the filter, output 0 < COV <=1.
;- NOTE: If the coverage COV is between 0.95 and 1, it's probably due to the difference in resolution between LIN and LFIL,
;	and there's probably 100% coverage, so COV is forced to 1.
;- INPUT: LIN - wavelength grid for the spectrum to be convolved
;	  LFIL, RFIL - filter wavelength grid and RSR
;- OUTPUT:COV - coverage between 0 and 1.
;	  KCOM - location of wavelengths in LIN that will be used in the common wavelength grid.
;		 This option avoids the interpolation of FIN which may be very noisy at small scales.
 ;lfilt=lfil[where(rfil gt 0.)] & minfil=min(lfilt,max=maxfil)
 lfilt=lfil[where(rfil gt .05*max(rfil))] & minfil=min(lfilt,max=maxfil)
 mask=lin gt minfil and lin lt maxfil
 if max(mask) ne 0 then begin
  ;;kcom1=where(mask) & kcom2=uniq(lin[kcom1],bsort(lin[kcom1])) & kcom=kcom1[kcom2]
  ;kcom1=where(mask) & kcom2=uniq2(lin[kcom1],bsort(lin[kcom1])) & kcom=kcom1[kcom2]
  kcom1=where(mask) & kcom2=uniq2(lin[kcom1]) & kcom=kcom1[kcom2]
  cov=(max(lin[kcom])-min(lin[kcom]))/(maxfil-minfil) ; what fraction of lfilt is covered by lin?
  if cov ge .95 then cov=1.
 endif else begin
  kcom=-1 & cov=0.
 endelse
 return,min([1.,cov])
END

FUNCTION bandfold,lin,fin,lfil,rfil,filter_name,flam=flam,verbose=verbose
;Updated July 14 2011
;	FIN can now either be 1 by n_elements(LIN) (single source) or NS x n_elements(LIN) (NS sources)
;Updated April 15 2011
;Written June 19 2009
;DESCRIPTION: Input spectrum (LIN,FIN) is folded through filter band (LFIL,RFIL)
;INPUTS: LIN,FIN - INPUT SPECTRUM, wavelength in microns, flux in F_nu or (if keyword FLAM is set) F_lam units
;        LFIL,RFIL - FILTER wavelengths and response in e-/energy. If response is in e-/photon, pass lambda*response to this module.
;OPTIONAL KEYWORDS:
;        FLAM - If set, the integral performed is INT(dlfil*rfil*fin)/INT(dlfil*rfil). Otherwise, extra lfil term in both integrands.
;        VERBOSE - Print warnings for partial or no overlap
;OUTPUTS:FOUT - output flux, same units as input flux
;NOTE:   COV<0.5 forces output flux to zero.
 cov=filtercoverage(lin,lfil,rfil,kcom) & if kcom[0] ne -1 then lcom=lin[kcom] else lcom=[0.]
 sy=size(fin) & for i=1,3 do junk=execute("str"+strtrim(string(i),2)+"=''")
 if sy[0] gt 1 then junk=execute("str1='transpose(cmreplicate(' & str2=',sy[1]))' & str3='*,'")
 junk=execute('fout='+str1+'cov gt 0.5'+str2) ;force output flux to zero if less than 50% overlap with filter
 if cov ne 0. and max(lcom) ne 0. then begin
  rfil1=interpol(rfil,lfil,lcom)
  junk=execute('fin1=fin['+str3+'kcom]')
  junk=execute('weight='+str1+'rfil1*lcom^(2.*keyword_set(flam)-1.)'+str2)
  fout=inttrap(lcom,weight*fin1)/inttrap(lcom,weight)
  if keyword_set(verbose) then begin
   if cov le 0. then print,'Spectrum outside '+filter_name+' filter' else if cov lt .5 then print,$
    'Not enough coverage to estimate filter-folded flux in '+filter_name+' filter' else if cov lt .8 then print,$
    'Warning: spectrum not entirely inside '+filter_name+' filter, output flux will be a lower limit.'
  endif
  if keyword_set(verbose) then print,'coverage for '+filter_name+' = ',string(cov,2)
 endif
 return,fout
END

PRO filterfold,lin,fin,fout,lout,flam=flam,verbose=verbose
;Updated July 14 2011
;       FIN can now either be 1 by n_elements(LIN) (single source) or NS x n_elements(LIN) (NS sources)
;Updated April 15 2011 to include MIPS 70,160, AKARI and WISE passbands
;       - Now reads in zeropoint and reference wavelength information from FITS file
;       - The IRAC and WISE bands are in e-/photon, so pass lambda*response to BANDFOLD
;Written June 19 2009
;This means the weight for IRAC Fnu calculation is R/lambda and for everything else, it is R/lambda^2.
;For Flambda calculation, the corresponding weights are R*lambda and R.
;DESCRIPTION:
;        INPUT SPECTRUM (LAMIN,FIN) is folded through all available filters
;INPUTS: LIN,FIN - Input spectrum. Wavelength in microns, flux in F_nu or (if keyword FLAM is set) F_lam units.
;OUTPUT: FOUT - Same number of elements as there are filters in the reference FITS file.
;        LOUT (optional) - Store the reference wavelengths from the FITS file in this variable.
;OPTIONAL KEYWORDS: FLAM - If set, FIN is in F_lam units instead of F_nu.
 common GRAMS_filters,s,path
 if n_elements(s) eq 0 then loadfits
 sy=size(fin) & str1='' & str2='' & if sy[0] gt 1 then junk=execute("str1+=',sy[1]' & str2+='*,'")
 lout=s.lamref & junk=execute('fout=replicate(0.'+str1+',n_elements(lout))')
 str='' & if keyword_set(verbose) then str+=',/verbose' & if keyword_set(flam) then str+=',/flam'
 for i=0,n_elements(lout)-1 do begin
  ;If RSR is in e-/photon, convert to e-/energy
  if strmatch(s[i].notes[0],'*energy*') then rsr=s[i].rsrfil else rsr=s[i].lamfil*s[i].rsrfil
  ;THIS BIT IS NEW: DOUBLE THE RESOLUTION OF THE INPUT SPECTRA BEFORE PASSING TO BANDFOLD
  lshift=(0.5*(lin+shift(lin,1)))[1:sy[2]-1]
  lin1=[lin,lshift] & lin1=lin1[bsort(lin1)] & lin1=lin1[uniq2(lin1)]
  fin1=dblarr(sy[1],n_elements(lin1))
  for j=0l,sy[1]-1 do begin
   fin1[j,*]=interpol(fin[j,*],lin,lin1)
  endfor
  junk=execute('f=bandfold(lin1,fin1,s[i].lamfil,rsr,s[i].filter_name'+str+') & fout['+str2+'i]=f')
  ;
  ;junk=execute('f=bandfold(lin,fin,s[i].lamfil,rsr,s[i].filter_name'+str+') & fout['+str2+'i]=f')
 endfor
END

FUNCTION filterfolderr,lspec,fspec,dfspec,options=str,median=fphotmed
;
;
 common GRAMS_filters,s,path

 nspec=n_elements(lspec) & nfil=n_elements(s)
 if n_elements(str) eq 0 then str=''
 Nrruns=1000
 fphot=fltarr([nfil,Nrruns])
 for i=0L,Nrruns-1 do begin
  fs=fspec+dfspec*randomn(seed,(size(fspec))[1:2])
  junk=execute('filterfold,lspec,fs,fp,lout'+str)
  fphot[0:nfil-1,i]=fp
 endfor
 fphotmed=median(fphot,dim=2,/ev,/do)
 return,1.4826*madm(fphot,dim=2)
END

PRO synthphot,l,f,fphot,mphot,flam=flam,dfspec=dfspec,dfphot=dfphot,dmphot=dmphot,$
        filter_name=filter_name,fitspath=fitspath,fitsfile=fitsfile,$
        verbose=verbose,help=help
;+
;Given an SED, f(_nu) vs l(micron), return the synthetic photometry in filters specified in GRAMS_filters.fits.
;PRO synthphot,l,f,fphot,mphot,flamm=flam,dfspec=dfspec,dfphot=dfphot,dmphot=dmphot,$
;   filter_name=filter_name,fitspath=fitspath,fitsfile=fitsfile,verbose=verbose,help=help,version=version
;   L is a NWx1 array of wavelengths in micron, F is a NSxNW array of fluxes where NS=#sources and NW=#wavelengths.
;       ALL NS FLUX ARRAYS HAVE TO BE ON THE SAME WAVELENGTH GRID!
;   Keywords:
;   FLAM: Set if input is in F_lambda units instead of F_nu.
;   FITSPATH: Path to GRAMS_filters.fits (or equivalent) file. If not set, the code searches through local folder.
;       If not found, user is prompted for a path.
;   FITSFILE: Name of filters file, if different from "GRAMS_filters".
;   DFSPEC/DFPHOT/DFMAG: If the flux uncertainties are provided in DFSPEC, the uncertainty in the synthetic photometry
;       is output in DFPHOT/DFMAG.
;   VERBOSE: Fires helpful bullets of information.
;Updated 2016-05-17: Keywords SAGE, HELP, and VERSION deprecated.
;Updated November 9 2014: Added DFSPEC, DFPHOT, DMPHOT keywords. The code now outputs an uncertainty estimate. 
;Updated March 7 2013: Added FITSFILE keyword, allows the use of an alternate GRAMS_filters file that may have non-standard filters
;Updated May 14 2012: FILTERFOLD and FILTERCOVERAGE routines modified to return input spectrum sampled at common wavelength points
;Updated July 25 2011: added FILTER_NAME keyword
;	FILTER_NAME can be set to an array of values included in GRAMS_filters.fits to output photometry only in those filters.
;Updated July 14 2011 to handle 2-D Y arrays
;       The format for Y is FIXED: NS by NX, where NX = n_elements(X)
;       The output will be a 1 by NS array of integrals over NX
;Updated July 13 2011 - added FITSPATH keyword
;Updated June 17 2011 - added HELP and SAGE keywords
;Written August 10 2009
 common GRAMS_filters,s,path

 str='' & if keyword_set(flam) then str+=',/flam' & if keyword_set(verbose) then str+=',/verbose'

 if n_elements(fitspath) ne 0 then path=fitspath
 str1='' & if keyword_set(fitspath) then str1+=',fitspath="'+makestr(fitspath)+'"' & if keyword_set(fitsfile) then $
	str1+=',fitsfile="'+makestr(fitsfile)+'"'
 if n_elements(s) eq 0 then junk=execute('loadfits'+str1)
 if n_elements(dfspec) eq 0 then dfspec=0.*f ;Checking for uncertainties, and setting them to zero if absent.

 junk=execute('filterfold,l,f,fphot,lout'+str)
 junk=execute('filterfold,l,f+dfspec,fphothi,lout'+str)
 junk=execute('filterfold,l,f-dfspec,fphotlo,lout'+str)
 dfphot=0.5*(fphothi-fphotlo)
 ;dfphot=filterfolderr(l,f,dfspec,options=str,median=fphot1) ;NEW 2016-05-05
 str1='' & str2='' & if (size(fphot))[0] gt 1 then junk=execute("str1+='transpose(cmreplicate(' & str2+=',string((size(fphot))[1],2)))'")
 if keyword_set(flam) then mphot=0.*fphot else junk=execute('mphot=-2.5*alog10(fphot/'+str1+'s.fnuref'+str2+')')
 if keyword_set(sage) then begin
  fphot=fphot[*,0:11] & mphot=mphot[*,0:11] 
 endif
 if keyword_set(filter_name) then begin
  match,strtrim(s.filter_name,2),strtrim(filter_name,2),a,b ;s[a].filter_name = filter_name[b]
  if a[0] eq -1 then begin
   print,'ERROR! None of the elements in FILTER_NAME are valid filters! Please see the GRAMS_filters.fits structure!'
  endif else begin
   if min(b[bsort(b)] ne indgen(n_elements(filter_name))) eq 0 then print,'Only some filters in FILTER_NAME were found!'
   fphot=fphot[*,a] & mphot=mphot[*,a]
  endelse
 endif
 dmphot=2.5/alog(10.)*dfphot/fphot
END
