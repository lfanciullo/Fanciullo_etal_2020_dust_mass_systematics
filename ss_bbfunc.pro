@rootnr
@ss_zeta
;+
;DESCRIPTION:
;	Return a [modified] blackbody function on a given wavelength[micron] (frequency[Hz]) grid, 
;		for the given temperature Teff[K].
;	The blackbody function B_nu is defined in terms of the dimensionless variable x = hc/(lambda*k*Teff) as
;		B_nu[x] = constant * x^3/[Exp(x) - 1] with constant = 2hc/(hc/(k*Teff))^3
;		Or, equivalently,
;		B_lambda[x] = constant * x^5/[Exp(x) - 1] with constant = 2hc^2/(hc/(k*Teff))^5
;	The mass absorption coefficient (MAC) is assumed to have a power-law dependence on the wavelength:
;		Kappa_lambda = KL * lambda^(-alpha) = KL2 * x^(alpha), with KL2 = KL * (hc/(k*Teff))^(-alpha)
;		Typically, 1 < alpha < 2
;		The normalisation KL depends on the value of Kappa at some reference wavelength. This needs to be
;			determined independently outside the program, so Kappa_lambda is fed into the modified
;			blackbody function in dimensionless form.
;	Therefore, the modified blackbody function output by the code is:
;		Bmod_nu[x] = constant * x^(3+alpha)/[Exp(x)-1] with constant = 2hc/(hc/(k*Teff))^3
;		Or, equivalently,
;		Bmod_lambda[x] = constant x^(5+alpha)/[Exp(x) - 1] with constant = 2hc^2/(hc/(k*Teff))^5
;HISTORY
;Written December 6 2011/Updated 2013-10-07 by Sundar Srinivasan
;Updated 2013-10-07: Added ALPHA and SCALED keywords
;Updated 2013-12-13: Added LUMINOSITY keyword
;Updated 2015-08-01: LUMINOSITY keyword changed to POWER
;USAGE
;FUNCTION SS_BBFUNC,xvar,teff,freq=freq,flam=flam,intensity=intensity,nufnu=nufnu,$
;	alpha=alpha,scaled=scaled,power=power,help=help
;	XVAR - by default, the wavelength in micron. If the keyword FREQ is set, XVAR must be the frequency in Hz.
;	TEFF - the temperature in K.
;	FREQ - if set, XVAR is the frequency in Hz.
;	ALPHA - power law index for the mass-absorption coefficient, such that Kappa_lambda ~ lambda^(-alpha). Defaults to zero.
;	FLAM - by default, the program returns a frequency-specific function. Setting FLAM outputs a wavelength-specific function
;		instead; i.e., ~x^(5+alpha) instead of ~x^(3+alpha)
;	INTENSITY - by default, the output is the flux in SI units. Setting this keyword returns the intensity in SI units instead.
;	NUFNU - result in nu*F_nu (SI) units.
;		Note1: To convert nu*F_nu into nu*L_nu in Lsun, multiply by (Rstar/Rsun)^2.
;		Note2: If both INTENSITY and NUFNU are set, output is nu*I_nu.
;		Note3: If FLAM is set, NUFNU is ignored.
;	SCALED - the peak of the blackbody function is normalised to unity.
;		Note1: the actual range of the output depends on how well the input wavelength grid samples the peak position of the blackbody function.
;			For x = hc/(k*Teff*lambda), x_max for B_nu = 2.82143937212, x_max for B_lambda = 4.965114231744276.
;		Note2: To compute the position of the peak value, the script solves for the Lambert W function using the Newton-Raphson script, ROOTNR.
;	POWER - The [modified] blackbody function integrated over the entire frequency/wavelength range. For alpha=0, the power radiated by the blackbody in W/m^2.
;	HELP - fires helpful bullets just like the Armoured Scorpion of Death (https://vimeo.com/64522190).
;-
FUNCTION solveLambertW,x,z
 return,x*exp(x)-z 
END

FUNCTION SS_BBFUNC,xvar,teff,freq=freq,flam=flam,intensity=intensity,nufnu=nufnu,$
        help=help,alpha=alpha,scaled=scaled,power=power
 if keyword_set(help) then begin
  doc_library,'bbfunc'
  return,0.
 endif else begin
  if n_elements(alpha) eq 0 then alpha=0d
  ;Some fundamental constant definitions
  doublepi = double(!pi) & speedoflight = 299792458d & plancksconstant = 6.62607004d-34 & kboltzmann = 1.38064853d-23
  hcbykSI=plancksconstant*speedoflight/kboltzmann & hcbykmicron = hcbykSI*1d6
  const=2*speedoflight*plancksconstant/(hcbykSI/teff)^3
  
  if keyword_set(freq) then intvar=1d6*speedoflight/xvar else intvar=xvar

  x=hcbykmicron/intvar/teff & bx=x^(3+alpha)/(exp(x)-1d) ;(3+alpha) to reflect a lambda^(-alpha) modification to the BB
  bnu=bx
  out=bnu

  p=0.
  if keyword_set(nufnu) then begin
   str='/nufnu'
   out*=x & const*=speedoflight/(hcbykSI/teff);new edit
   p=1.
  endif
  
  str=''
  if keyword_set(flam) then begin
   str='/flam'
   out*=x^2 & const*=speedoflight/(hcbykSI/teff)^2.
   p=2.
  endif

  power=doublepi*const*speedoflight/(hcbykSI/teff)*gamma(4+alpha+p)*ss_zeta(4+alpha+p) ;extra factors from converting \nu to x
  if keyword_set(intensity) then begin
   str+=',/intensity'
  endif else begin
   const*=doublepi
  end

  rootnr,'solveLambertW',w,exp(-1d),-(3+alpha+p)/exp(3+alpha+p)

  if ~keyword_set(scaled) then begin
   out*=const
  endif else begin
   xmax=w+3+alpha+p & lammax=hcbykmicron/teff/xmax & junk=execute("outmax=bbfunc(lammax,teff"+str+",alpha="+strtrim(string(alpha),2)+")")
   out*=const/outmax
  endelse

  return,out
 endelse
END
