FUNCTION MADM,X,DIMENSION=dimension,STDDEV=STDDEV,MEDIAN=median
;+Rewritten 2016-05-08 for 2+dimensional arrays.
; Given an array X, return the median absolute deviation from the median (MADM) along dimension DIMENSION.
;   if keyword STDDEV is set, then return 1.4826*MADM, which is a robust estimate for the STDDEV (exact for Gaussian).
;   if keyword MEDIAN is set, also return the median in the variable MEDIAN.
 s=size(x)
 if n_elements(dimension) ne 0 then begin
  median=median(x,dimension=dimension,/double,/even)
  if dimension eq 1 then p1=1. else p1=product(s[1:dimension-1])
  if dimension ge s[0] then p2=1. else p2=product(s[dimension+1:s[0]])
  dimlist=strjoin(strtrim(string(s[1:s[0]]),2),',')
  ;p1=product(s[1:dimension-1]) & p2=product(s[dimension+1:s[0]]) & dimlist=strjoin(strtrim(string(s[1:s[0]]),2),',')

  junk=execute("xxmed=reform(transpose(reform(cmreplicate(transpose(reform(median,p1,p2)),s[dimension]),p2,p1*s[dimension])),"+dimlist+")")
  madm_=median(abs(x-xxmed),/double,/even,dimension=dimension)
 endif else begin
  median=median(x,/double,/even) & madm_=median(abs(x-median),/double,/even)
 endelse
 if keyword_set(stddev) then return,1.4826*madm_ else return,madm_
END
