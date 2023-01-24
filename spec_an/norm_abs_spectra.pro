function rej_out,y,fact,TMP=tmp

  ;--------function rej_out----------
  ;-------this function exclude outliers from an array doing a simple "fact" sigma clipping
  ;-------
  ;-------input:
  ;----------y - 1D array of values
  ;----------fact - number of sigma for clipping
  ;-------output:
  ;----------y_clip (y[tmp]) = "sigma-clipped" array
  ;----------
  ;-------optional output parameters:
  ;----------tmp - indexes of "good" values 
  ;-------last modification: January 19th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  med = mean(y) ;-----calculating mean----
  sig = stddev(y) ;------calculating standard deviation----
  
  condition = (y gt med - fact*sig) and (y lt med + fact*sig) ;-----sigma clipping condition----
  tmp = where(condition)
  
  return, y[tmp]

end

pro reject_outliers, y,ynew,NITER=niter,TMP0=tmp0,FACT=fact


  ;--------procedure reject_outliers----------
  ;-------this procedure performs sigma clipping outlier rejection using mean "niter" times
  ;-------input:
  ;----------y - array of values
  ;-------output:
  ;----------ynew - array of "good" values
  ;-------optional input parameters
  ;----------NITER - number of iteration for sigma clipping
  ;-------optional output parameters:
  ;----------tmp - indexes of "good" values
  ;-------
  ;-------last modification: January 19th 2023
  ;
  ;
  ;------------------------------------------
  
  on_error, 2


  ;-------setting default values for optional input parameters------
  if not keyword_set(fact) then fact = 3.
  if not keyword_set(niter) then niter=5
  
  ynew = y
  ynew = rej_out(ynew,fact,TMP=tmp0)
  ;----------cycle for sigma clipping--------
  for i = 1, niter-1 do begin

    ynew = rej_out(ynew,fact,TMP=tmp)
    tmp0 = tmp0[tmp]

  endfor

end

;------------------------COMPLETARE---------------------------
pro norm_order,wocc,focc,fn,deg,coeffs,IFPLOT=ifplot,IFSTOP=ifstop,TMP=tmp,FACT=fact,N_ITER=niter,NITER_FIT=niter_fit

  ;--------procedure norm_order----------
  ;-------this procedure performs absorption spectra normalization for the continuum
  ;-------using "simple sigma clipping" on derivatives and lower sigma clipping to the spectrum
  ;-------to identify and fit continuum
  ;-------input:
  ;----------wocc - array of wavelengths (or equivalent indipendent variable)
  ;----------focc - array of fluxes
  ;----------deg - degree of polynomial for continuum fit (default = 6)
  ;-------output:
  ;----------fn - array containing normalized flux
  ;----------coeffs - coefficient of the fitted polynomial
  ;-------optional input parameters
  ;----------IFPLOT - set this keyword to plot the spectrum with fitted continuum
  ;----------IFSTOP - set this keyword to stop at the end of the procedure
  ;----------FACT - number of sigmas for derivative sigma clipping
  ;----------NITER - number of iteration for derivative sigma clipping
  ;----------NITER_FIT - number of iteration for fit residuals sigma clipping
  ;-------optional output parameters:
  ;----------tmp - indexes of "good" values
  ;-------
  ;-------last modification: January 20th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  ;-------setting default values for input parameters
  if n_elements(deg) eq 0 then deg = 6 ;-----usually a value between 5 and 7 is reasonable
  if not keyword_set(fact) then fact=1. ;----number of sigma for sigma clipping of derivatives--
  if not keyword_set(n_iter) then n_iter=3 
  if not keyword_set(niter_fit) then niter_fit=2
  ;------DARE IN INPUT UNO SPETTRO SENZA NAN-----------------------------
  ;--------COME IDENTIFICARE IL CONTINUO PER NORMALIZZARE: DERIVATA---------

  ;--------excluding possible NaNs---------------
  mask = finite(focc,/NAN)
  tmp = where(mask eq 0,complement=ntmp)
  if tmp[0] ge 0 then begin
    woc = wocc[tmp]
    foc = focc[tmp]
  endif else begin
    fn=focc
    return
  endelse
  ;------------calculatin derivative with IDL "DERIV" function-----
  FD = DERIV(FOC)
  
  ;------------rejecting outliers of derivatives (selecting "strongly non zero" values)---
  REJECT_OUTLIERS,FD,FDNEW,TMP0=TMP,FACT=fact,NITER=n_iter
  ;--------fitting WOC[TMP] ED FOC[TMP] with a polynomial-----
  ;------If a lot of points are available also an interpolation would be possible but is risky------

  
  wcont = woc[tmp]
  fcont = foc[tmp]
  ;---------subtracting the central value to the array of wavelengths---
  ;---------this is done in order to make the fits with "lighter" numbers-----
  ;xbar = mean(wcont)
  xbar = (woc[-1]+woc[0])/2d
  xcont = wcont-xbar

  ;-------calculating the first fit-------
  ;-------hopefully "robust_poly_fit" is more robust-------
  coeffs = robust_poly_fit(xcont,fcont,deg,fcont_fit)

  for i=1,niter_fit-1 do begin
    ;------since we assume that this is an absorption spectrum
    ;------we exclude all the "lower outliers"------
    ;------lower outliers defined as points with residual lower than zero-----
    resfit = fcont-fcont_fit
    tmp = where(resfit gt 0)
    xcont = xcont[tmp]
    fcont = fcont[tmp]
    ;------performing again the fit--------
    coeffs = robust_poly_fit(xcont,fcont,deg,fcont_fit)
    ;print, tmp
  endfor


  ;------deriving fitted continuum----------
  foc_cont = poly(woc-xbar,coeffs)
  focc_cont = poly(wocc-xbar,coeffs)
  if keyword_set(ifplot) then begin
    plot,woc,foc
    oplot,woc,foc_cont,col=180
  endif
  
  ;--------normalizing--------
  fn = focc/focc_cont

  if keyword_set(ifstop) then stop
  
end

