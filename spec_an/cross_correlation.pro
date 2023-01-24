pro crosscor,w,f,fteor,vrad_max,vrad,V=v,CC=cc

  ;--------procedure crosscor----------
  ;------- based on idl "c_correlate" procedure, this procedure/module
  ;------- calculates cross correlation of two identically sampled emission spectra
  ;------- it follows the treatment of "Hilditch - An Introduction to Close Binary Stars"
  ;-------
  ;-------input:
  ;----------w - array of wavelengths (equal for the two spectra)
  ;----------f - first array of fluxes
  ;----------fteor - second array of fluxes
  ;-------output:
  ;----------vrad_max - radial velocity taken from maximum (peak) value
  ;----------vrad - radial velocity taken from the peak fit
  ;-------optional output parameters:
  ;----------v - array of velocities
  ;----------cc - array of cross correlation function
  ;-------
  ;-------last modification: January 20th 2023
  ;
  ;
  ;------------------------------------------

  ;light velocity in km/s
  vellight=299792.458 ;in km/s

  ;-------number of points for the cross correlation function
  n = 700
  enne = n_elements(w)

  ;------array of "lags" - lags for which cc has to be calculated
  lag = dindgen(n+1,start = -n/2)

  ;------shifting the wavelengths to central wavelength
  lmid = (w[0]+w[-1])/2.
  ;------constructing the "log" scale for the wavelengths 
  logl1 = alog(w/lmid)
  samp = dindgen(enne,start = logl1[0],increment=(logl1[-1]-logl1[0])/enne)
  ;-------interpolating the flux values to the log scale of wavelengths--------
  nf = interpol(f,logl1,samp)
  nfteor = interpol(fteor,logl1,samp)
  ;-------calculating cross correlation using "C_CORRELATE"
  cc = c_correlate(nf,nfteor,lag)
  dlog = samp-shift(samp,1)
  dlog = dlog[1]
  ;------converting ccf in velocity--------
  v = vellight*(exp(dlog*lag)-1.)
  ;------fitting peak of ccf-----------
  ccfit = gaussfit(v,cc,coeffs)
  plot,v,cc
  vrad = coeffs[1]
  mmm = max(cc,locmas)
  vrad_max = v[locmas]
  oplot,v,ccfit,col=180
  print, "Vrad_max = ",vrad_max
  print, "Vrad = ",vrad

end

pro calc_crosscorr,wo,fo,wt,ft,vrad,WI=wi,WF=wf,V=v,CC=cc

  ;--------procedure calc_crosscor----------
  ;------- based on idl "c_correlate" procedure, this procedure/module
  ;------- calculates cross correlation of two spectra
  ;------- the procedure takes as input spectra already normalized for the continuum
  ;------- it follows the treatment of "Hilditch - An Introduction to Close Binary Stars"
  ;-------
  ;-------input:
  ;----------wo - array of wavelengths
  ;----------fo - first array of fluxes
  ;----------wt - second array of wavelengths
  ;----------ft - second array of fluxes
  ;-------output:
  ;----------vrad - radial velocity taken from the peak fit
  ;-------optional input parameters:
  ;----------wi - initial wavelength
  ;----------wf - final wavelength
  ;-------optional output parameters:
  ;----------v - array of velocities
  ;----------cc - array of cross correlation function
  ;-------
  ;-------last modification: January 20th 2023
  ;
  ;
  ;------------------------------------------

  if not keyword_set(wi) then wi=wo[0]
  if not keyword_set(wf) then wf=wo[-1]

  ;-------cropping the two spectra-------
  crop_spectrum,wo,fo,wi,wf,woc,foc
  crop_spectrum,wt,ft,wi,wf,wtc,ftc

  ;-------interpolating the second spectrum to the wavelengths of the first
  ftco = spline(wtc,ftc,woc)

  ;-------converting the two spectra in "emission"
  yoc = 1.-foc
  ytc = 1.-ftco
  
  ;-------calculating cross correlation function
  crosscor,woc,yoc,ytc,vrad_max,vrad,V=v,CC=cc


end