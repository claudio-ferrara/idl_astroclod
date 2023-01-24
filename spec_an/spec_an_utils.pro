function app_vr_correction, w, vrad

  ;--------function app_vr_correction----------
  ;-------this function applies radial velictiy correction to a wavelength w
  ;-------
  ;-------
  ;-------input:
  ;----------w - wavelength or array of wavelengths
  ;----------vrad - radial velocity in km/s
  ;-------output:
  ;----------l' - radial velocity shifted wavelength(s)  (= (1-vrad/vellight)*l)
  ;-------
  ;;---vellight = 299792.458 - value of the velocity of light in km/s
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  return, (1.-vrad[0]/299792.458d)*w

end


function app_vr_correction_list,lam,vrad

  ;--------function app_vr_correction_list----------
  ;-------this function applies radial velictiy correction to a list of wavelengths lam
  ;-------
  ;-------
  ;-------input:
  ;----------lam - list of (array of) wavelengths
  ;----------vrad - radial velocity in km/s
  ;-------output:
  ;----------lam2 - list of radial velocity shifted wavelength(s)  (= (1-vrad/vellight)*lam)
  ;-------
  ;;---vellight = 299792.458 - value of the velocity of light in km/s
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  lam2 = list()
  foreach l,lam do begin
    lam2.add,(1.-vrad[0]/299792.458d)*l
  endforeach

  return, lam2

end

pro totspectrum_mid,lam,flx,wtot,ftot

  ;--------procedure totspectrum_mid----------
  ;-------this procedure generates a "mosaic" spectrum (wtot,ftot) starting
  ;-------from a 2D spectrum (lam,flx)
  ;-------the mosaic is made changing between an order and the other in the mid-point of the superposition region
  ;-------
  ;-------
  ;-------input:
  ;----------lam - 2D array of wavelengths
  ;----------flx - 2D array of fluxes
  ;-------output:
  ;----------wtot - overall array of wavelengths
  ;----------ftot - overall array of fluxes
  ;----------
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  ;-------getting number of orders
  nords = n_elements(lam(0,*))
  ;---sorting orders by initial wavelengths----
  srt = sort(lam(0,*))
  lamm = lam(*,srt)
  flxm = flx(*,srt)

  wtot = lamm(*,0)
  ftot = flxm(*,0)

  ;-----for cycle to compose the mosaic--------
  for i = 1,nords-1 do begin
    wlim = (wtot[-1]+lamm(0,i))/2.
    tmptot = where(wtot lt wlim)
    tmpord = where(lamm(*,i) gt wlim)
    wtot = [wtot[tmptot],lamm(tmpord,i)]
    ftot = [ftot[tmptot],flxm(tmpord,i)]
  endfor

end

pro get_order,lam,flx,kord,wo,fo

  ;--------procedure get_order----------
  ;-------this procedure extracts an order (kord) from a 2D (lam,flx) array
  ;-------input:
  ;----------lam - 2D array of wavelengths
  ;----------flx - 2D array of fluxes
  ;----------kord - number of the order (of the aperture)
  ;-------output:
  ;----------wo - order array of wavelengths
  ;----------fo - order array of fluxes
  ;----------
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  ;-------checking if lam is of type "list"-----
  if isa(lam,"list") then begin
    wo = lam[kord]
    fo = flx[kord]
  endif else begin
    wo = lam(*,kord)
    fo = flx(*,kord)
  endelse

end

function find_order,lam,l,IFFIRST=iffirst

  ;--------function find_order----------
  ;-------this function finds the order for a given wavelength 
  ;-------returns -1 if the wavelength is not covered by the spectrum
  ;-------input:
  ;----------lam - 2D array of wavelengths
  ;----------l - wavelength to be checked for
  ;-------output:
  ;----------n_ord = number of the order "containing" wavelength l
  ;----------
  ;-------optional input parameters:
  ;----------iffirst - set to return only the first order in case of multiple occurrence of l (order superposition)
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  ;--------creating arrays of initials and final wavelengths
  ;--------case of list
  if isa(lam,"list") then begin
    
    ins = list()
    fins = list()
    for i=0,n_elements(lam)-1 do begin
      wo = lam[i]
      ins.add, wo[0]
      fins.add,wo[-1]
    endfor
    ins = ins.toarray()
    fins = fins.toarray()
    
  endif else begin ;----case of vector
    
    ins = lam(0,*)
    fins = lam(-1,*)
    
  endelse

  tmp = where((ins lt l) and (fins gt l))
  
  if n_elements(tmp) gt 1 then begin
    if keyword_set(iffirst) then n_ord = tmp[0] else n_ord = tmp
  endif else n_ord = tmp
    
  return, n_ord

end

pro crop_spectrum,w,f,w1,w2,wc,fc,TMP=tmp

  ;--------procedure crop_spectrum----------
  ;-------this procedure crops a spectrum given a region of (w1,w2) boundary wavelengths
  ;-------input:
  ;----------w - array of wavelengths
  ;----------f - array of fluxes
  ;----------w1 - initial wavelength
  ;----------w2 - final wavelength
  ;-------output:
  ;----------wc - wavelengths of subspectrum
  ;----------fc - fluxes of subspectrum
  ;-------optional output parameters:
  ;----------tmp - indexes array of the subspectrum
  ;-------
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  tmp = where((w gt w1) and (w lt w2))
  if tmp[0] ne -1 then begin
    wc = w[tmp]
    fc = f[tmp]  
  endif else return
  
end

pro conf_lists_2,w1,w2,tvals,inds,inds_2,inds_mult

  ;--------procedure conf_lists_2----------
  ;-------this procedure compares two array of wavelengths (of values in general)
  ;-------selects the indices of the common value according to a given threshold
  ;-------THIS PROCEDURE IS NOT VECTORIZED--------
  ;-------input:
  ;----------w1 - first array of wavelengths
  ;----------w2 - second array of wavelengths
  ;----------tvals - array of threshold values for comparison
  ;-------output:
  ;----------inds - indexes of the common values of w1
  ;----------inds_2 - indexes of the common values of w2
  ;----------inds_mult - indexes of the values with multiple occurrences
  ;-------
  ;-------
  ;-------
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------
  
  
  ;------procedura che identifica i punti in w1 che hanno punti in w2 a meno di tval
  ;----restituisce inds, che sono gli indici degli elementi di w1 che hanno coincidenze in w2
  
  
  n1 = n_elements(w1) ;---number of elements of w1
  if not isarray(tvals) then tvals = make_array(n1,value=tvals)

  inds = 0
  inds_2 = 0
  inds_mult = 0
  ;-------iterating over all values of w1---
  for i=0,n1-1 do begin
    w = w1[i]
    tval = tvals[i]
    dw = abs(w2-w)
    tmp = where(dw lt tval)
    if n_elements(tmp) gt 1 then begin
      inds_mult = [inds_mult,i] ;---storing all the indexes of v1 of values with multiple occurrences in w2
      print,"MULT!"
    endif
    if tmp[0] ne -1 then begin ;---storing the two indexes of the "common" values
      inds = [inds,i*1l]
      inds_2 = [inds_2,tmp[0]]
    endif
  endfor
  if n_elements(inds) gt 1 then begin
    inds = inds[1:-1]
    inds_2 = inds_2[1:-1]
  endif else begin
    inds = -1
    inds_2 = -1
  endelse
  if n_elements(inds_mult) gt 1 then inds_mult = inds_mult[1:-1] else inds_mult = -1
  
end


