pro my_read_iraf_fits,filename,lam,flx,SILENT=silent

  ;--------PROCEDURE my_read_iraf_fits----------
  ;---------------------------------------------
  ;---------------------------------------------
  ;---------------------------------------------
  ;---procedure to read "multispec" format fits (WAT2 wavelength solution)
  ;---reduced with Image Reduction and Analysis Facility (IRAF)--------- 
  ;---
  ;---input parameters:
  ;-----filename - name of .fits file (including extension) to be read----
  ;---output parameters:
  ;-----lam - nords*npix matrix of wavelengths
  ;-----flx - nords*pix matrix of fluxes
  ;-----where nords is the number of orders and npix is the number of pixels (points) of the spectrum
  ;---optional keywords:
  ;-----silent - to perform "silent" (non verbose) reading of fits file
  ;---last modification - January 16th, 2023
  ;-----
  ;-----
  ;-----
  ;----------------------------------------------

  on_error,2

  keyroot = "WAT2" ;-------root of the keyword for wavelength solution
                   ;------written in iraf format-----

  ;---------reading fits image---------------------
  flx = readfits(filename,h,silent=silent) ;-------flx is the array of fluxes-------
  npix = n_elements(flx(*,0)) ;-----getting number of pixels------
  
  tmp = where(strpos(h,keyroot) ge 0) ;------extracting all the header rows containing wavelength solution
  if tmp[0] eq -1 then begin
    print, "Not a multispec format"
    print, "Returning -1"
    flx = -1
    lam = -1
    return
  endif
  
  h_sol = h[tmp]
  ;---------------deriving wavelength solution information-------
  strhsol = ''
  for i=1,n_elements(h_sol) do strhsol = strhsol +sxpar(h_sol,keyroot+"_"+string(i,format="(I3.3)"))
  
  strhext = strsplit(strhsol,'"',/extract)
  strhext = strhext[1:-1]
  ;------n_elements of strhext is always even------
  strhext = strhext[indgen(n_elements(strhext)/2,start = 0,increment=2)]
  
  
  ;-------------reconstructing arrays of wavelengths------------------
  lam = list()
  for i=0,n_elements(strhext)-1 do begin

    cols = strsplit(strhext[i]," ",/EXTRACT)
    lam.add,dindgen(npix,start=double(cols[3]),increment=double(cols[4]))

  endfor
  lam = transpose(lam.toarray())
  

end
