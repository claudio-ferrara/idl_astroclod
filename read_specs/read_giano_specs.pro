pro read_giano_s1d,filename,w,f,SILENT=silent

  ;--------procedure read_giano_s1d----------
  ;-------procedure which reads GIANO *s1d* format spectra
  ;-------input:
  ;----------filename - name of fits file to be read
  ;-------output:
  ;----------w - array of wavelengths
  ;----------f - array of fluxes
  ;-------optional input parameters:
  ;----------silent - keyword to perform "silent" (non verbose) reading of fits file
  ;-------------
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  f = readfits(filename,h,silent=silent)
  
  w0 = sxpar(h,"CRVAL1")
  dw = sxpar(h,"CDELT1")
  
  w = dindgen(n_elements(f),start=w0,increment=dw)

end


pro read_giano_ms1d,filename,lam,flx,snr,ords,SILENT=silent

  ;--------procedure read_giano_ms1d----------
  ;-------procedure which reads GIANO *ms1d* format spectra
  ;-------this procedure reads ms1d spectra and sorts the wavelength arrays
  ;-------input:
  ;----------filename - name of fits file to be read
  ;-------output:
  ;----------lam - 2D array of wavelengths
  ;----------flx - 2D array of fluxes
  ;----------snr - 2D array of snr
  ;----------ords - number of orders
  ;-------optional input parameters:
  ;----------silent - keyword to perform "silent" (non verbose) reading of fits file
  ;-------------
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  im = mrdfits(filename,1,SILENT=silent)
  
  lam = im.wave
  flx = im.flux
  snr = im.snr
  ords = im.order
  
  for i = 0, n_elements(lam(0,*))-1 do begin
    srt = sort(lam(*,i))
    lam(*,i) = lam(srt,i)
    flx(*,i) = flx(srt,i)
    snr(*,i) = snr(srt,i)
  endfor

end


pro read_half_giano_ms1d,filename,lamh,flxh,snrh,ordsh,SILENT=silent

  ;--------procedure read_half_giano_ms1d----------
  ;-------procedure which reads GIANO *ms1d* format spectra
  ;-------splits the orders at mid-point to take into account different sensor areas
  ;-------this procedure reads ms1d spectra and sorts the wavelength arrays
  ;-------input:
  ;----------filename - name of fits file to be read
  ;-------output:
  ;----------lamh - 2D array of wavelengths
  ;----------flxh - 2D array of fluxes
  ;----------snrh - 2D array of snr
  ;----------ordsh - number of orders
  ;-------optional input parameters:
  ;----------silent - keyword to perform "silent" (non verbose) reading of fits file
  ;-------------
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  im = mrdfits(filename,1,silent=silent)
  lam = im.wave
  flx = im.flux
  snr = im.snr
  ords = im.order

  lamh = list()
  flxh = list()
  snrh = list()
  ordsh = list()

  nords = n_elements(lam(0,*)) ;-----getting number of orders-----
  npix = n_elements(lam(*,0)) ;------getting number of pixels-----
  hn = npix/2 

  for i=0,nords-1 do begin

    ;------sorting wavelengths--------
    srt = sort(lam(*,i))
    lam(*,i) = lam(srt,i)
    flx(*,i) = flx(srt,i)
    snr(*,i) = snr(srt,i)

    ;------taking the first half------
    wo1 = lam(0:hn-1,i)
    fo1 = flx(0:hn-1,i)
    snro1 = snr(0:hn-1,i)
    ordso1 = ords[i]

    lamh.add,wo1
    flxh.add,fo1
    snrh.add,snro1
    ordsh.add,ordso1
    
    ;------taking the second half------
    wo2 = lam(hn:-1,i)
    fo2 = flx(hn:-1,i)
    snro2 = snr(hn:-1,i)
    ordso2 = ords[i]

    lamh.add,wo2
    flxh.add,fo2
    snrh.add,snro2
    ordsh.add,ordso2

  endfor

  ;------transposing arrays to standardize-----
  lamh = transpose(lamh.toarray())
  flxh = transpose(flxh.toarray())
  snrh = transpose(snrh.toarray())
  ordsh = ordsh.toarray()

end

