function extract_head_val,filename,value

  ;--------function extract_head_val----------
  ;-------procedure which gets a keyword from a fits file header
  ;-------it works with header files in which the comment starts with '/'
  ;-------if the keyword is found twice or more returns the first occurrence-----
  ;-------input:
  ;----------filename - name of fits file to be read for the header
  ;-------output:
  ;----------value: keyword value 
  ;-------optional input parameters:
  ;-------------
  ;-------------
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------
  
  on_error, 2

  header = headfits(filename) ;------reading the fits header
  ;-------searching the keyword in header--------
  ind = strpos(header,value)
  tmp = where(ind ge 0)

  ;------returns "-1" if the keyword is not found------------
  ;------if the keyword is found twice or more returns the first occurrence-----
  if tmp[0] ne -1 then valstr = header[tmp[0]] else return,"-1"
  i1 = strpos(valstr,"=")
  i2 = strpos(valstr,"/")

  val = strtrim(strmid(valstr,i1+1,i2-i1-1),2)

  return, val

end


pro read_harpn_e2ds,filename,lam,flx,SILENT=silent

  ;--------procedure read_harpn_e2ds----------
  ;-------procedure which reads HARPS-N *e2ds* 2D format spectra
  ;-------input:
  ;----------filename - name of fits file to be read
  ;-------output:
  ;----------lam - 2D array of wavelengths
  ;----------flx - 2D array of fluxes
  ;-------optional input parameters:
  ;----------silent - keyword to perform "silent" (non verbose) reading of fits file
  ;-------------
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------
  
  on_error, 2


  flx = readfits(filename,h,SILENT=silent) ;------reading fits file------

  n = n_elements(flx(0,*))
  npix = n_elements(flx(*,0))

  ;-------------reconstructing wavelength solution--------------
  deg = extract_head_val(filename,"CAL TH DEG LL")
  deg = fix(deg[0])

  coeff_key = "CAL TH COEFF LL"
  ind = strpos(h,coeff_key)
  tmp = where(ind ge 0)
  h_coeff = h[tmp]

  m = n_elements(h_coeff)
  coeffs = 0
  coeff_numbers = 0
  
  for i=0,m-1 do begin
    
    hstring = h_coeff[i]
    i1 = strpos(hstring,"=")
    i2 = strpos(hstring,"/")
    val = strtrim(strmid(hstring,i1+1,i2-i1-1),2)
    cof = double(val)
    coeffs = [coeffs,cof]
    nstring = strmid(hstring,0,i1)

    coeff_numbers = [coeff_numbers,fix(strmid(nstring,nstring.IndexOf("LL")+2))]

  endfor
  
  coeffs = coeffs[1:-1]
  coeff_numbers = coeff_numbers[1:-1]
  srt = sort(coeff_numbers)
  coeffs = coeffs[srt]
  x = indgen(npix)

  ;-------constructing array of wavelengths----------------
  lam = fltarr(npix)

  for i=0,m-1,deg+1 do begin
    lam = [[lam],[poly(x,coeffs[i:i+deg])]]
  endfor
  
  lam = lam(*,1:-1)

end

pro read_harpn_e2ds_filelist,filelist,lams,flxs,SILENT=silent

  ;--------procedure read_harpn_e2ds_filelist----------
  ;-------procedure which reads filelists of HARPS-N *e2ds* 2D format spectra
  ;-------input:
  ;----------filelist - list of fits files to be read
  ;-------output:
  ;----------lams - 3D array of wavelengths (npix,nords,nfiles)
  ;----------flxs - 3D array of fluxes (npix,nords,nfiles)
  ;-------optional input parameters:
  ;----------silent - keyword to perform "silent" (non verbose) reading of fits file
  ;-------------
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  nfiles = n_elements(filelist)

  lams = list()
  flxs = list()
  for i=0,nfiles-1 do begin

    read_harpn_e2ds,filelist[i],lam2,flx2,silent=silent
    lams.add,lam2
    flxs.add,flx2

  endfor
  
  lams = transpose(lams.toarray(),[1,2,0])
  flxs = transpose(flxs.toarray(),[1,2,0])

end

pro read_harpn_s1d,filename,w,f,SILENT=silent

  ;--------procedure read_harpn_s1d----------
  ;-------procedure which reads HARPS-N *s1d* format spectra
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


