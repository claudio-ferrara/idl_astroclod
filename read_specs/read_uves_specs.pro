pro read_uves_spec,filename,w,f,ef,SILENT=silent

  ;--------procedure read_uves_spec----------
  ;-------procedure which reads UVES reduced spectra from eso science archive (ADP.)
  ;-------input:
  ;----------filename - name of fits file to be read
  ;-------output:
  ;----------w - array of wavelengths
  ;----------f - array of fluxes
  ;----------ef - array of errors of fluxes
  ;-------optional input parameters:
  ;----------silent - keyword to perform "silent" (non verbose) reading of fits file
  ;-------------
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  im = mrdfits(filename,1,SILENT=silent) ;------reading fits file-----
  
  w = im.WAVE
  f = im.FLUX_REDUCED
  ef = im.ERR_REDUCED

end