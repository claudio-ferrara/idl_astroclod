pro read_espadons_simple,filename,w,f,ef,SILENT=silent,UNNORM=unnorm

  ;--------procedure read_espadons_simple----------
  ;-------procedure which reads ESPADONS *i* format spectra
  ;-------this procedure doesn not reorder wavelengths from the mosaic
  ;-------wavelengths are in NANOMETERS
  ;-------input:
  ;----------filename - name of fits file to be read
  ;-------output:
  ;----------w - array of wavelengths
  ;----------f - array of fluxes
  ;----------ef - array of flux errors
  ;-------optional input parameters:
  ;----------silent - keyword to perform "silent" (non verbose) reading of fits file
  ;----------unnorm - if set, reads the unnormalized spectrum from the file, default is normalized
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  im = readfits(filename,SILENT=silent)
  
  if keyword_set(unnorm) then begin
    w = im(*,3)
    f = im(*,4)
    ef = im(*,5)  
  endif else begin
    w = im(*,0)
    f = im(*,1)
    ef = im(*,2)  
  endelse
  
end


pro read_espadons_orders,filename,lam,flx,wtot,ftot,IFAVG=ifavg,UNNORM=unnorm,SILENT=silent,E_FLX=e_flx,EFTOT=eftot


  ;--------procedure read_espadons_orders----------
  ;-------procedure which reads ESPADONS *i* format spectra
  ;-------gives as output two lists of orders (lam,flx)
  ;-------and a mosaic composed with a cut or average
  ;-------wavelengths are in NANOMETERS
  ;-------input:
  ;----------filename - name of fits file to be read
  ;-------output:
  ;----------lam - list of wavelength arrays
  ;----------flx - list of flux arrays
  ;----------wtot - composed array of wavelengths
  ;----------ftot - composed array of fluxes
  ;-------optional input parameters:
  ;----------ifavg - if set constructs the mosaic taking an average on the common regions
  ;----------      - otherwise is a simple mosaic
  ;----------silent - keyword to perform "silent" (non verbose) reading of fits file
  ;----------unnorm - if set, reads the unnormalized spectrum from the file, default is normalized
  ;-------optional output parameters:
  ;----------e_flx - list of flux errors
  ;----------eftot - composed array of flux errors
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------

  on_error, 2

  read_espadons_simple,filename,w,f,ef,SILENT=silent,UNNORM=unnorm
  
  ;-----finding the superimposition regions-----
  d_w = w-shift(w,1)
  tmp = where(d_w lt -1 or d_w gt 1)
  
  lam = list()
  flx = list()
  e_flx = list()
  
  for i=1,n_elements(tmp)-1 do begin
    
    word = w[tmp[i-1]:tmp[i]-1]
    ford = f[tmp[i-1]:tmp[i]-1]
    eford = ef[tmp[i-1]:tmp[i]-1]

    lam.add,word
    flx.add,ford
    e_flx.add,eford

    ;------------------------------------------------------
    ;---------merging with average in common regions---------
    if keyword_set(ifavg) then begin
      if i eq 1 then begin
        wtot = word
        ftot = ford
        eftot = eford
      endif else begin
        tmp1 = where(word lt lam(i-2,-1),complement=ntmp1)
        tmp2 = where(lam(i-2) gt word[0],complement=ntmp2)
        if (tmp1[0] ge 0 and tmp2[0] ge 0) then begin
          
          wavg = word[tmp1]
          f2_1 = spline(lam(i-2,tmp2),flx(i-2,tmp2),word[tmp1])
          favg = (ford[tmp1]+f2_1)/2.
          word2 = word
          ford2 = ford
          eford2 = eford
          ford2[tmp1] = favg
          tmpf = where(wtot lt word2[0])
          wtot = [wtot[tmpf],word2]
          ftot = [ftot[tmpf],ford2]
          eftot = [eftot[tmpf],eford2]

        endif else begin
          wtot = [wtot,word]
          ftot = [ftot,ford]
          eftot = [eftot,eford]
        endelse
      endelse
    endif else begin
      ;--------merging composing a single mosaic--------------------
      if i eq 1 then begin
        wtot = word
        ftot = ford
        eftot = eford
      endif else begin
        tmpp = where(word lt lam(i-2,-1),complement=ntmpp)
        ;print, tmpp
        if tmpp[0] ge 0 then begin
          wtot = [wtot,word[ntmpp]]
          ftot = [ftot,ford[ntmpp]]
          eftot = [eftot,eford[ntmpp]]
          ;stop
        endif else begin
          wtot = [wtot,word]
          ftot = [ftot,ford]
          eftot = [eftot,eford]
        endelse
      endelse
    endelse
    ;---------------------------
  endfor


end

