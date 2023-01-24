pro oplot_lines,wlt,XR=xr,COL=col,LINESTYLE=linestyle,ORIENTATION=orientation, $ 
                    LINES_LABEL=lines_label,LABEL_HEIGHT=label_height, $ 
                    D_OFFS=d_offs

  
  ;--------procedure oplot_lines----------
  ;-------this procedure "over-plots" (oplots) a set of spectral lines on a 
  ;-------high resolution (typically stellar) spectrum
  ;-------
  ;-------
  ;-------input:
  ;----------wlt: array of wavelengths to be plotted on the spectrum
  ;-------optional input parameters:
  ;----------col, linestyle - go directly to "oplot"
  ;----------xr - range for the plot, it selects the wavelengths to be plotted
  ;----------orientation - goes directly to xyouts, sets the rotation angle of the labels
  ;----------lines_label - label to be plotted next to spectral lines marks
  ;----------label_height - height of the label in data units
  ;----------d_offs - offset (in data units) to alternate label heights of consecutive labels
  ;-------last modification: January 16th 2023
  ;
  ;
  ;------------------------------------------


    
  if keyword_set(xr) then dw_label = 0.02*(xr[-1]-xr[0]) else dw_label = 0.1

  if not keyword_set(orientation) then orientation = 0
  if not keyword_set(label_height) then label_height=1.1
  if not keyword_set(d_offs) then d_offs = label_height*0.1


  if keyword_set(xr) then begin
    
    tmp = where(wlt gt xr[0] and wlt lt xr[1])
    wl = wlt[tmp]   
    if keyword_set(lines_label) then lines_label_c = lines_label[tmp]
 
  endif else begin
    
    wl = wlt
    if keyword_set(lines_label) then lines_label_c = lines_label
  
  endelse


  for i=0,n_elements(wl)-1 do begin
    
    oplot,[wl[i],wl[i]],[-1d7,1d7],col=col,linestyle=linestyle
    
    if keyword_set(lines_label) then begin
      if keyword_set(offs) then begin
        lablh = label_height +(-1)^i*d_offs
      endif else begin
        lablh = label_height
      endelse
      xyouts,wl[i]+dw_label,lablh,lines_label_c[i],col=col,ORIENTATION=orientation
    endif
    
  endfor

end
