pro read_iraf_fits,filename,l,f,l0_vect,Dl_vect
tic
;im=readfits(filename+'.fits',h,/silent)
im = readfits(filename,h,/silent)
n_orders=n_elements(im(0,*))
head=strjoin(h)
; k=strpos(head,'MJD-OBS')
;JD=double(strmid(head,k+9,23))+0.5D
;date = head[strpos(head,'DATE_OBS')]
date = strmid(head,strpos(head,'DATE_OBS')+9,26)
if(n_orders eq 1) then begin
 flx=readfits(filename+'.fits',h)
 Npix=n_elements(flx)
 head=strjoin(h)
 k=strpos(head,'CRVAL1')
 l0=double(strmid(head,k+9,22))
 k=strpos(head,'CDELT1')
 Dl=double(strmid(head,k+9,22))

 Npix=findgen(Npix)
 lam=l0+Dl*Npix
endif

if(n_orders gt 1) then begin
imx=im(*,0)
imy=im(0,*)
Nord=n_elements(imy)
Npix=n_elements(imx)
l=dblarr(Npix,Nord)
f=dblarr(Npix,Nord)
Npix=findgen(Npix)
str=''
ord=strarr(Nord)
for i=0,180 do begin
 k=strpos(head,'WAT2')
if (k gt 0) then begin
 str=str+strmid(head,k+11,68)
 head=strmid(head,k+11)
endif
endfor
head=str
 k=strpos(head,'spec1')
 head=strmid(head,k)

l0_vect = 0
Dl_vect = 0
for i=0,Nord-1 do begin
k=strpos(head,'s')
head=strmid(head,k)

k1=strpos(head,'"')
sub_head=strmid(head,k1+1)
k2=strpos(sub_head,'"')
tmp=strmid(head,k1+1,k2-1)
k=strpos(tmp,' ')
tmp=strmid(tmp,k+1)
k=strpos(tmp,' ')
tmp=strmid(tmp,k+1)
k=strpos(tmp,' ')
tmp=strmid(tmp,k+1)
k=strpos(tmp,' ')
l0=strmid(tmp,0,k)*1.D
tmp=strmid(tmp,k+1)
k=strpos(tmp,' ')
Dl=strmid(tmp,0,k)*1.D

head=strmid(head,k+1)

l(*,i)=l0+Dl*Npix
f(*,i)=im(*,i)

l0_vect = [l0_vect,l0]
Dl_vect = [Dl_vect,Dl]
endfor
l0_vect = l0_vect[1:-1]
Dl_vect = Dl_vect[1:-1]

lam=l(*,0)
flx=f(*,0)
for i=0,Nord-1 do lam = [lam,l[*,i]]
for i=0,Nord-1 do flx = [flx,f[*,i]]
endif 
toc
end
