pro readimages_jesper,list,images,nim,headers,nbin=nbin,silent=silent,noshow=noshow,raw=raw,getnames=getnames,iover=iover

; names contains the list of images read or to be read
; images contain the images
; headers contain the headers
; nbin gives the desired image size. Must divide into the original image
; size (4200 for raw images, 4096 for ones with the cross removed. No binning
; is done if not set.
; /silent makes the program shut up during reading. Otherwise the image
; number (counting from zero), the time part of the image name, the FSN,
; the exposure time, the focus, wavelength and polarization table indices
; and the mean dark value is printed.
; /noshow makes the program not show thumbnails as it goes
; /raw makes the program return raw 4200^2 images
; /getnames makes the program query for the image name instead of using
; those in names
; iover gives the intensity in the overscan pixels


if (n_elements(silent) eq 0) then silent=0
if (n_elements(noshow) eq 0) then noshow=0
if (n_elements(raw) eq 0) then raw=0
if (n_elements(getnames) eq 0) then getnames=0

;if (getnames ne 0) then begin
;  selectnames,names,nim
;end
n='qqq'
names=STRARR(nim)
OPENR,1,list
FOR i=0,nim-1 DO BEGIN
    READF,1,n
    names[i]=n
ENDFOR
CLOSE,1

SET_PLOT,'x'
!P.MULTI=0
WINDOW,0,retain=2

nim=n_elements(names)
headers=ptrarr(nim)
iover=fltarr(nim)

nx=4096
if (raw ne 0) then nx=4200
if (n_elements(nbin) eq 0) then nbin=nx
nsmall=64
if (raw ne 0) then nsmall=60

if (noshow eq 0) then nshow=(!d.x_size/nsmall)*(!d.y_size/nsmall)

images=fltarr(nbin,nbin,nim)
ix=[2098-2047+indgen(2048),2101+indgen(2048)]
iy=[2067-2047+indgen(2048),2132+indgen(2048)]
iy0=2068+indgen(64)

if (silent eq 0) then print,'   #  Name    FSN   EXP CAL  WL POL   Dark     Over'
for i=0,nim-1 do begin
  name=names(i)
  im=readfits(name,header,/silent)
  q=im(ix,*)
  io=rebin(q(*,iy0)+0.0,1)
  if (raw eq 0) then begin
    if (nx eq nbin) then images(*,*,i)=q(*,iy)+0.0 else images(*,*,i)=rebin(q(*,iy)+0.0,nbin,nbin)
  endif else begin
    if (nx eq nbin) then images(*,*,i)=im+0.0 else images(*,*,i)=rebin(im+0.0,nbin,nbin)
  end
  ln=strlen(names(0))
  sname=strmid(name,ln-10,6)
  iover(i)=io
  if (silent eq 0) then print,i,sname,sxpar(header,'HSQFGSN'),sxpar(header,'HSHIEXP'),sxpar(header,'HCFTID'),sxpar(header,'HWLTID'),sxpar(header,'HPLTID'),rebin(q(0:63,iy(0:63))+0.0,1),io,format='(i4,a7,i6,i6,3i4,2f9.2)'
  headers(i)=ptr_new(header)
  if (noshow eq 0) then tvscl,rebin(images(*,*,i),nsmall,nsmall),i mod nshow
end

end

