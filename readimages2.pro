; program of J. Schou, slightly modified by S. Couvidat
; FOR DATA AFTER 01/01/2007


pro readimages2,namelist,images,imx,headers,nbin=nbin,silent=silent,noshow=noshow,raw=raw,getnames=getnames,iover=iover,pad=pad,fixcenter=fixcenter

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
; /raw makes the program return raw 4200^2 images for RAL
; /raw makes the program return raw 4096^2 images for CIF
; /getnames makes the program query for the image name instead of using
; those in names
; iover gives the intensity in the overscan pixels
; pad will pad images if CIF and overscan (imcfg>2)
; fixcenter will interpolate center pixel if CIF and imcfg=2


; ADDED BY ME
WINDOW,0,RETAIN=2
nim=30;13
n='qqq'
names=STRARR(nim)
OPENR,1,namelist
FOR i=0,nim-1 DO BEGIN
    READF,1,n
    names[i]=n
ENDFOR
CLOSE,1
plpos=LONARR(3,nim)
pldel=plpos
sign=LONARR(3,nim)

if (n_elements(silent) eq 0) then silent=0
if (n_elements(noshow) eq 0) then noshow=0
if (n_elements(raw) eq 0) then raw=0
if (n_elements(getnames) eq 0) then getnames=0
if (n_elements(pad) eq 0) then pad=1
if (n_elements(fixcenter) eq 0) then fixcenter=1

if (getnames ne 0) then begin
  selectnames,names,nim
end

; Use first image to find type
im=readfits(names(0),header,/silent)
nx=sxpar(header,'NAXIS1')
ny=sxpar(header,'NAXIS2')
config=sxpar(header,'CONFIG')
sz=size(config,/type)
type=-1
if (sz eq 7) then begin ; Keyword exists and is of type string
  if (config eq 'RAL') then type=0
  if (config eq 'CIF') then type=1
endif else begin ; Otherwise guess based on size
  if (nx eq 4200) then type=0
  if (nx eq 4096) then type=1
end
if (type lt 0) then begin
  print,'Unknown image type'
  stop
end

nim=n_elements(names)
headers=ptrarr(nim)
iover=fltarr(nim)

if (type eq 0) then begin ; RAL
  nn=4096
  if (raw ne 0) then nn=4200
  if (n_elements(nbin) eq 0) then nbin=nn
  nsmall=64
  if (raw ne 0) then nsmall=60
  
  if (noshow eq 0) then nshow=(!d.x_size/nsmall)*(!d.y_size/nsmall)
  
  images=fltarr(nbin,nbin,nim)
  ix=[2098-2047+indgen(2048),2101+indgen(2048)]
  iy=[2067-2047+indgen(2048),2132+indgen(2048)]
  iy0=2068+indgen(64)
  
  if (silent eq 0) then print,'   #  Name    FSN   EXP CAL  WL POL   Dark     Over'

  jd      = DBLARR(nim)       ; AJOUT SEBASTIEN

  for i=0,nim-1 do begin
    name=names(i)
    im=readfits(name,header,/silent)

    ; AJOUT SEBASTIEN
    date = SXPAR(header,'T_OBS') ; we extract the date and time of the series
    yyyy = LONG(STRMID(date,0,4))
    mm   = LONG(STRMID(date,5,2))
    dd   = LONG(STRMID(date,8,2))
    heure= LONG(STRMID(date,11,2))
    minute= LONG(STRMID(date,14,2))
    second= LONG(STRMID(date,17,2))
    jd[i]= SXPAR(header,'SHS')
   ;jd[i]=86400*JULDAY(mm,dd,yyyy,heure,minute,second)
    PRINT,'date= ',date
   ;PRINT,jd[i],mm,dd,yyyy,heure,minute,second
    IF(sxpar(header,'HSHIEXP') NE 0) THEN exposuretime=sxpar(header,'HSHIEXP')
    plpos[0,i]=sxpar(header,'HWL1POS')
    plpos[1,i]=sxpar(header,'HWL2POS')
    plpos[2,i]=sxpar(header,'HWL4POS')
    pldel[0,i]=sxpar(header,'H0428') ; Commanded delays
    pldel[1,i]=sxpar(header,'H042A')
    pldel[2,i]=sxpar(header,'H042C')
    FSN=sxpar(header,'FSN')
    move=readx('/home/schou/hmi/suntest/mech_vaccal.move')
    sign[0,i]=(move(3,FSN)<1)>(-1)
    sign[1,i]=(move(4,FSN)<1)>(-1)
    sign[2,i]=(move(6,FSN)<1)>(-1)


    q=im(ix,*)
    io=rebin(q(*,iy0)+0.0,1)
    if (raw eq 0) then begin
      if (nn eq nbin) then images(*,*,i)=q(*,iy)+0.0 else images(*,*,i)=rebin(q(*,iy)+0.0,nbin,nbin)
    endif else begin
      if (nn eq nbin) then images(*,*,i)=im+0.0 else images(*,*,i)=rebin(im+0.0,nbin,nbin)
    end
    ln=strlen(names(0))
    sname=strmid(name,ln-10,6)
    iover(i)=io
    if (silent eq 0) then print,i,sname,sxpar(header,'HSQFGSN'),sxpar(header,'HSHIEXP'),sxpar(header,'HCFTID'),sxpar(header,'HWLTID'),sxpar(header,'HPLTID'),rebin(q(0:63,iy(0:63))+0.0,1),io,format='(i4,a7,i6,i6,3i4,2f9.2)'
    headers(i)=ptr_new(header)
    if (noshow eq 0) then tvscl,rebin(images(*,*,i),nsmall,nsmall),i mod nshow
  end
endif else begin ; CIF
  nn=4096
  if (n_elements(nbin) eq 0) then nbin=nn
  nsmall=64
  
  if (noshow eq 0) then nshow=(!d.x_size/nsmall)*(!d.y_size/nsmall)
  
  images=fltarr(nbin,nbin,nim)
  ix=[2098-2047+indgen(2048),2101+indgen(2048)]
  iy=[2067-2047+indgen(2048),2132+indgen(2048)]
  iy0=2068+indgen(64)
  
  if (silent eq 0) then print,'   #  Name    FSN   EXP CAL  WL POL   Dark     Over'
  jd      = DBLARR(nim)         ; AJOUT SEBASTIEN
  for i=0,nim-1 do begin
    name=names(i)
    im=readfits(name,header,/silent)

    ; AJOUT SEBASTIEN
    date = SXPAR(header,'T_OBS') ; we extract the date and time of the series
    yyyy = LONG(STRMID(date,0,4))
    mm   = LONG(STRMID(date,5,2))
    dd   = LONG(STRMID(date,8,2))
    heure= LONG(STRMID(date,11,2))
    minute= LONG(STRMID(date,14,2))
    second= LONG(STRMID(date,17,2))
    jd[i]= SXPAR(header,'SHS')
    ;jd[i]= 86400*JULDAY(mm,dd,yyyy,heure,minute,second)
    PRINT,'date= ',date
    ;PRINT,jd[i],mm,dd,yyyy,heure,minute,second
    IF(sxpar(header,'HSHIEXP') NE 0) THEN exposuretime=sxpar(header,'HSHIEXP')
    plpos[0,i]=sxpar(header,'HWL1POS')
    plpos[1,i]=sxpar(header,'HWL2POS')
    plpos[2,i]=sxpar(header,'HWL4POS')
    pldel[0,i]=sxpar(header,'H0428') ; Commanded delays
    pldel[1,i]=sxpar(header,'H042A')
    pldel[2,i]=sxpar(header,'H042C')
    FSN=sxpar(header,'FSN')
    move=readx('/home/schou/hmi/suntest/mech_vaccal.move')
    sign[0,i]=(move(3,FSN)<1)>(-1)
    sign[1,i]=(move(4,FSN)<1)>(-1)
    sign[2,i]=(move(6,FSN)<1)>(-1)

    if (fixcenter ne 0) then im(2047,2048)=im(2047,2047)
    imcfg=sxpar(header,'H0149')
    if ((imcfg lt 2) or (imcfg gt 4)) then stop
    im=im and 32767
    if (imcfg eq 2) then begin
      io=0.0
    end
    if (imcfg eq 3) then begin
      io=(mean(im(0:2046,2047:2048))+mean(im(2049:4095,2047:2048)))/2
      im(1:2047,*)=im(0:2046,*)
      im(2048:4094,*)=im(2049:4095,*)
      im(*,1:2047)=im(*,0:2046)
      im(*,2048:4094)=im(*,2049:4095)
      if (pad eq 0) then begin
        im(0,*)=-32767
        im(4095,*)=-32767
        im(*,0)=-32767
        im(*,4095)=-32767
      end
    end
    if (imcfg eq 4) then begin
      io=(mean(im(0:2045,2046:2049))+mean(im(2050:4095,2046:2049)))/2
      im(2:2047,*)=im(0:2045,*)
      im(2048:4093,*)=im(2050:4095,*)
      im(*,2:2047)=im(*,0:2045)
      im(*,2048:4093)=im(*,2050:4095)
      if (pad eq 0) then begin
        im(0:1,*)=-32767
        im(4094:4095,*)=-32767
        im(*,0:1)=-32767
        im(*,4094:4095)=-32767
      end
    end

    if (nn eq nbin) then begin
      images(*,*,i)=im+0.0
    endif else begin
      images(*,*,i)=rebin(im+0.0,nbin,nbin)
    end
    ln=strlen(names(0))
    sname=strmid(name,ln-10,6)
    iover(i)=io
    if (silent eq 0) then print,i,sname,sxpar(header,'HSQFGSN'),sxpar(header,'HSHIEXP'),sxpar(header,'HCFTID'),sxpar(header,'HWLTID'),sxpar(header,'HPLTID'),rebin(im(0:63,iy(0:63))+0.0,1),io,format='(i4,a7,i6,i6,3i4,2f9.2)'
    headers(i)=ptr_new(header)
    if (noshow eq 0) then tvscl,rebin(images(*,*,i),nsmall,nsmall),i mod nshow
  end
end

; AJOUT DE SEBASTIEN
imx   = images
nelem = nim
nx    = nbin
; REMOVE THE DARK CURRENT (BASED ON CLEAN.PRO FROM JESPER SCHOU)
jd    = jd-TOTAL(REBIN(jd,1))
wdark = [1,nelem-1]             ; location of the 2 dark frames (WARNING : VALID ONLY FOR DETUNES)
dark  = REBIN(imx[*,*,wdark],nx,nx,1)
c1    = REBIN(imx[0:nx/64-1,0:nx/64-1,*],1,1,nelem)
c2    = REBIN(imx[nx-nx/64:nx-1,0:nx/64-1,*],1,1,nelem)
c3    = REBIN(imx[0:nx/64-1,nx-nx/64:nx-1,*],1,1,nelem)
c4    = REBIN(imx[nx-nx/64:nx-1,nx-nx/64:nx-1,*],1,1,nelem)
cav0  = REFORM(c1+c2+c3+c4)/4
; Try to separate gain variations and dark current variations
pcav  = POLY_FIT(jd,cav0,2)
cavfit= POLY(jd,pcav)
cavres= cav0-cavfit
cav   = cavfit
cav   = cav/TOTAL(REBIN(cav(wdark),1))
FOR i=0,nelem-1 DO imx[*,*,i]=imx[*,*,i]-cav[i]*dark
intin = 6*TOTAL(REBIN(imx,1))

c1    = REBIN(imx[0:nx/64-1,0:nx/64-1,*],1,1,nelem)
c2    = REBIN(imx[nx-nx/64:nx-1,0:nx/64-1,*],1,1,nelem)
c3    = REBIN(imx[0:nx/64-1,nx-nx/64:nx-1,*],1,1,nelem)
c4    = REBIN(imx[nx-nx/64:nx-1,nx-nx/64:nx-1,*],1,1,nelem)
q     = [[c1,c2],[c3,c4]]
imx   = imx - REBIN(q,nx,nx,nelem,/sample)

PRINT,exposuretime
SAVE,exposuretime,imx,pldel,plpos,sign,file='SEQUENCE_'+STRTRIM(namelist,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'


end

