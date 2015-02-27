FUNCTION readfile_lindop,header,nx,base,list
; THE ROUTINE ROTATION.PRO MUST BE COMPILED PRIOR
; TO THIS CODE
; ROTATION.PRO RETURNS THE ROTATION ANGLE DUE TO THE HELIOSTAT
; suninfo PROVIDES THE p-angle

SET_PLOT,'x'
n      = 5 ; numbers of filtergrams: the 1st and last ones are just dark current
ndop   = 17 ; number of dopplergrams
nkey   = 91
header = STRARR(n*ndop,nkey)
I0     = FLTARR(ndop*n)
expo   = FLTARR(ndop*n)

imx    = FLTARR(nx,nx,ndop*n)
name   = 'qqq'
ix     = [2098-2047+INDGEN(2048),2101+INDGEN(2048)]
iy     = [2067-2047+INDGEN(2048),2132+INDGEN(2048)]

WINDOW,0,retain=2,xsize=800,ysize=600
!p.multi=0

; WE READ THE 2 DARK CURRENTS
;-----------------------------------------
OPENR,1,list
READF,1,name
PRINT,'1st DARK CURRENT= ',name
im     = READFITS(base+name,head,/silent)
dark1  = im[ix,*] + 0.d0        ; to suppress vertical central line
dark1  = dark1[*,iy]                ; to suppress horizontal central line
year   = LONG(STRMID(head[8],19,4))
month  = LONG(STRMID(head[8],24,2))
day    = LONG(STRMID(head[8],27,2))
heure  = LONG(STRMID(head[9],21,2)) ; !!!!!!!!!!!!! I ASSUME THE HOUR IS UT !!!!!!!!!!!
minute = LONG(STRMID(head[9],24,2))
seconde= LONG(STRMID(head[9],27,2))
timedark1   = heure*3600.d0 + minute*60.d0 + seconde
dark1 = REBIN(dark1,nx,nx)
PRINT,'OBSERVATION DATE',day,month,year
PRINT,'OBSERVATION TIME',heure,minute,seconde
PRINT,'EXPOSURE TIME',head[68]

; READ A CALMODE IMAGE
;-----------------------------------------

READF,1,name

; WE READ THE FILTERGRAMS
;-----------------------------------------
READF,1,name
FOR i=0,n*ndop-1 DO BEGIN
    READF,1,name
    PRINT,' '
    PRINT,'FILE READ= ',name
    im = READFITS(base+name,head,/silent)

    year   = LONG(STRMID(head[8],19,4))
    month  = LONG(STRMID(head[8],24,2))
    day    = LONG(STRMID(head[8],27,2))
    heure  = LONG(STRMID(head[9],21,2)) ; !!!!!!!!!!!!! I ASSUME THE HOUR IS UT !!!!!!!!!!!
    minute = LONG(STRMID(head[9],24,2))
    seconde= LONG(STRMID(head[9],27,2))
    time   = heure*3600.d0 + minute*60.d0 + seconde

    header[i,*] = head

    q      = im[ix,*] + 0.d0   ; to suppress vertical central line
    q      = q[*,iy]           ; to suppress horizontal central line

    imx[*,*,i] = REBIN(q,nx,nx); we rebin the filtergram

    imx[*,*,i] = imx[*,*,i]-dark1

    temp = REFORM(imx[*,*,i])
    a = WHERE(temp LT 0.d0)
    IF (a[0] NE -1) THEN temp[a]= 0.d0
    imx[*,*,i]=temp


    TVIM,imx[*,*,i],/scale

  ; we read the dye laser intensity
    I0[i]  = FLOAT(STRMID(head[74],26,3))

    PRINT,'OBSERVATION DATE',day,month,year
    PRINT,'OBSERVATION TIME',heure,minute,seconde
    PRINT,'LASER WAVELENGTH',head[73]
    PRINT,'LASER INTENSITY',head[74]
    PRINT,'EXPOSURE TIME',head[68]
    expo[i] = TOTAL(TOTAL(imx[0:128.*nx/4096.,0:128.*nx/4096.,i],1),1)
    PRINT,'BAD EXPOSURE ?',expo[i] ; we average intensity over a small square in the lower left corner
    
ENDFOR
CLOSE,1

;  I NORMALIZE EACH IMAGE BY THE RELATIVE LASER INTENSITY
;norm = MEAN(I0)
;FOR i=0,26 DO imx[*,*,i] = imx[*,*,i] / I0[i]*norm

; WE FREE MEMORY
dark1 = 0.d0
dark2 = 0.d0
q     = 0.d0
im    = 0.d0

SAVE,imx,header,expo,FILE='SEQUENCE_LINDOP'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

RETURN,imx

END
