;----------------------------------------------------------------------
;
; SUBROUTINE TO 1) READ THE .FITS FILES CONTAINING THE FILTERGRAMS
; 2) READ THE DYE LASER WAVELENGTH  
; 3) REBINS THE FILTERGRAMS
; 4) NORMALIZE EACH IMAGE BY THE DYE LASER RELATIVE INTENSITY
; based on a routine provided by J. Schou
;
;----------------------------------------------------------------------


FUNCTION readfile2,header,nx,ny,lam0,base,list
; THE ROUTINE ROTATION.PRO MUST BE COMPILED PRIOR
; TO THIS CODE
; ROTATION.PRO RETURNS THE ROTATION ANGLE DUE TO THE HELIOSTAT
; suninfo PROVIDES THE p-angle

SET_PLOT,'x'
n      = 29 ; numbers of filtergrams: the 1st and last ones are just dark current
nkey   = 91
header = STRARR(n-2,nkey)
wavelen= DBLARR(n-2)
I0     = FLTARR(27)
lam0   = FLTARR(27)
expo   = FLTARR(27)

imx    = FLTARR(nx,ny,n-2)
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
dark1 = REBIN(dark1,nx,ny)
PRINT,'OBSERVATION DATE',day,month,year
PRINT,'OBSERVATION TIME',heure,minute,seconde
PRINT,'EXPOSURE TIME',head[68]

FOR i=1,n-2 DO READF,1,name
READF,1,name
PRINT,'2nd DARK CURRENT= ',name
im     = READFITS(base+name,head,/silent)
dark2  = im[ix,*] + 0.d0        ; to suppress vertical central line
dark2  = dark2[*,iy]
year   = LONG(STRMID(head[8],19,4))
month  = LONG(STRMID(head[8],24,2))
day    = LONG(STRMID(head[8],27,2))
heure  = LONG(STRMID(head[9],21,2)) ; !!!!!!!!!!!!! I ASSUME THE HOUR IS UT !!!!!!!!!!!
minute = LONG(STRMID(head[9],24,2))
seconde= LONG(STRMID(head[9],27,2))
timedark2   = heure*3600.d0 + minute*60.d0 + seconde
dark2 = REBIN(dark2,nx,ny)
PRINT,'OBSERVATION DATE',day,month,year
PRINT,'OBSERVATION TIME',heure,minute,seconde
PRINT,'EXPOSURE TIME',head[68]

CLOSE,1

; WE READ THE FILTERGRAMS
;-----------------------------------------
OPENR,1,list
READF,1,name
FOR i=0,n-3 DO BEGIN
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

  ; dye laser wavelength
    lam0[i]= FLOAT(STRMID(head[73],25,4)) ;*10.d0 ; in Angstrom

    header[i,*] = head

    q      = im[ix,*] + 0.d0   ; to suppress vertical central line
    q      = q[*,iy]           ; to suppress horizontal central line


;    FOR j=0,4095 DO BEGIN
;        dark[0,j]=MEAN( (dark1[0:200,j]+dark2[0:200,j])/2.0)
;        dark[1,j]=MEAN( (dark1[4095-200:4095,j]+dark2[4095-200:4095,j])/2.0)
;        q[0:2047,j]=q[0:2047,j]-dark[0,j]
;        q[2048:4095,j]=q[2048:4095,j]-dark[1,j]
;    ENDFOR

    imx[*,*,i] = REBIN(q,nx,ny); we rebin the filtergram

;   LINEAR INTERPOLATION AS SUGGESTED BY JESPER
    FOR j=0,ny-1 DO BEGIN
        FOR jj=0,nx-1 DO BEGIN
            imx[jj,j,i] = imx[jj,j,i] - INTERPOL([dark1[jj,j],dark2[jj,j]],[timedark1,timedark2],time)
        ENDFOR
    ENDFOR
    temp = REFORM(imx[*,*,i])
    a    = WHERE(temp LT 0.d0);    FOR j=0,4095 DO BEGIN
;        dark[0,j]=MEAN( (dark1[0:200,j]+dark2[0:200,j])/2.0)
;        dark[1,j]=MEAN( (dark1[4095-200:4095,j]+dark2[4095-200:4095,j])/2.0)
;        q[0:2047,j]=q[0:2047,j]-dark[0,j]
;        q[2048:4095,j]=q[2048:4095,j]-dark[1,j]
;    ENDFOR
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
    expo[i] = TOTAL(TOTAL(imx[0:128*nx/4096.,0:128*nx/4096.,i],1),1)
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

SAVE,imx,header,expo,FILE='SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

RETURN,imx

END
