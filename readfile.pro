;----------------------------------------------------------------------
;
; SUBROUTINE TO 1) READ THE .FITS FILES CONTAINING THE FILTERGRAMS
; 2) ROTATE THE IMAGES TO ALIGN THE SOLAR ROTATION AXIS
; 3) DETERMINE THE SOLAR Fe I LINE CENTRAL WAVELENGTH
;    THIS WAVELENGTH IS ASSUMED TO BE CONSTANT DURING THE ENTIRE
;    DETUNE SEQUENCE    

; THIS ROUTINE CORRECTS FOR BAD EXPOSURES 
; THIS ROUTINE INDICATES THE OVERSCANS

; based on a routine provided by J. Schou
;
;----------------------------------------------------------------------


FUNCTION readfile,nx,ny,lam0,list
; THE ROUTINE ROTATION.PRO MUST BE COMPILED PRIOR
; TO THIS CODE
; ROTATION.PRO RETURNS THE ROTATION ANGLE DUE TO THE HELIOSTAT
; suninfo PROVIDES THE p-angle

n      = 29 ; numbers of filtergrams: the 1st and last ones are just dark current
wavelen= FLTARR(n-2)
overscan=FLTARR(n) ; to check if there are overscans
exposure=FLTARR(n) ; to check if there are bad exposures

imx    = FLTARR(nx,ny,n-2)
name   = 'qqq'
ix     = [2098-2047+INDGEN(2048),2101+INDGEN(2048)]
iy     = [2067-2047+INDGEN(2048),2132+INDGEN(2048)]
iy0    =  2068+INDGEN(64)

WINDOW,0,retain=2,xsize=800,ysize=600
!p.multi=0

; WE READ THE 2 DARK CURRENTS
;-----------------------------------------
;OPENR,1,'list'
OPENR,1,list
READF,1,name
PRINT,'1st DARK CURRENT= ',name
im     = READFITS(name,head,/silent)
dark1  = im[ix,*]                   ; to suppress vertical central line
overscan[0]= REBIN(dark1[*,iy0]+0.0,1)
dark1  = dark1[*,iy]+ 0.0          ; to suppress horizontal central line
year   = LONG(STRMID(head[8],19,4))
month  = LONG(STRMID(head[8],24,2))
day    = LONG(STRMID(head[8],27,2))
heure  = LONG(STRMID(head[9],21,2)) ; !!!!!!!!!!!!! I ASSUME THE HOUR IS UT !!!!!!!!!!!
minute = LONG(STRMID(head[9],24,2))
seconde= LONG(STRMID(head[9],27,2))
timedark1   = heure*3600.d0 + minute*60.d0 + seconde
corner = REBIN(dark1[0:63,0:63],1,1)
exposure[0]= corner
dark1  = dark1-corner[0]
dark1  = REBIN(dark1,nx,ny)
PRINT,'EXPOSURE TIME',head[68]

FOR i=1,n-2 DO READF,1,name
READF,1,name

PRINT,'2nd DARK CURRENT= ',name
im     = READFITS(name,head,/silent)
dark2  = im[ix,*]       ; to suppress vertical central line
overscan[n-1]= REBIN(dark2[*,iy0]+0.0,1)
dark2  = dark2[*,iy] + 0.0
year   = LONG(STRMID(head[8],19,4))
month  = LONG(STRMID(head[8],24,2))
day    = LONG(STRMID(head[8],27,2))
heure  = LONG(STRMID(head[9],21,2)) ; !!!!!!!!!!!!! I ASSUME THE HOUR IS UT !!!!!!!!!!!
minute = LONG(STRMID(head[9],24,2))
seconde= LONG(STRMID(head[9],27,2))
timedark2   = heure*3600.d0 + minute*60.d0 + seconde
corner = REBIN(dark2[0:63,0:63],1,1)
exposure[n-1]= corner
dark2  = dark2-corner[0]
dark2 = REBIN(dark2,nx,ny)
PRINT,'EXPOSURE TIME',head[68]
CLOSE,1

; WE READ THE FILTERGRAMS
;-----------------------------------------
;OPENR,1,'list'
OPENR,1,list
READF,1,name
FOR i=0,n-3 DO BEGIN
    READF,1,name
    PRINT,' '
    PRINT,'FILE READ= ',name
    im = READFITS(name,head,/silent)
    year   = LONG(STRMID(head[8],19,4))
    month  = LONG(STRMID(head[8],24,2))
    day    = LONG(STRMID(head[8],27,2))
    heure  = LONG(STRMID(head[9],21,2)) ; !!!!!!!!!!!!! I ASSUME THE HOUR IS UT !!!!!!!!!!!
    minute = LONG(STRMID(head[9],24,2))
    seconde= LONG(STRMID(head[9],27,2))
    time   = heure*3600.d0 + minute*60.d0 + seconde



  ; WE REMOVE THE DARK CURRENT
    q      = im[ix,*]                 ; to suppress vertical central line
    overscan[i+1]=REBIN(q[*,iy0]+0.0,1)
    q      = q[*,iy]+0.0              ; to suppress horizontal central line
    corner = REBIN(q[0:63,0:63],1,1)
    exposure[i+1]= corner
    q      = q-corner[0]

    imx[*,*,i] = REBIN(q,nx,ny)       ; we rebin the filtergram
    TVIM,imx[*,*,i],/scale

    FOR j=0,ny-1 DO BEGIN
        FOR jj=0,nx-1 DO BEGIN
            imx[jj,j,i] = imx[jj,j,i] - INTERPOL([dark1[jj,j],dark2[jj,j]],[timedark1,timedark2],time)
        ENDFOR
    ENDFOR
    temp = REFORM(imx[*,*,i])
    a    = WHERE(temp LT 0.d0)
    IF (a[0] NE -1) THEN temp[a]= 0.d0
    imx[*,*,i]=temp


  ; HELIOSTAT ROTATION ANGLE
;    angle  = ROTATION(year,month,day,time) ; RETURNS ANGLE IN DEGREES, POSITIVE ANGLE FOR COUNTER-CLOCKWISE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    PRINT,'OBSERVATION DATE',day,month,year
    PRINT,'OBSERVATION TIME',heure,minute,seconde
    PRINT,'EXPOSURE TIME',head[68]
;    PRINT,'ROTATION ANGLE',angle


  ; P-ANGLE
;    command = '/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day),1)+'.'+STRTRIM(STRING(heure),1)+'.'+STRTRIM(STRING(minute),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -17 | head -1 | cut -b19-32'
;    SPAWN,command,res
;    pangle = DOUBLE(res)
;    PRINT,'P-ANGLE',pangle ; TO E FROM N: MEANING IN DEGREES COUNTER-CLOCKWISE

;    angle = angle + pangle   ; !!!!!! NEED TO ADD CONSTANT DEPENDING ON THE HELIOSTAT SETTING !!!!!!
;    q[*,*] = ROT(q[*,*],angle[0]) ; ROTATES IN DEGREES CLOCKWISE: POSITIVE ANGLES ROTATE FROM LEFT TO RIGHT




  ; THE LINE-OF-SIGHT VELOCITY DUE TO THE EARTH MOVEMENT HAS TWO COMPONENTS
  ; THE EARTH SPIN + EARTH RADIAL VELOCITY TOWARD/AWAY TO/FROM THE SUN
  ; WE ASSUME THESE 2 COMPONENTS ARE CONSTANT DURING A DETUNE SEQUENCE
  ; (SHOULD CHANGE BY LESS THAN 1 m/s FOR AN EXPOSURE TIME OF 8s)

  ; Sun-Earth distance 5 minutes before the observation time
    heure1 = heure
    IF(minute GE 5) THEN BEGIN
        minute1 = minute - 5 
        day1    = day
    ENDIF ELSE BEGIN
        IF (minute EQ 4) THEN minute1 = 59
        IF (minute EQ 3) THEN minute1 = 58
        IF (minute EQ 2) THEN minute1 = 57
        IF (minute EQ 1) THEN minute1 = 56
        IF (minute EQ 0) THEN minute1 = 55
        heure1  = heure - 1
        IF(heure1 LT 0) THEN BEGIN
            heure1 = 23
            day1=day-1
        ENDIF ELSE day1 = day
    ENDELSE
    command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day1),1)+'.'+STRTRIM(STRING(heure1),1)+'.'+STRTRIM(STRING(minute1),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -13 | head -1 | cut -b13-25'
    SPAWN,command,res1
    res1 = DOUBLE(res1)
    
  ; Sun-Earth distance 5 minutes after the observation time
    heure1 = heure
    IF(minute LE 54) THEN BEGIN
        minute1 = minute + 5 
        day1    = day
    ENDIF ELSE BEGIN
        IF (minute EQ 55) THEN minute1 = 0
        IF (minute EQ 56) THEN minute1 = 1
        IF (minute EQ 57) THEN minute1 = 2
        IF (minute EQ 58) THEN minute1 = 3
        IF (minute EQ 59) THEN minute1 = 4
        heure1  = heure + 1
        IF(heure1 GE 24) THEN BEGIN
            heure1 = 0
            day1   = day +1
        ENDIF ELSE day1 = day
    ENDELSE
    command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day1),1)+'.'+STRTRIM(STRING(heure1),1)+'.'+STRTRIM(STRING(minute1),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -13 | head -1 | cut -b13-25'
    SPAWN,command,res2
    res2 = DOUBLE(res2)
    
  ; radial velocity in m/s (v>0 FOR MOVEMENTS AWAY FROM SUN)
    radial = ( res2 - res1 )*149597870691.d0/10.d0/60.d0 ; Sun-Earth distance in AU taken 6-minute apart

  ; line-of-sight velocity to the Earth spin in m/s 
  ; NB: THE VELOCITY RETURNED BY suninfo COMES FROM THE FUNCTION delta_v    
  ; WHICH RETURNS A SPEED dv SUCH AS V_ACTUAL = V_OBS + dv
  ; WE HAVE: V_ACTUAL + V_OBSERVER = V_OBSERVED
  ; WE WANT V_OBSERVED TO BE ABLE TO PREDICT THE SOLAR LINE OBSERVED WAVELENGTH

    command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day),1)+'.'+STRTRIM(STRING(heure),1)+'.'+STRTRIM(STRING(minute),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -1 | cut -b14-28'
    SPAWN,command,res3
    spin   = DOUBLE(res3)

    los    = radial[0] - spin[0]      ; line-of-sight velocity in m/s
    PRINT,'LINE-OF-SIGHT VELOCITY',los,radial[0],spin[0]




  ; REBIN


    wavelen[i] = los
    
ENDFOR
CLOSE,1

lam0 = MEAN(wavelen)
PRINT,'FINAL LINE-OF-SIGHT VELOCITY',lam0,SIGMA(wavelen)

SAVE,imx,dark1,dark2,lam0,overscan,exposure,FILE='SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

RETURN,imx

END
