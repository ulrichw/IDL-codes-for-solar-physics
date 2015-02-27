; ROUTINE THAT GIVES THE SUN-EARTH VELOCITY AT A SPECIFIC DATE
; USED TO KNOW THE LOCATION OF THE SOLAR IRON LINE ON A SPECTRUM
; TAKEN FROM EARTH
;
; NB: THIS CODE RUNS ON n00
;
; NUMERICAL DERIVATIVE COMPUTED BY A FIVE POINT ALGORITHM
;----------------------------------------------------------------------


; NB: the program sunearthvel has been updated
; the new version is in /auto/home0/wso/bin/_linux4/suninfo
; and contain, on top of the earth spin velocity, the total velocity sun-earth
; the sign error in the earth spin has also been corrected
; times are UT

; THE PARAMETER SHOULD ALL BE ENTERED AS INTEGER (NO FLOATS)

PRO sunearthvel,year,month,day,heure,minute,seconde

; TIME MUST BE UT

;time   = heure*3600.d0 + minute*60.d0 + seconde

PRINT,'OBSERVATION DATE',day,month,year
PRINT,'OBSERVATION TIME',heure,minute,seconde

  ; THE LINE-OF-SIGHT VELOCITY DUE TO THE EARTH MOVEMENT HAS TWO COMPONENTS
  ; THE EARTH SPIN + EARTH RADIAL VELOCITY TOWARD/AWAY TO/FROM THE SUN
  ; WE ASSUME THESE 2 COMPONENTS ARE CONSTANT DURING A DETUNE SEQUENCE
  ; (SHOULD CHANGE BY LESS THAN 1 m/s FOR AN EXPOSURE TIME OF 8s)

  ; Sun-Earth distance 6 minutes before the observation time
heure1 = heure
IF(minute GE 6) THEN BEGIN
    minute1 = minute - 6 
    day1    = day
ENDIF ELSE BEGIN
    IF (minute EQ 5) THEN minute1 = 59
    IF (minute EQ 4) THEN minute1 = 58
    IF (minute EQ 3) THEN minute1 = 57
    IF (minute EQ 2) THEN minute1 = 56
    IF (minute EQ 1) THEN minute1 = 55
    IF (minute EQ 0) THEN minute1 = 54
    heure1  = heure - 1
    IF(heure1 LT 0) THEN BEGIN
        heure1 = 23
        day1=day-1
    ENDIF ELSE day1 = day
ENDELSE
command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day1),1)+'.'+STRTRIM(STRING(heure1),1)+'.'+STRTRIM(STRING(minute1),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -13 | head -1 | cut -b13-25'
SPAWN,command,res1
res1 = DOUBLE(res1)

  ; Sun-Earth distance 3 minutes before the observation time
heure1 = heure
IF(minute GE 3) THEN BEGIN
    minute1 = minute - 3 
    day1    = day
ENDIF ELSE BEGIN
    IF (minute EQ 2) THEN minute1 = 59
    IF (minute EQ 1) THEN minute1 = 58
    IF (minute EQ 0) THEN minute1 = 57
    heure1  = heure - 1
    IF(heure1 LT 0) THEN BEGIN
        heure1 = 23
        day1=day-1
    ENDIF ELSE day1 = day
ENDELSE
command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day1),1)+'.'+STRTRIM(STRING(heure1),1)+'.'+STRTRIM(STRING(minute1),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -13 | head -1 | cut -b13-25'
SPAWN,command,res2
res2 = DOUBLE(res2)

  ; Sun-Earth distance at the observation time
heure1 = heure
minute1= minute
day1   = day
command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day1),1)+'.'+STRTRIM(STRING(heure1),1)+'.'+STRTRIM(STRING(minute1),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -13 | head -1 | cut -b13-25'
SPAWN,command,res3
res3 = DOUBLE(res3)
PRINT,'RADIAL DISTANCE (a.u.)=',res3

  ; Sun-Earth distance 3 minutes after the observation time
heure1 = heure
IF(minute LE 56) THEN BEGIN
    minute1 = minute + 3 
    day1    = day
ENDIF ELSE BEGIN
    IF (minute EQ 57) THEN minute1 = 0
    IF (minute EQ 58) THEN minute1 = 1
    IF (minute EQ 59) THEN minute1 = 2
    heure1  = heure + 1
    IF(heure1 GE 24) THEN BEGIN
        heure1 = 0
        day1   = day +1
    ENDIF ELSE day1 = day
ENDELSE
command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day1),1)+'.'+STRTRIM(STRING(heure1),1)+'.'+STRTRIM(STRING(minute1),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -13 | head -1 | cut -b13-25'
SPAWN,command,res4
res4 = DOUBLE(res4)

  ; Sun-Earth distance 6 minutes after the observation time
heure1 = heure
IF(minute LE 53) THEN BEGIN
    minute1 = minute + 6 
    day1    = day
ENDIF ELSE BEGIN
    IF (minute EQ 54) THEN minute1 = 0
    IF (minute EQ 55) THEN minute1 = 1
    IF (minute EQ 56) THEN minute1 = 2
    IF (minute EQ 57) THEN minute1 = 3
    IF (minute EQ 58) THEN minute1 = 4
    IF (minute EQ 59) THEN minute1 = 5
    heure1  = heure + 1
    IF(heure1 GE 24) THEN BEGIN
        heure1 = 0
        day1   = day +1
    ENDIF ELSE day1 = day
ENDELSE
command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day1),1)+'.'+STRTRIM(STRING(heure1),1)+'.'+STRTRIM(STRING(minute1),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -13 | head -1 | cut -b13-25'
SPAWN,command,res5
res5 = DOUBLE(res5)


 ; radial velocity in m/s (v>0 FOR MOVEMENTS AWAY FROM SUN)
;radial = ( res2 - res1 )*149597870691.d0/10.d0/60.d0 ; Sun-Earth distance in AU taken 6-minute apart
 radial = 1.d0/24.d0/(3.d0*60.d0)*(2.d0*res1-16.d0*res2+16.d0*res4-2.d0*res5)*149597870691.d0
   
  ; line-of-sight velocity to the Earth spin in m/s 
  ; NB: THE VELOCITY RETURNED BY suninfo COMES FROM THE FUNCTION delta_v    
  ; WHICH RETURNS A SPEED dv SUCH AS V_ACTUAL = V_OBS + dv
  ; WE HAVE: V_ACTUAL + V_OBSERVER = V_OBSERVED
  ; WE WANT V_OBSERVED TO BE ABLE TO PREDICT THE SOLAR LINE OBSERVED WAVELENGTH

command='/home/wso/bin/_linux/suninfo -puv '+STRTRIM(STRING(year),1)+'.'+STRTRIM(STRING(month),1)+'.'+STRTRIM(STRING(day),1)+'.'+STRTRIM(STRING(heure),1)+'.'+STRTRIM(STRING(minute),1)+'.'+STRTRIM(STRING(seconde),1)+'| tail -1 | cut -b14-28'
SPAWN,command,res3
spin   = DOUBLE(res3)

los    = radial[0] - spin[0]    ; line-of-sight velocity in m/s
PRINT,'LINE-OF-SIGHT VELOCITY',los,radial[0],-spin[0]
PRINT,'WAVELENGTH SHIFT',los*6173.3433/299792458. 



END
