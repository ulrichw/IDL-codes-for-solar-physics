; NB: TO AUTOMATICALLY FIND THE CENTER OF THE FIELD STOP ON THE
; OBSMODE IMAGES, USE:
; d[2047,2048]=0.0 & grid=fltarr(4096,4096) & for i=0,4095 do grid[*,i]=FINDGEN(4096) & a=where(d GT 600,complement=b) &  print,(MIN(grid[a])+MAX(grid[a]))/2.0
; and:
; d[2047,2048]=0.0 & grid=fltarr(4096,4096) & for i=0,4095 do grid[i,*]=FINDGEN(4096) & a=where(d GT 600,complement=b) &  print,(MIN(grid[a])+MAX(grid[a]))/2.0


;------------------------------------------------------------------------
; THIS PROGRAM IS FOR THE ANGULAR DEPENDENCE TEST IN CALMODE WITH 9
; FIELD STOPS POSITIONS, AND IN OBSMODE WITH 9 APERTURE STOPS
;------------------------------------------------------------------------

PRO HMI_angulardependence

nx       = 256
anglim   = 930.
dpi      = 2.d0*!dpi
distance = SHIFT(DIST(nx,ny),nx/2,nx/2)*0.5d0*4096.d0/nx

nim      = 9 ; number of angles
list     = STRARR(nim)
center   = FLTARR(nim,2)
OBSMODE  = 1.0

GOTO,next  ; 6173.340 A
list[0]='listLaser070225_85408' ; central field stop
list[1]='listLaser070225_85439'
list[2]='listLaser070225_85470'
list[3]='listLaser070225_85532'
list[4]='listLaser070225_85500'
list[5]='listLaser070225_85624'
list[6]='listLaser070225_85593'
list[7]='listLaser070225_85686'
list[8]='listLaser070225_85655'
center[0,*] = [2088,2173]; using algorithm
center[1,*] = [2662,2144];
center[2,*] = [3490,1960];
center[3,*] = [1470,2053];
center[4,*] = [586,2148];
center[5,*] = [2363,1392];
center[6,*] = [1890,374];
center[7,*] = [2189,2770];
center[8,*] = [2210,3521];

next:  ;6172.331 A
GOTO,next2
list[0]='listLaser070225_86978'
list[1]='listLaser070225_87009'
list[2]='listLaser070225_87040'
list[3]='listLaser070225_87071'
list[4]='listLaser070225_87102'
;list[5]='listLaser070225_87133' ; first darks missing, and then shitty data...
list[5]='listLaser070225_87164'
list[6]='listLaser070225_87164'
list[7]='listLaser070225_87195'
list[8]='listLaser070225_87226'
center[0,*] = [1845,1839]; using 4 points
center[1,*] = [1136,1693];
center[2,*] = [427,1898];
center[3,*] = [2756,2095];
center[4,*] = [3637,2191];
center[5,*] = [1938,3693];
center[6,*] = [1938,3693];
center[7,*] = [1991,1091];
center[8,*] = [1921,408];
next2: ; 6174.328 A
;GOTO,endd
list[0]='listLaser070225_87342'
list[3]='listLaser070225_87373'
list[4]='listLaser070225_87404'
list[1]='listLaser070225_87435'
list[2]='listLaser070225_87466'
list[5]='listLaser070225_87497'
list[6]='listLaser070225_87528'
list[7]='listLaser070225_87559'
list[8]='listLaser070225_87590'
center[0,*] = [2121,2056]; using automatic algorithm
center[3,*] = [1112,2018];
center[4,*] = [423,1948];
center[1,*] = [3043,2142];
center[2,*] = [3609,2090];
center[5,*] = [2038,1216];
center[6,*] = [1956,436];
center[7,*] = [2016,3034];
center[8,*] = [2117,3625];
endd:  ; 6173.343 OBSMODE !!!!
GOTO,jesus
OBSMODE=-1.
list[0]='listLaser070225_85780' ; 85808 missing + laser wavelength drifted
list[1]='listLaser070225_85882'
list[2]='listLaser070225_86018'
list[3]='listLaser070225_85814' 
list[4]='listLaser070225_85950'
;list[5]='listLaser070225_85848' ; crappy sequence
list[7]='listLaser070225_85984'
list[8]='listLaser070225_85984'
list[5]='listLaser070225_85916' 
list[6]='listLaser070225_86052'
center[0,*] = [2112,1988]
center[1,*] = [2925,1988]
center[2,*] = [3655,1989]
center[3,*] = [1300,1989]
center[4,*] = [505,1192]
center[7,*] = [2116,3589]
center[8,*] = [2116,3589]
center[5,*] = [2110,1174]
center[6,*] = [2110,442]
jesus:

center      = center*nx/4096.0

residuNB    = FLTARR(nim,nx,nx) ; phase difference between central aperture/field stop and other angles
residuWB    = FLTARR(nim,nx,nx)
residuE1    = FLTARR(nim,nx,nx)


; CALCULATE SHIFTS
NBimage = 5.81;6.1 ; image radius at full field (1026'') in mm
WBimage = 9.83;10.4
E1image = 11.94;12.6
NBcone  = 8.69 ; cone radius in the middle of the Michelson
WBcone  = 4.61 
E1cone  = 1.96

SHIFTSNB  = FLTARR(nim,2)
SHIFTSWB  = FLTARR(nim,2)
SHIFTSE1  = FLTARR(nim,2)

FOR i=1,nim-1 DO BEGIN
    SHIFTSNB[i,0] =  OBSMODE*(center[i,0]-center[0,0])*NBcone/NBimage*0.863027 ; not sure about the sign
    SHIFTSNB[i,1] =  OBSMODE*(center[i,1]-center[0,1])*NBcone/NBimage*0.863207
    SHIFTSWB[i,0] = -OBSMODE*(center[i,0]-center[0,0])*WBcone/WBimage*0.863207 ; empirical correction
    SHIFTSWB[i,1] = -OBSMODE*(center[i,1]-center[0,1])*WBcone/WBimage*0.863207 ; empirical correction
    SHIFTSE1[i,0] = -OBSMODE*(center[i,0]-center[0,0])*E1cone/E1image*0.863207
    SHIFTSE1[i,1] = -OBSMODE*(center[i,1]-center[0,1])*E1cone/E1image*0.863207 
ENDFOR


; WE RESTORE THE PHASE MAPS OF THE CENTRAL APERTURE/FIELD STOP
;-----------------------------------------------------------------------------------

RESTORE,'RESULTS/RESULTS_'+STRTRIM(list[0],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
Phig0central = Phig0
a  = WHERE(distance GT anglim)

FOR i=0,2 DO BEGIN
    temp=REFORM(Phig0central[*,*,i])
    temp[a]=-10000.d0
    phig0central[*,*,i]=temp[*,*]
ENDFOR
; ADD 180 DEGREES TO PHASES < -180
FOR i=0,2 DO BEGIN
temp = REFORM(Phig0central[*,*,i])
    a=WHERE(temp NE -10000.d0,COMPLEMENT=b)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 100.*!dpi/180.0)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-dpi
        Phig0central[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -100.*!dpi/180.0)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+dpi
        Phig0central[*,*,i]=temp
    ENDELSE
ENDFOR

FOR i=0,2 DO BEGIN
    temp = REFORM(Phig0central[*,*,i])
    a=WHERE(temp NE -10000.d0)
    PRINT,'MEAN=',MEAN(temp[a])*180./!pi
ENDFOR


SET_PLOT,'ps'
DEVICE,FILE='yo2.ps',xoffset=0,yoffset=0,xsize=26,ysize=20,/color,bits=24
!P.MULTI=[0,3,2]
LOADCT,4


; WE RESTORE THE PHASE MAPS OF ALL OTHER APERTURE/FIELD STOPS
;------------------------------------------------------------------------

FOR ii=1,nim-1 DO BEGIN

    RESTORE,'RESULTS/RESULTS_'+STRTRIM(list[ii],1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
    ;phig0=bg0

    IF SHIFTSNB[ii,0] LE 0.0 THEN Phig0[0:-SHIFTSNB[ii,0],*,0] = -10000.d0 ELSE Phig0[256-SHIFTSNB[ii,0]:255,*,0] = -10000.d0
    IF SHIFTSWB[ii,0] LE 0.0 THEN Phig0[0:-SHIFTSWB[ii,0],*,1] = -10000.d0 ELSE Phig0[256-SHIFTSWB[ii,0]:255,*,1] = -10000.d0
    IF SHIFTSE1[ii,0] LE 0.0 THEN Phig0[0:-SHIFTSE1[ii,0],*,2] = -10000.d0 ELSE Phig0[256-SHIFTSE1[ii,0]:255,*,2] = -10000.d0

    IF SHIFTSNB[ii,1] LE 0.0 THEN Phig0[*,0:-SHIFTSNB[ii,1],0] = -10000.d0 ELSE Phig0[*,256-SHIFTSNB[ii,1]:255,0] = -10000.d0
    IF SHIFTSWB[ii,1] LE 0.0 THEN Phig0[*,0:-SHIFTSWB[ii,1],1] = -10000.d0 ELSE Phig0[*,256-SHIFTSWB[ii,1]:255,1] = -10000.d0
    IF SHIFTSE1[ii,1] LE 0.0 THEN Phig0[*,0:-SHIFTSE1[ii,1],2] = -10000.d0 ELSE Phig0[*,256-SHIFTSE1[ii,1]:255,2] = -10000.d0

    Phig0[*,*,0] = SHIFT(Phig0[*,*,0],SHIFTSNB[ii,0],SHIFTSNB[ii,1])
    Phig0[*,*,1] = SHIFT(Phig0[*,*,1],SHIFTSWB[ii,0],SHIFTSWB[ii,1])
    Phig0[*,*,2] = SHIFT(Phig0[*,*,2],SHIFTSE1[ii,0],SHIFTSE1[ii,1])

    a            = WHERE(SHIFT(distance,SHIFTSNB[ii,0],SHIFTSNB[ii,1]) GT anglim,COMPLEMENT=ba)
    temp=REFORM(Phig0[*,*,0])
    temp[a]=-10000.d0
    phig0[*,*,0] = temp[*,*]
    
    a            = WHERE(SHIFT(distance,SHIFTSWB[ii,0],SHIFTSWB[ii,1]) GT anglim,COMPLEMENT=ba)
    temp=REFORM(Phig0[*,*,1])
    temp[a]=-10000.d0
    phig0[*,*,1] = temp[*,*]
    
    a            = WHERE(SHIFT(distance,SHIFTSE1[ii,0],SHIFTSE1[ii,1]) GT anglim,COMPLEMENT=ba)
    temp=REFORM(Phig0[*,*,2])
    temp[a]=-10000.d0
    phig0[*,*,2] = temp[*,*]

  ; ADD OR SUBTRACT 180 DEGREES TO PHASES < -180
  ; E1 HAS NEGATIVE PHASE CLOSE TO -180 DEGREES
    temp = REFORM(Phig0[*,*,2])
    a=WHERE(temp NE -10000.d0)
    b = WHERE(temp[a] GT 0.*!dpi/180.0)
    IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-dpi
    Phig0[*,*,2]=temp
  ; NB HAS NEGATIVE PHASES CLOSE TO 180 DEGREES
    temp = REFORM(Phig0[*,*,0])
    a=WHERE(temp NE -10000.d0)
    ;b = WHERE(temp[a] GT 150.*!dpi/180.0)
    ;IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-dpi
    ;Phig0[*,*,0]=temp
    ;b = WHERE(temp[a] LT -150.*!dpi/180.0)
    ;IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+dpi 
    ;Phig0[*,*,0]=temp
    IF (OBSMODE EQ 1) THEN BEGIN
        b = WHERE(temp[a] GT 0.*!dpi/180.0)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-dpi 
        phig0[*,*,0]=temp
    ENDIF
    IF(OBSMODE EQ -1) THEN BEGIN 
       b = WHERE(temp[a] LT 0.*!dpi/180.0)
       IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+dpi 
       Phig0[*,*,0]=temp
   ENDIF
   
    a = WHERE(Phig0[*,*,0] NE -10000.d0 AND Phig0central[*,*,0] NE -10000.d0,COMPLEMENT=b)
    temp=REFORM(Phig0[*,*,0]-Phig0central[*,*,0])
    temp[b]=-10000.d0
    residuNB[ii,*,*]=temp

    a = WHERE(Phig0[*,*,1] NE -10000.d0 AND Phig0central[*,*,1] NE -10000.d0,COMPLEMENT=b)
    temp=REFORM(Phig0[*,*,1]-Phig0central[*,*,1])
    temp[b]=-10000.d0
    residuWB[ii,*,*]=temp

    a = WHERE(Phig0[*,*,2] NE -10000.d0 AND Phig0central[*,*,2] NE -10000.d0,COMPLEMENT=b)
    temp=REFORM(Phig0[*,*,2]-Phig0central[*,*,2])
    temp[b]=-10000.0
    residuE1[ii,*,*]=temp

    TVIM,Phig0central[*,*,0],/scale
    TVIM,Phig0central[*,*,1],/scale
    TVIM,Phig0central[*,*,2],/scale
    TVIM,residuNB[ii,*,*],/scale
    TVIM,residuWB[ii,*,*],/scale
    TVIM,residuE1[ii,*,*],/scale

ENDFOR
DEVICE,/close


; WE DISPLAY THE RESULTS
;---------------------------------------------------------------------------------

SET_PLOT,'ps'
DEVICE,FILE='yo.ps',xoffset=0,yoffset=0,xsize=21,ysize=26,/color,bits=24
!P.MULTI=[0,2,2]
LOADCT,4

; RESULTS NB
NB   = [0.,0.,0.]
NB1=NB
NB2=NB
NB3=NB
angle= [0.,0.,0.] ; incidence angles
angle1=angle
angle2=angle
angle3=angle
a    = WHERE(phig0central[*,*,0] NE -10000.d0 AND RESIDUNB[1,*,*] NE -10000.d0 AND RESIDUNB[2,*,*] NE -10000.d0)
disc = FLTARR(nx,nx) & disc[a]=1.0
temp = phig0central[*,*,0]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
temp = RESIDUNB[1,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
PRINT,'NB ---RIGHT'
NB[1]=-MEAN(temp[a])*172./dpi
angle[1]=sqrt( (center[1,0]-center[0,0])^2.d0+(center[1,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'NB',NB[1],SIGMA(temp[a])*172./dpi
temp = RESIDUNB[2,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
NB[2]=-MEAN(temp[a])*172./dpi
angle[2]=sqrt( (center[2,0]-center[0,0])^2.d0+(center[2,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'NB',NB[2],SIGMA(temp[a])*172./dpi
PLOT,angle/3600.,NB,psym=4
res=POLY_FIT(angle,NB,2,yfit=y1)
OPLOT,angle/3600.,y1


a = WHERE(phig0central[*,*,0] NE -10000.d0 AND RESIDUNB[3,*,*] NE -10000.d0 AND RESIDUNB[4,*,*] NE -10000.d0)
disc = FLTARR(nx,nx)
disc[a]=1.0
temp = phig0central[*,*,0]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
temp = RESIDUNB[3,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
PRINT,'NB ---LEFT'
NB1[1]=-MEAN(temp[a])*172./dpi
angle1[1]=sqrt( (center[3,0]-center[0,0])^2.d0+(center[3,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'NB',NB1[1],SIGMA(temp[a])*172./dpi
temp = RESIDUNB[4,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
NB1[2]=-MEAN(temp[a])*172./dpi
angle1[2]=sqrt( (center[4,0]-center[0,0])^2.d0+(center[4,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'NB',NB1[2],SIGMA(temp[a])*172./dpi
PLOT,angle1/3600.,NB1,psym=4
res=POLY_FIT(angle1,NB1,2,yfit=y2)
OPLOT,angle1/3600.,y2


a = WHERE(phig0central[*,*,0] NE -10000.d0 AND RESIDUNB[7,*,*] NE -10000.d0 AND RESIDUNB[8,*,*] NE -10000.d0)
disc = FLTARR(nx,nx)
disc[a]=1.0
temp = phig0central[*,*,0]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
temp = RESIDUNB[7,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
PRINT,'NB ---UP'
NB2[1]=-MEAN(temp[a])*172./dpi
angle2[1]=sqrt( (center[7,0]-center[0,0])^2.d0+(center[7,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'NB',NB2[1],SIGMA(temp[a])*172./dpi
temp = RESIDUNB[8,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
NB2[2]=-MEAN(temp[a])*172./dpi
angle2[2]=sqrt( (center[8,0]-center[0,0])^2.d0+(center[8,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'NB',NB2[2],SIGMA(temp[a])*172./dpi
PLOT,angle2/3600.,NB2,psym=4
res=POLY_FIT(angle2,NB2,2,yfit=y3)
OPLOT,angle2/3600.,y3

a = WHERE(phig0central[*,*,0] NE -10000.d0 AND RESIDUNB[5,*,*] NE -10000.d0 AND RESIDUNB[6,*,*] NE -10000.d0)
disc = FLTARR(nx,nx)
disc[a]=1.0
temp = phig0central[*,*,0]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
temp = RESIDUNB[5,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
PRINT,'NB ---DOWN'
NB3[1]=-MEAN(temp[a])*172./dpi
angle3[1]=sqrt( (center[5,0]-center[0,0])^2.d0+(center[5,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'NB',NB3[1],SIGMA(temp[a])*172./dpi
temp = RESIDUNB[6,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
NB3[2]=-MEAN(temp[a])*172./dpi
angle3[2]=sqrt( (center[6,0]-center[0,0])^2.d0+(center[6,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'NB',NB3[2],SIGMA(temp[a])*172./dpi
PLOT,angle3/3600.,NB3,psym=4
res=POLY_FIT(angle3,NB3,2,yfit=y4)
OPLOT,angle3/3600.,y4


; RESULTS WB
angle4=[0.,0.,0.]
angle5=angle
angle6=angle
angle7=angle
WB=[0.,0.,0.]
WB1=WB
WB2=WB
WB3=WB
a = WHERE(phig0central[*,*,1] NE -10000.d0 AND RESIDUWB[1,*,*] NE -10000.d0 AND RESIDUWB[2,*,*] NE -10000.d0)
disc = FLTARR(nx,nx)
disc[a]=1.0
temp = phig0central[*,*,1]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
temp = RESIDUWB[1,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
PRINT,'WB ---RIGHT'
WB[1]=-MEAN(temp[a])*344./dpi
angle4[1]=sqrt( (center[1,0]-center[0,0])^2.d0+(center[1,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'WB',WB[1],SIGMA(temp[a])*344./dpi
temp = RESIDUWB[2,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
WB[2]=-MEAN(temp[a])*344./dpi
angle4[2]=sqrt( (center[2,0]-center[0,0])^2.d0+(center[2,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'WB',WB[2],SIGMA(temp[a])*344./dpi
PLOT,angle4/3600.,WB,psym=4
res=POLY_FIT(angle4,WB,2,yfit=y5)
OPLOT,angle4/3600.,y5

a = WHERE(phig0central[*,*,1] NE -10000.d0 AND RESIDUWB[3,*,*] NE -10000.d0 AND RESIDUWB[4,*,*] NE -10000.d0)
disc = FLTARR(nx,nx)
disc[a]=1.0
temp = phig0central[*,*,1]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
temp = RESIDUWB[3,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
PRINT,'WB ---LEFT'
WB1[1]=-MEAN(temp[a])*344./dpi
angle5[1]=sqrt( (center[3,0]-center[0,0])^2.d0+(center[3,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'WB',WB1[1],SIGMA(temp[a])*344./dpi
temp = RESIDUWB[4,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
WB1[2]=-MEAN(temp[a])*344./dpi
angle5[2]=sqrt( (center[4,0]-center[0,0])^2.d0+(center[4,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'WB',WB1[2],SIGMA(temp[a])*344./dpi
PLOT,angle5/3600.,WB1,psym=4
res=POLY_FIT(angle5,WB1,2,yfit=y6)
OPLOT,angle5/3600.,y6

a = WHERE(phig0central[*,*,1] NE -10000.d0 AND RESIDUWB[7,*,*] NE -10000.d0 AND RESIDUWB[8,*,*] NE -10000.d0)
disc = FLTARR(nx,nx)
disc[a]=1.0
temp = phig0central[*,*,1]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
temp = RESIDUWB[7,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
PRINT,'WB ---UP'
WB2[1]=-MEAN(temp[a])*344./dpi
angle6[1]=sqrt( (center[7,0]-center[0,0])^2.d0+(center[7,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'WB',WB2[1],SIGMA(temp[a])*344./dpi
temp = RESIDUWB[8,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
WB2[2]=-MEAN(temp[a])*344./dpi
angle6[2]=sqrt( (center[8,0]-center[0,0])^2.d0+(center[8,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'WB',WB2[2],SIGMA(temp[a])*344./dpi
PLOT,angle6/3600.,WB2,psym=4
res=POLY_FIT(angle6,WB2,2,yfit=y7)
OPLOT,angle6/3600.,y7

a = WHERE(phig0central[*,*,1] NE -10000.d0 AND RESIDUWB[5,*,*] NE -10000.d0 AND RESIDUWB[6,*,*] NE -10000.d0)
disc = FLTARR(nx,nx)
disc[a]=1.0
temp = phig0central[*,*,1]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
temp = RESIDUWB[5,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
PRINT,'WB ---DOWN'
WB3[1]=-MEAN(temp[a])*344./dpi
angle7[1]=sqrt( (center[5,0]-center[0,0])^2.d0+(center[5,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'WB',WB3[1],SIGMA(temp[a])*344./dpi
temp = RESIDUWB[6,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
WB3[2]=-MEAN(temp[a])*344./dpi
angle7[2]=sqrt( (center[6,0]-center[0,0])^2.d0+(center[6,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'WB',WB3[2],SIGMA(temp[a])*344./dpi
PLOT,angle7/3600.,WB3,psym=4
res=POLY_FIT(angle7,WB3,2,yfit=y8)
OPLOT,angle7/3600.,y8


; RESULTS E1
E1=[0.,0.,0.]
E11=E1
E12=E1
E13=E1
angle8=[0.,0.,0.]
angle9=angle8
angle10=angle8
angle11=angle8
a = WHERE(phig0central[*,*,2] NE -10000.d0 AND RESIDUE1[1,*,*] NE -10000.d0 AND RESIDUE1[2,*,*] NE -10000.d0)
disc = FLTARR(nx,nx)
disc[a]=1.0
temp = phig0central[*,*,2]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
temp = RESIDUE1[1,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
PRINT,'E1 ---RIGHT'
E1[1]=-MEAN(temp[a])*690./dpi
angle8[1]=sqrt( (center[1,0]-center[0,0])^2.d0+(center[1,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'E1',E1[1],SIGMA(temp[a])*690./dpi
temp = RESIDUE1[2,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
E1[2]=-MEAN(temp[a])*690./dpi
angle8[2]=sqrt( (center[2,0]-center[0,0])^2.d0+(center[2,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'E1',E1[2],SIGMA(temp[a])*690./dpi
PLOT,angle8/3600.,E1,psym=4
res=POLY_FIT(angle8,E1,2,yfit=y9)
OPLOT,angle8/3600.,y9

a = WHERE(phig0central[*,*,2] NE -10000.d0 AND RESIDUE1[3,*,*] NE -10000.d0 AND RESIDUE1[4,*,*] NE -10000.d0)
disc = FLTARR(nx,nx)
disc[a]=1.0
temp = phig0central[*,*,2]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
temp = RESIDUE1[3,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
PRINT,'E1 ---LEFT'
E11[1]=-MEAN(temp[a])*690./dpi
angle9[1]=sqrt( (center[3,0]-center[0,0])^2.d0+(center[3,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'E1',E11[1],SIGMA(temp[a])*690./dpi
temp = RESIDUE1[4,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
E11[2]=-MEAN(temp[a])*690./dpi
angle9[2]=sqrt( (center[4,0]-center[0,0])^2.d0+(center[4,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'E1',E11[2],SIGMA(temp[a])*690./dpi
PLOT,angle9/3600.,E11,psym=4
res=POLY_FIT(angle9,E11,2,yfit=y10)
OPLOT,angle9/3600.,y10


a = WHERE(phig0central[*,*,2] NE -10000.d0 AND RESIDUE1[7,*,*] NE -10000.d0 AND RESIDUE1[8,*,*] NE -10000.d0)
disc = FLTARR(nx,nx)
disc[a]=1.0
temp = phig0central[*,*,2]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
temp = RESIDUE1[7,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
PRINT,'E1 ---UP'
E12[1]=-MEAN(temp[a])*690./dpi
angle10[1]=sqrt( (center[7,0]-center[0,0])^2.d0+(center[7,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'E1',E12[1],SIGMA(temp[a])*690./dpi
temp = RESIDUE1[8,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
E12[2]=-MEAN(temp[a])*690./dpi
angle10[2]=sqrt( (center[8,0]-center[0,0])^2.d0+(center[8,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'E1',E12[2],SIGMA(temp[a])*690./dpi
PLOT,angle10/3600.,E12,psym=4
res=POLY_FIT(angle10,E12,2,yfit=y11)
OPLOT,angle10/3600.,y11

a = WHERE(phig0central[*,*,2] NE -10000.d0 AND RESIDUE1[5,*,*] NE -10000.d0 AND RESIDUE1[6,*,*] NE -10000.d0)
disc = FLTARR(nx,nx)
disc[a]=1.0
temp = phig0central[*,*,2]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
temp = RESIDUE1[5,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
PRINT,'E1 ---DOWN'
E13[1]=-MEAN(temp[a])*690./dpi
angle11[1]=sqrt( (center[5,0]-center[0,0])^2.d0+(center[5,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'E1',E13[1],SIGMA(temp[a])*690./dpi
temp = RESIDUE1[6,*,*]*disc
TVIM,temp*180./!pi,range=[MIN(temp[a]),MAX(temp[a])]*180./!pi,/scale,barwidth=0.5
E13[2]=-MEAN(temp[a])*690./dpi
angle11[2]=sqrt( (center[6,0]-center[0,0])^2.d0+(center[6,1]-center[0,1])^2.d0)*1024.*2./FLOAT(nx)
PRINT,'E1',E13[2],SIGMA(temp[a])*690./dpi
PLOT,angle11/3600.,E13,psym=4
res=POLY_FIT(angle11,E13,2,yfit=y12)
OPLOT,angle11/3600.,y12

DEVICE,/CLOSE

!P.MULTI=[0,1,3]
DEVICE,FILE='yo3.ps',xoffset=0,yoffset=0.5,xsize=20,ysize=26

PLOT,angle/3600.,NB,psym=4,yrange=[-15,7],xrange=[0,0.3],xst=1,yst=1,charsize=1.5,tit='!17Narrow-Band Michelson',xtit='Angle (degrees)',ytit='Wavelength shift (mA)'
OPLOT,angle/3600.,y1
OPLOT,angle1/3600.,NB1,psym=4,linestyle=2
OPLOT,angle2/3600.,NB2,psym=4,linestyle=3
OPLOT,angle3/3600.,NB3,psym=4,linestyle=4
OPLOT,angle1/3600.,y2,linestyle=2
IF (OBSMODE NE -1) THEN OPLOT,angle2/3600.,y3,linestyle=3
OPLOT,angle3/3600.,y4,linestyle=4

LEGEND,['right','left','up','down'],linestyle=[0,2,3,4],/bottom
;angle=[angle,angle1,angle2,angle3]
;NB=[NB,NB1,NB2,NB3]
;a=SORT(angle)
;angle=angle[a]
;NB=NB[a]
;res=POLY_FIT(angle,NB,2,yfit=y)
;OPLOT,angle/3600.,y

PLOT,angle4/3600.,WB,psym=4,yrange=[-15,7],xrange=[0,0.3],xst=1,yst=1,charsize=1.5,tit='!17Broad-Band Michelson',xtit='Angle (degrees)',ytit='Wavelength shift (mA)'
OPLOT,angle4/3600.,y5
OPLOT,angle5/3600.,WB1,psym=4,linestyle=2
OPLOT,angle6/3600.,WB2,psym=4,linestyle=3
OPLOT,angle7/3600.,WB3,psym=4,linestyle=4
OPLOT,angle5/3600.,y6,linestyle=2
IF (OBSMODE NE -1) THEN OPLOT,angle6/3600.,y7,linestyle=3
OPLOT,angle7/3600.,y8,linestyle=4
LEGEND,['right','left','up','down'],linestyle=[0,2,3,4],/bottom
;angle=[angle4,angle5,angle6,angle7]
;WB=[WB,WB1,WB2,WB3]
;a=SORT(angle)
;angle=angle[a]
;WB=WB[a]
;res=POLY_FIT(angle,WB,2,yfit=y)
;OPLOT,angle/3600.,y

PLOT,angle8/3600.,E1,psym=4,yrange=[-6,60],xrange=[0,0.3],xst=1,yst=1,charsize=1.5,tit='!17Lyot Element E1',xtit='Angle (degrees)',ytit='Wavelength shift (mA)'
OPLOT,angle8/3600.,y9
OPLOT,angle9/3600.,E11,psym=4,linestyle=2
OPLOT,angle10/3600.,E12,psym=4,linestyle=3
OPLOT,angle11/3600.,E13,psym=4,linestyle=4
OPLOT,angle9/3600.,y10,linestyle=2
IF (OBSMODE NE -1) THEN OPLOT,angle10/3600.,y11,linestyle=3
OPLOT,angle11/3600.,y12,linestyle=4
LEGEND,['right','left','up','down'],linestyle=[0,2,3,4],/top
;angle=[angle8,angle9,angle10,angle11]
;E1=[E1,E11,E12,E13]
;a=SORT(angle)
;angle=angle[a]
;E1=E1[a]
;res=POLY_FIT(angle,E1,2,yfit=y)
;OPLOT,angle/3600.,y

DEVICE,/close

SET_PLOT,'X'
read,pause


END
