; TO GET 9 DIFFERENT DETERMINATIONS OF THE NB MICHELSON
; PHASE AND CONTRAST MAPS

; PROGRAM TO PERFORM THE WAVELENGTH AND SPATIAL DEPENDENCE CALIBRATION
; TEST, WITH THE DYE LASER AS A SOURCE, AND FOR THE DE-TUNE SEQUENCE
; OBJECTIVE: TO LEARN HOW TO CO-TUNE THE ELEMENTS
; AND THE PHASES AND CONTRASTS OF THE TUNABLE ELEMENTS
;
; THIS PROGRAMS MENTIONS THE PRESENCE OF BAD EXPOSURES AND OVERSCANS
;
; ver 1.3 June 12, 2006
;
;----------------------------------------------------------------------

;----------------------------------------------------------------------
;
; MAIN PROGRAM
;
;----------------------------------------------------------------------

PRO HMI_laser1temp

time0       = SYSTIME(1)
dpi         = 2.d0*!dpi
lamref      = 6173.3433d0 ; solar Fe I 6173 line central wavelength in air
; when phases are at 0, we are centered on this wavelength

; VARIOUS PARAMETERS AND VARIABLES DEFINITION
;--------------------------------------------

nx          = 128              ; number of rows
ny          = 128              ; number of columns

Bg0         = FLTARR(nx,ny,9)  ; contrasts of NB MICHELSON, 9 TIMES
Phig0       = FLTARR(nx,ny,9)  ; phases of NB MICHELSON, 9 TIMES
nseq        = 27               ; number of positions in the de-tune sequence
Inten       = FLTARR(nx,ny,nseq);measured output intensities (measured on a HMI CCD)
anglim      = 1024.;980.;960.

xcenter     = nx/2;70;63;70 ; depends on the image !
ycenter     = ny/2;58;64;60

; DATA PROVIDED BY LOCKHEED-MARTIN (R. SHINE)
; e-mail 10/25/2005
; e-mail 11/18/2005
;--------------------------------------------

FSR         = 0.172457         ; FSR for the narrow-band Michelson



; INFORMATION TO BE ADDED BY THE USER
;----------------------------------------

;list         = 'listLaser060131_014331' ; in Calmode
;list         = 'listLaser060210_230522' ; in Calmode
;list         = 'listLaser060210_225517' ; in Obsmode
;list         = 'listLaser060210_221457' ; in Calmode
;list         = 'listLaser060303_181526' ; js_pcu_detune180 in Obsmode
;list         = 'listLaser060303_182759' ; js_pcu_detune180 in Obsmode
;list         = 'listLaser060303_184032' ; js_pcu_detune180 in Obsmode
;list         = 'listLaser060303_185253' ; js_pcu_detune180 in Obsmode
;list         = 'listLaser060303_190514' ; js_pcu_detune180 in Obsmode
;list         = 'listLaser060303_191737' ; js_pcu_detune180 in Obsmode
;list         = 'listLaser060620_232532' ; in Calmode
;list         = 'listLaser060620_000903' ; in Calmode LOT OF PROBLEMS
;list         = 'listLaser060620_234618' ; in Calmode
;list         = 'listLaser060622_195601' ; in Calmode
;list         = 'listLaser060622_182748' ; in Obsmode
;list         = 'listLaser060622_194336' ; in Obsmode
;list         = 'listLaser060622_201041' ; in Obsmode with diffuser
;list         = 'listLaser060622_202816' ; in Obsmode with diffuser
;list         = 'listLaser060622_221251' ; in Obsmode 
;list         = 'listLaser060622_214137' ; in Obsmode with diffuser
;list         = 'listLaser060622_223947' ; in Obsmode
;list         = 'listLaser060622_233055' ; in Obsmode
 list         = 'listSun060616_211450'
;list ='listSun060224_225806'
;wgood=INDGEN(29)

;wgood        = 0.0
;READIMAGES,list,images,headers,time,29,wgood,nbin=nx,iover=iover
;CLEAN,images,imx,headers
;SAVE,imx,time,wgood,FILE='SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
RESTORE,'SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'
;Inten = imx[*,*,1:27]

inten = imx
;goto,fuck

;orner        = REBIN(images(0:nx/64-1,0:nx/64-1,*),1,1,29)
;OR i         = 0,28 DO images[*,*,i] = images[*,*,i]-corner[i]
;
;arks         = FLTARR(nx,ny,2)
;arks[*,*,0]  = images[*,*,0]
;arks[*,*,1]  = images[*,*,28]
;nten         = images[*,*,1:27]
;imeimages    = FLTARR(nseq)
;imeimages    = time[1:27]
;imedarks     = [time[0],time[28]]
;
; WE REMOVE THE DARK FRAME
;RINT,'DARK FRAME REMOVAL'
;   FOR j=0,ny-1 DO BEGIN
;       FOR jj=0,nx-1 DO BEGIN
;           Inten[jj,j,*] = Inten[jj,j,*] - INTERPOL(REFORM(Darks[jj,j,*]),timedarks,timeimages[*])
;       ENDFOR
;   ENDFOR
;    = WHERE(Inten LT 0.d0)   
;F (a[0] NE -1) THEN Inten[a]= 0.d0


; LASER CENTRAL WAVELENGTH

CASE list OF
    'listLaser060131_014331': lam0          = 6173.538d0-lamref 
    'listLaser060210_221457': lam0          = 6173.585d0-lamref 
    'listLaser060210_225517': lam0          = 6173.463d0-lamref 
    'listLaser060210_230522': lam0          = 6173.463d0-lamref
    'listLaser060303_181526': lam0          = 6173.467d0-lamref 
    'listLaser060303_182759': lam0          = 6173.467d0-lamref
    'listLaser060303_184032': lam0          = 6173.467d0-lamref
    'listLaser060303_185253': lam0          = 6173.467d0-lamref
    'listLaser060303_190514': lam0          = 6173.467d0-lamref
    'listLaser060303_191737': lam0          = 6173.467d0-lamref
    'listLaser060620_232532': lam0          = 6173.227d0-lamref
    'listLaser060620_000903': lam0          = 6173.227d0-lamref
    'listLaser060620_234618': lam0          = 6173.227d0-lamref
    'listLaser060622_195601': lam0          = 6173.397d0-lamref
    'listLaser060622_194336': lam0          = 6173.397d0-lamref
    'listLaser060622_182748': lam0          = 6173.450d0-lamref
    'listLaser060622_201041': lam0          = 6173.397d0-lamref
    'listLaser060622_202816': lam0          = 6173.243d0-lamref
    'listLaser060622_221251': lam0          = 6173.512d0-lamref
    'listLaser060622_214137': lam0          = 6173.512d0-lamref
    'listLaser060622_223947': lam0          = 6173.640d0-lamref
    'listLaser060622_233055': lam0          = 6173.080d0-lamref
    'listSun060616_211450': lam0  = 0.005149
    'listSun060224_225806': lam0  = 0.013049383
ENDCASE

distance = SHIFT(DIST(nx,ny),xcenter,ycenter)*0.5d0*4096.d0/nx ; distance in arcseconds from the image center


; BEGINNING OF THE ITERATIONS
;----------------------------------------

mat      =  DBLARR(2,2)
vec      =  DBLARR(1,2)
mat[0,0] =  COS(dpi/FSR*lam0)
mat[1,0] = -SIN(dpi/FSR*lam0)
mat[0,1] =  COS(dpi/FSR*lam0 + dpi/3.d0)
mat[1,1] = -SIN(dpi/FSR*lam0 + dpi/3.d0)
mat[*,*] =  LA_INVERT(mat,/DOUBLE)


; BEGINNING OF COMPUTATIONS
In = FLTARR(nseq) 
FOR jjj=0,ny-1 DO BEGIN
    FOR iii=0,nx-1 DO BEGIN
  ; GUESS OF THE PHASES AND CONTRASTS OF TUNABLE ELEMENTS 
  ; ASSUMING THE LASER LINE PROFILE IS A DELTA FUNCTION
  ; we use the fact that cos(a)+cos(a+2pi/3)+cos(a+4pi/3) = 0
  ; I0 = [ cos(2\pi l0/FSR0)         -sin(2\pi l0/FSR0)         ] * B cos(\phi) * T/2 + T/2
  ; I1   [ cos(2\pi l0/FSR0+2\pi/3)  -sin(2\pi l0/FSR0+2\pi/3)  ]   B sin(\phi) * T/2 + T/2
  ; with 3 T/2 = I0+I1+I2
  ; USING THE MEASURED INTENSITIES
  ; THE FOLLOWING PROCEDURE IS WRITTEN FOR [1.5,1.5,-1.5]
  ; (WHEN ONE OF THE COEFFICIENT IS -1.5, THEN I1=I2 AND I2=I1) 
  ;----------------------------------------------------------

        IF(distance[iii,jjj] LE anglim) THEN BEGIN

            FOR kkk=0,8 DO BEGIN


                I00= Inten[iii,jjj,0+kkk*3]
                I1 = Inten[iii,jjj,1+kkk*3]
                I2 = Inten[iii,jjj,2+kkk*3]

                tot      = I00+I1+I2
                vec[0,0] = I00*3.d0/tot-1.d0
                vec[0,1] = I1 *3.d0/tot-1.d0
                vec      = mat##vec
                Bg0[iii,jjj,kkk]   = SQRT(vec[0,0]^2.d0+vec[0,1]^2.d0)

               ;WE CLEAN THE SOLUTION
                IF(Bg0[iii,jjj,kkk] GT 1.6d0) THEN Bg0[iii,jjj,kkk] = 1.6d0
                Phig0[iii,jjj,kkk] = ATAN(vec[0,1],vec[0,0])
                
            ENDFOR
            
        ENDIF ELSE BEGIN
            
            Bg0[iii,jjj,*]  = 0.d0
            Phig0[iii,jjj,*]= 0.d0
            
        ENDELSE
           
    ENDFOR
ENDFOR


; WE APPLY A MASK TO KEEP ONLY THE PHASES AND CONTRASTS ON THE ELEMENT
; APERTURE
a        = WHERE(distance GT anglim,COMPLEMENT=ba)

FOR i=0,8 DO BEGIN
    temp=REFORM(Phig0[*,*,i])
    temp[a]=-10000.d0
    phig0[*,*,i]=temp[*,*]
    temp=REFORM(Bg0[*,*,i])
    temp[a]=-10000.d0
    Bg0[*,*,i]=temp[*,*]
   ; TO GET RID OF THE STRIPES (LEFT AND UPPER EDGES)
   ;Phig0[0:5,*,i]=-10000.d0
   ;Phig0[*,122:127,i]=-10000.d0
   ;Bg0[0:5,*,i]=-10000.d0
   ;Bg0[*,122:127,i]=-10000.d0
ENDFOR

a = WHERE(FINITE(Phig0) EQ 0)
IF(a[0] NE -1) THEN Phig0[a]= -10000.0
a = WHERE(FINITE(Bg0) EQ 0)
IF(a[0] NE -1) THEN Bg0[a]  = -10000.0

a        = WHERE(REFORM(Bg0[*,*,0]) EQ -10000.d0,COMPLEMENT=ba)

; ADD 180 DEGREES TO PHASES < -180
FOR i=0,8 DO BEGIN
temp = REFORM(Phig0[*,*,i])
    a=WHERE(temp NE -10000.d0,COMPLEMENT=b)
    IF(MEAN(temp[a]) LT 0.0) THEN BEGIN
        b = WHERE(temp[a] GT 120.*!dpi/180.0)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]-dpi
        Phig0[*,*,i]=temp
    ENDIF ELSE BEGIN
        b = WHERE(temp[a] LT -120.*!dpi/180.0)
        IF(b[0] NE -1) THEN temp[a[b]]=temp[a[b]]+dpi
        Phig0[*,*,i]=temp
    ENDELSE
ENDFOR


; WE COMPUTE THE AVERAGE VALUE
moy=DBLARR(9)
mini=DBLARR(9)
maxi=DBLARR(9)
mini2=DBLARR(9)
maxi2=DBLARR(9)
FOR i=0,8 DO BEGIN
    temp=REFORM(Phig0[*,*,i])
    a=WHERE(temp NE -10000.d0,COMPLEMENT=b)
    moy[i]=MEAN(temp[a])
    mini[i]=MIN(temp[a])*180./!dpi
    maxi[i]=MAX(temp[a])*180./!dpi
    temp[b]=moy[i]
    Phig0[*,*,i]=temp
    temp=REFORM(Bg0[*,*,i])
    a=WHERE(temp NE -10000.d0,COMPLEMENT=b)
    temp[b]=MEAN(temp[a])
    mini2[i]=MIN(temp[a])
    maxi2[i]=MAX(temp[a])
    Bg0[*,*,i]=temp
ENDFOR

; WE PLOT THE RESULT
SET_PLOT,'ps'
!p.multi=[0,2,3]
device,file='yo.ps',bits=24,xoffset=0,yoffset=0,xsize=20,ysize=27,/color
loadct,4

FOR i=0,8 DO tvim,phig0[*,*,i]*180.d0/!dpi,/scale,tit=STRING(3*i)+STRING(3*i+1)+STRING(3*i+2),xtit='!17pixels',ytit='!17pixels',stit='!17Relative Phase (in degrees)',range=[mini[i],maxi[i]],pcharsize=1.5

FOR i=0,8 DO tvim,Bg0[*,*,i],/scale,tit=STRING(3*i)+STRING(3*i+1)+STRING(3*i+2),xtit='!17pixels',ytit='!17pixels',stit='!17Contrast',range=[mini2[i],maxi2[i]],pcharsize=1.5
temp = Bg0[*,*,0]

FOR i=0,8 DO TVIM,TOTAL(Inten[*,*,i*3:i*3+2],3),/scale,tit='!17Laser Intensity * n.t. transmittance',xtit='!17pixels',ytit='!17pixels'

DEVICE,/close
PRINT,'AVERAGES'

FOR i=0,8 DO PRINT,moy[i]*180.d0/!dpi


SET_PLOT,'x'
!p.multi=0
READ,pause

END
