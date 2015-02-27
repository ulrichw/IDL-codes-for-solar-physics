;----------------------------------------------------------------------
;
; SUBROUTINE TO 1) READ THE .FITS FILES CONTAINING THE FILTERGRAMS
; 2) READ THE DYE LASER WAVELENGTH  
; 3) REBINS THE FILTERGRAMS
; 4) NORMALIZE EACH IMAGE BY THE DYE LASER RELATIVE INTENSITY
; based on a routine provided by J. Schou
;
;----------------------------------------------------------------------


FUNCTION readfile_detune180,header,nx,base,list
; THE ROUTINE ROTATION.PRO MUST BE COMPILED PRIOR
; TO THIS CODE
; ROTATION.PRO RETURNS THE ROTATION ANGLE DUE TO THE HELIOSTAT
; suninfo PROVIDES THE p-angle

SET_PLOT,'x'
n      = 29 ; numbers of filtergrams: the 1st and last ones are just dark current
nkey   = 91
nhole  = 6  ; number of holes for the angular dependence test
header = STRARR((n-2)*nhole,nkey)
imx    = FLTARR(nx,nx,(n-2)*nhole)
calmode= FLTARR(nx,nx,nhole)
name   = 'qqq'
ix     = [2098-2047+INDGEN(2048),2101+INDGEN(2048)]
iy     = [2067-2047+INDGEN(2048),2132+INDGEN(2048)]

WINDOW,0,retain=2,xsize=800,ysize=600
!p.multi=0

OPENR,1,list
FOR hole=0,nhole-1 DO BEGIN

   ; WE READ THE CALMODE IMAGE
   ;-----------------------------------------
 
    READF,1,name
    im     = READFITS(base+name,head,/silent)
    im     = im[ix,*] + 0.d0   
    im     = im[*,iy]
    calmode[*,*,hole] = REBIN(im,nx,nx)

   ; WE READ THE FIRST DARK CURRENT
   ;-----------------------------------------
    READF,1,name
    PRINT,'1st DARK CURRENT= ',name
    im     = READFITS(base+name,head,/silent)
    dark1  = im[ix,*] + 0.d0    ; to suppress vertical central line
    dark1  = dark1[*,iy]        ; to suppress horizontal central line
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
    
    calmode[*,*,hole]=calmode[*,*,hole]-dark1
    TVIM,calmode[*,*,hole]

  ; WE READ THE FILTERGRAMS
  ;-----------------------------------------
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
        
        header[i+hole*(n-2),*] = head
        
        q      = im[ix,*] + 0.d0 ; to suppress vertical central line
        q      = q[*,iy]        ; to suppress horizontal central line
        
        imx[*,*,i+hole*(n-2)] = REBIN(q,nx,nx) ; we rebin the filtergrams

        PRINT,'OBSERVATION DATE',day,month,year
        PRINT,'OBSERVATION TIME',heure,minute,seconde
        PRINT,'LASER WAVELENGTH',head[73]
        PRINT,'LASER INTENSITY',head[74]
        PRINT,'EXPOSURE TIME',head[68]
        
    ENDFOR

   ; WE READ THE SECOND DARK CURRENT
   ;----------------------------------------------

    READF,1,name
    PRINT,'2nd DARK CURRENT= ',name
    im     = READFITS(base+name,head,/silent)
    dark2  = im[ix,*] + 0.d0    ; to suppress vertical central line
    dark2  = dark2[*,iy]
    year   = LONG(STRMID(head[8],19,4))
    month  = LONG(STRMID(head[8],24,2))
    day    = LONG(STRMID(head[8],27,2))
    heure  = LONG(STRMID(head[9],21,2)) ; !!!!!!!!!!!!! I ASSUME THE HOUR IS UT !!!!!!!!!!!
    minute = LONG(STRMID(head[9],24,2))
    seconde= LONG(STRMID(head[9],27,2))
    timedark2   = heure*3600.d0 + minute*60.d0 + seconde
    dark2 = REBIN(dark2,nx,nx)
    PRINT,'OBSERVATION DATE',day,month,year
    PRINT,'OBSERVATION TIME',heure,minute,seconde
    PRINT,'EXPOSURE TIME',head[68]

;   LINEAR INTERPOLATION AS SUGGESTED BY JESPER
    FOR i=0,n-3 DO BEGIN
        FOR j=0,nx-1 DO BEGIN
            FOR jj=0,nx-1 DO BEGIN
                imx[jj,j,i+hole*(n-2)] = imx[jj,j,i+hole*(n-2)] - INTERPOL([dark1[jj,j],dark2[jj,j]],[timedark1,timedark2],time)
            ENDFOR
        ENDFOR
        temp = REFORM(imx[*,*,i+hole*(n-2)])
        a    = WHERE(temp LT 0.d0) 
        IF (a[0] NE -1) THEN temp[a]= 0.d0
        imx[*,*,i+hole*(n-2)]=temp

        TVIM,imx[*,*,i+hole*(n-2)],/scale
    ENDFOR

ENDFOR
CLOSE,1


; WE GET RID OF THE BAD EXPOSURES
corner = REBIN(imx(0:nx/64-1,0:nx/64-1,*),1,1,(n-2)*nhole)
medc   = MEDIAN(corner)
wgood  = WHERE(corner LE medc+100)


; WE FREE MEMORY
dark1 = 0.d0
dark2 = 0.d0
q     = 0.d0
im    = 0.d0

SAVE,imx,header,calmode,wgood,FILE='SEQUENCE_DETUNE180_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

RETURN,imx

END
