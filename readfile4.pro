;----------------------------------------------------------------------
;
; SUBROUTINE TO 1) READ THE .FITS FILES CONTAINING THE REBINNED FILTERGRAMS
; 2) READ THE DYE LASER WAVELENGTH
; 3) NORMALIZE EACH IMAGE BY THE DYE LASER RELATIVE INTENSITY
; based on a routine provided by J. Schou
;
;----------------------------------------------------------------------


FUNCTION readfile4,header,nx,ny,base,list,nhole
; THE ROUTINE ROTATION.PRO MUST BE COMPILED PRIOR
; TO THIS CODE
; ROTATION.PRO RETURNS THE ROTATION ANGLE DUE TO THE HELIOSTAT
; suninfo PROVIDES THE p-angle

SET_PLOT,'x'
n      = 29 ; numbers of filtergrams: the 1st and last ones are just dark current
nkey   = 91
header = STRARR((n-2)*nhole,nkey)

imx    = FLTARR(nx,ny,nhole*(n-2))
calmode= FLTARR(nx,ny,nhole)
name   = 'qqq'
ix     = [2098-2047+INDGEN(2048),2101+INDGEN(2048)]
iy     = [2067-2047+INDGEN(2048),2132+INDGEN(2048)]

WINDOW,0,retain=2,xsize=1200,ysize=900
!p.multi=0


OPENR,1,list

FOR iii=0,nhole-1 DO BEGIN
     !p.multi=[0,5,6]

    ; WE READ THE CALMODE IMAGE
    ;-----------------------------------------

    READF,1,name
    PRINT,'CALMODE IMAGE= ',name
    im = READFITS(base+name,head,/silent)
    im = im[ix,*] + 0.d0        ; to suppress vertical central line
    im = im[*,iy]
    calmode[*,*,iii] = REBIN(im,nx,ny)
    
    year   = LONG(STRMID(head[8],19,4))
    month  = LONG(STRMID(head[8],24,2))
    day    = LONG(STRMID(head[8],27,2))
    heure  = LONG(STRMID(head[9],21,2)) ; !!!!!!!!!!!!! I ASSUME THE HOUR IS UT !!!!!!!!!!!
    minute = LONG(STRMID(head[9],24,2))
    seconde= LONG(STRMID(head[9],27,2))
    PRINT,'OBSERVATION DATE',day,month,year
    PRINT,'OBSERVATION TIME',heure,minute,seconde

   ; WE READ THE FIRST DARK FRAME
   ;-----------------------------------------
   
    READF,1,name
    PRINT,'1st DARK CURRENT= ',name
    dark1  = READFITS(base+name,head,/silent)
    dark1  = dark1[ix,*] + 0.d0 ; to suppress vertical central line
    dark1  = dark1[*,iy]        ; to suppress horizontal central line
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

   ; WE READ THE FILTERGRAMS
   ;-----------------------------------------

    FOR i=0,n-3 DO BEGIN
        READF,1,name
        PRINT,' '
        PRINT,'FILE READ= ',name
        q = READFITS(base+name,head,/silent)
        q = q[ix,*] + 0.d0      ; to suppress vertical central line
        q = q[*,iy]
        
        year   = LONG(STRMID(head[8],19,4))
        month  = LONG(STRMID(head[8],24,2))
        day    = LONG(STRMID(head[8],27,2))
        heure  = LONG(STRMID(head[9],21,2)) ; !!!!!!!!!!!!! I ASSUME THE HOUR IS UT !!!!!!!!!!!
        minute = LONG(STRMID(head[9],24,2))
        seconde= LONG(STRMID(head[9],27,2))
        time   = heure*3600.d0 + minute*60.d0 + seconde
        
        IF(i EQ 0) THEN IF (heure LT 12) THEN PRINT,'TIME=',time+24.*3600. ELSE PRINT,'TIME=',time
        
        header[iii*(n-2)+i,*] = head
        
        imx[*,*,iii*(n-2)+i]  = REBIN(q,nx,ny) ; we rebin the filtergram
        
        PRINT,'OBSERVATION DATE',day,month,year
        PRINT,'OBSERVATION TIME',heure,minute,seconde
        PRINT,'EXPOSURE TIME',head[68]
        PRINT,'EXPOSURE PROBLEM ?',TOTAL(imx[0:15,0:15,iii*(n-2)+i])
        
    ENDFOR
    
   ; WE READ THE SECOND DARK FRAME
   ;--------------------------------------------

    READF,1,name
    PRINT,'2nd DARK CURRENT= ',name
    dark2  = READFITS(base+name,head,/silent)
    dark2  = dark2[ix,*] + 0.d0 ; to suppress vertical central line
    dark2  = dark2[*,iy]        ; to suppress horizontal central line
    
    year   = LONG(STRMID(head[8],19,4))
    month  = LONG(STRMID(head[8],24,2))
    day    = LONG(STRMID(head[8],27,2))
    heure  = LONG(STRMID(head[9],21,2)) ; !!!!!!!!!!!!! I ASSUME THE HOUR IS UT !!!!!!!!!!!
    minute = LONG(STRMID(head[9],24,2))
    seconde= LONG(STRMID(head[9],27,2))
    timedark2   = heure*3600.d0 + minute*60.d0 + seconde
    dark2 = REBIN(dark2,nx,ny)
    PRINT,'WAVELENGTH NUMBER:',iii
    PRINT,'OBSERVATION DATE',day,month,year
    PRINT,'OBSERVATION TIME',heure,minute,seconde
    PRINT,'EXPOSURE TIME',head[68]
    PRINT,'LASER WAVELENGTH',head[73]
    
    FOR i=0,n-3 DO BEGIN

    ;   LINEAR INTERPOLATION AS SUGGESTED BY JESPER
        FOR j=0,ny-1 DO BEGIN
            FOR jj=0,nx-1 DO BEGIN
                imx[jj,j,iii*(n-2)+i] = imx[jj,j,iii*(n-2)+i] - INTERPOL([dark1[jj,j],dark2[jj,j]],[timedark1,timedark2],time)
            ENDFOR
        ENDFOR
        temp = REFORM(imx[*,*,iii*(n-2)+i])
        a    = WHERE(temp LT 0.d0)
        IF (a[0] NE -1) THEN temp[a]= 0.d0
        imx[*,*,iii*(n-2)+i]=temp
        
        
        TVIM,imx[*,*,iii*(n-2)+i],/scale,pcharsize=2
        
    ENDFOR

    calmode[*,*,iii] = calmode[*,*,iii] - dark1
    TVIM,calmode[*,*,iii],/scale,pcharsize=2


    tot = TOTAL(TOTAL(imx[*,*,iii*(n-2):(iii+1)*(n-2)-1],1),1)
    a   = WHERE(tot EQ MAX(tot))
    b   = WHERE(tot EQ MIN(tot))
    PRINT,'WAVELENGTH NUMBER',iii
    PRINT,'MAXIMUM INTENSITY',a
    PRINT,'VALUES',tot[a],tot[b]

    READ,pause
    
ENDFOR

CLOSE,1

; WE FREE MEMORY
dark1 = 0.d0
dark2 = 0.d0
q     = 0.d0
im    = 0.d0

SAVE,imx,calmode,header,FILE='SEQUENCE_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

RETURN,imx

END
