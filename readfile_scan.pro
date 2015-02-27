FUNCTION readfile_scan,header,nx,base,list
; THE ROUTINE ROTATION.PRO MUST BE COMPILED PRIOR
; TO THIS CODE
; ROTATION.PRO RETURNS THE ROTATION ANGLE DUE TO THE HELIOSTAT
; suninfo PROVIDES THE p-angle

n      = 13 ; numbers of filtergrams: the 1st and last ones are just dark current
nkey   = 91
header = STRARR(n-2,nkey)

imx    = FLTARR(nx,nx,n-2)
name   = 'qqq'
ix     = [2098-2047+INDGEN(2048),2101+INDGEN(2048)]
iy     = [2067-2047+INDGEN(2048),2132+INDGEN(2048)]

WINDOW,0,retain=2,xsize=800,ysize=600
!p.multi=0

; WE READ THE 2 DARK CURRENTS
;-----------------------------------------
;OPENR,1,'list'
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
dark2 = REBIN(dark2,nx,nx)
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
    header[i,*] = head

  ; WE REMOVE THE DARK CURRENT
    q      = im[ix,*] + 0.d0  ; to suppress vertical central line
    q      = q[*,iy]          ; to suppress horizontal central line

    imx[*,*,i] = REBIN(q,nx,nx)       ; we rebin the filtergram
    TVIM,imx[*,*,i],/scale
    FOR j=0,nx-1 DO BEGIN
        FOR jj=0,nx-1 DO BEGIN
            imx[jj,j,i] = imx[jj,j,i] - INTERPOL([dark1[jj,j],dark2[jj,j]],[timedark1,timedark2],time)
        ENDFOR
    ENDFOR
    temp = REFORM(imx[*,*,i])
    a    = WHERE(temp LT 0.d0)
    IF (a[0] NE -1) THEN temp[a]= 0.d0
    imx[*,*,i]=temp


    PRINT,'OBSERVATION DATE',day,month,year
    PRINT,'OBSERVATION TIME',heure,minute,seconde
    PRINT,'EXPOSURE TIME',head[68]

ENDFOR
CLOSE,1

SAVE,imx,dark1,dark2,FILE='SEQUENCE_SCAN_'+STRTRIM(list,1)+'_'+STRTRIM(STRING(LONG(nx)),1)+'.BIN'

RETURN,imx

END
