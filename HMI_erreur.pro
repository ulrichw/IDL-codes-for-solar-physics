FUNCTION line,long
npar  = 3
nlamp = 1223
lamp  = -0.3D0+DINDGEN(nlamp)*0.6D/(nlamp-1)
nroger = 98
lines  = FLTARR(3,nlamp)
lon    = [0.5,0.707107,1.0] ;longitudes
lineint= FLTARR(nlamp) ; interpolated line at the specific longitude long

FOR i=0,2 DO BEGIN
    IF i EQ 0 THEN OPENR,1,'Ulrich_Fe_60.txt' ; Line profile from Roger Ulrich
    IF i EQ 1 THEN OPENR,1,'Ulrich_Fe_45.txt'
    IF i EQ 2 THEN OPENR,1,'Ulrich_Fe_0.txt'
    roger  = FLTARR(2,98)
    READF,1,roger
    CLOSE,1

    rlam   = REFORM(roger(0,*))
    rp     = REFORM(roger(1,*))
    rp1    = rp/INTERPOL([rp(0),rp(nroger-1)],[rlam(0),rlam(nroger-1)],rlam)
    lines[i,*]   = INTERPOL(rp,rlam,lamp)
    lines[i,*]   = lines[i,*]/INTERPOL([lines[i,0],lines[i,nlamp-1]],[lamp(0),lamp(nlamp-1)],lamp)
ENDFOR

FOR i=0,nlamp-1 DO lineint[i]=INTERPOL(lines[*,i],lon,long)

SAVE,lineint,lamp,FILE='interpolatedsolarline.bin'

RETURN,0

END



PRO HMI_erreur

RESTORE,'RESULTS/RESULTS_listLaser070108_40939_256.BIN'
erreur=FLTARR(256,256,3)
a=0.0
b=0.0
c=0.0

WINDOW,0,RETAIN=2
FOR i=0,255 DO BEGIN
    PRINT,i
    TVIM,erreur[*,*,0],/scale
    FOR j=0,255 DO BEGIN
        long=1.0-SQRT((i-128.0)^2.0+(j-123.0)^2.)/120.0
        IF(Bg0[i,j,0] NE 0.0) THEN BEGIN 
            d=line(long)
            HMI,REFORM(Phig0[i,j,*])-[80.02,178.23,-52.17]*!pi/180.,REFORM(Bg0[i,j,*]),a,b,c,long
            erreur[i,j,*] = [a,b,c]
        ENDIF
        PRINT,erreur[i,j,*],long,REFORM(Phig0[i,j,*])-[80.02,178.23,-52.17]*!pi/180.
    ENDFOR
ENDFOR

SAVE,erreur,FILE='erreur_line2.bin' ;'erreur.bin'
; some results in:
;erreur_line_j=123.bin
;erreur_line_i=128.bin
READ,pause

END
