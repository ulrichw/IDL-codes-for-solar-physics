; THIS PROGRAM PRODUCES A 2D POWER SPECTRUM FROM A 3D FOURIER
; TRANSFORM OF A DATACUBE


PRO spectrum,data,spectre

nx = N_ELEMENTS(data[*,0,0])
nt = N_ELEMENTS(data[0,0,*])

spectre = FLTARR(nx/2,nt/2)

FOR t=0,nt/2-1 DO BEGIN
    FOR l=0,nx/2-1 DO BEGIN
        PRINT,t,l
        mask  = dist(nx,nx)
        FOR i = 0,nx-1 DO FOR j=0,nx-1 DO mask[i,j]=EXP(-(mask[i,j]-l)^2.d0/(2.d0*0.5d0^2.d0))
        mask  = mask/TOTAL(mask)
        spectre[l,t] = TOTAL(mask[*,*]*data[*,*,t])
    ENDFOR
ENDFOR

;tvim,alog10(spectre),xrange=[0,128./256./0.826*696.*2.*!pi],yrange=[0,256./512./60.*1000.]

END

