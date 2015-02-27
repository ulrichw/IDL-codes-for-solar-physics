; to azimutahlly average a 2D or 3D datacube


FUNCTION azimutaverage2D,center,data

nx   = N_ELEMENTS(data[*,0])
mask = SHIFT(DIST(nx,nx),center[0],center[1])
sep  = FLTARR(nx/2)

FOR ii=0,nx/2-1 DO BEGIN
    print,ii
    mask2=EXP(-(mask-ii)^2.d0)  ;!thickness
    tot  = TOTAL(mask2)
    mask2=mask2/FLOAT(tot)
    sep[ii] = TOTAL(mask2*REFORM(data))
ENDFOR    

    RETURN,sep
END

FUNCTION azimutaverage,center,data

;nx   = N_ELEMENTS(data[*,0])
nx   = N_ELEMENTS(data[*,0,0])
nt   = N_ELEMENTS(data[0,0,*])
mask = SHIFT(DIST(nx,nx),center[0],center[1])
sep  = FLTARR(nx/2,nt)

FOR ii=0,nx/2-1 DO BEGIN
    print,ii
    mask2=EXP(-(mask-ii)^2.d0)  ;!thickness
    tot  = TOTAL(mask2)
    mask2=mask2/FLOAT(tot)
    FOR jj=0,nt-1 DO BEGIN
        sep[ii,jj] = TOTAL(mask2*REFORM(data[*,*,jj]))
    ENDFOR
ENDFOR    

    RETURN,sep
END
