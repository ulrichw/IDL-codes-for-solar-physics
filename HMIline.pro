; produce the Fe I line profile for different center-to-limb distances
; using line profiles of R. Ulrich

PRO HMIline

distance = 20 ; center-to-limb distance in degrees

OPENR,1,'Ulrich_Fe_0.txt'
roger  = FLTARR(2,98)
READF,1,roger
CLOSE,1
nroger = 98
rlam0  = REFORM(roger[0,*])
rp0    = REFORM(roger[1,*])

OPENR,1,'Ulrich_Fe_45.txt'
READF,1,roger
CLOSE,1
nroger = 98
rlam45 = REFORM(roger[0,*])
rp45   = REFORM(roger[1,*])

OPENR,1,'Ulrich_Fe_60.txt'
READF,1,roger
CLOSE,1
nroger = 98
rlam60 = REFORM(roger[0,*])
rp60   = REFORM(roger[1,*])

rp     = FLTARR(98)

FOR i=0,97 DO rp[i]=INTERPOL([rp0[i],rp45[i],rp60[i]],[0.,45.,60.],distance)

READ,pause

END
