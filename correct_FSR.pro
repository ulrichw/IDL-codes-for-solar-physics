PRO correct_FSR

vel1=-2307.55d0  ;6 velocity 1st detune
vel2=2285.40d0  ;31 veloctiy 2nd detune 

;NB
phase1 = -132.310d0; in degrees
phase2 = -128.930d0
;WB
phase3= 11.77d0
phase4= 14.02d0
;E1
phase5= -141.238d0
phase6= -140.920d0


dlamdv = 2.059205672212074294d-5
dl=dlamdv*(vel2-vel1)

FSRNB=0.172d0-0.0010576d0
FSRWB=0.344d0-0.00207683d0
FSRE1=0.693d0+0.000483467d0

dFSRNB=(phase2-phase1)/360.d0/dl*FSRNB^2.d0
dFSRWB=(phase4-phase3)/360.d0/dl*FSRWB^2.d0
dFSRE1=(phase6-phase5)/360.d0/dl*FSRE1^2.d0



PRINT,dFSRNB,FSRNB-dFSRNB
PRINT,dFSRWB,FSRWB-dFSRWB
PRINT,dFSRE1,FSRE1-dFSRE1



END
