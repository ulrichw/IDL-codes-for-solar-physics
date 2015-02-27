; PRODUCE A COTUNE TABLE FOR SPECIFIC PHASES OF THE TUNABLE ELEMENTS

PRO cotunetable

NBphase  = -122.28;-125.;-130.00098;-131.524;-106.26;-22.61 ; in degrees !!!
WBphase  =  51.;39.;34.;21.799087;42.35;-87.05
E1phase  = -137.28;-136.7;-138.339;-138.924216;-140.32;-140.19

HCMNB    = ROUND(-NBphase/6.0)+60
HCMWB    = ROUND(-WBphase/6.0)+60
HCME1    = ROUND( E1phase/6.0)+60


table=FLTARR(4,20)
table= [ [-30,0,0,0], $
         [-27,-6,0,-12], $
         [-24,-12,0,-24], $
         [-21,-18,0,24], $
         [-18,-24,0,12], $
         [-15,-30,0,0], $
         [-12,24,0,-12], $
         [-9,18,0,-24], $
         [-6,12,0,24], $
         [-3,6,0,12], $
         [0,0,0,0], $
         [3,-6,0,-12], $
         [6,-12,0,-24], $
         [9,-18,0,24], $
         [12,-24,0,12], $
         [15,-30,0,0], $
         [18,24,0,-12], $
         [21,18,0,-24], $
         [24,12,0,24], $
         [27,6,0,12] ]

FOR i=0,19 DO table[*,i] = table[*,i] + [HCME1,HCMWB,0,HCMNB]

PRINT,table


READ,pause

END
