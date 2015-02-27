; test the artifacts that can appear when focus blocks are changed

PRO HMI_artifactcheck

nx=256
RESTORE,'CPT/SEQUENCE_listLamp070903_539006_512.BIN'
i1=images
RESTORE,'CPT/SEQUENCE_listLamp070903_539074_512.BIN'
i2=images
RESTORE,'CPT/SEQUENCE_listLamp070903_539142_512.BIN'
i3=images
RESTORE,'CPT/SEQUENCE_listLamp070903_539210_512.BIN'
i4=images
RESTORE,'CPT/SEQUENCE_listLamp070903_539278_512.BIN'
i5=images
RESTORE,'CPT/SEQUENCE_listLamp070903_539346_512.BIN'
i9=images
RESTORE,'CPT/SEQUENCE_listLamp070903_539414_512.BIN'
i13=images

flat = READFITS('~richard/public_html/flatfield_front.fits') ; flat field
flat1= REBIN(flat[*,*,0],512,512,1)
flat13=REBIN(flat[*,*,6],512,512,1)
temp1= REFORM(REBIN(i1,1,1,34))
temp13=REFORM(REBIN(i13,1,1,34))
FOR i=0,33 DO i1[*,*,i] =i1[*,*,i]/flat1
FOR i=0,33 DO i13[*,*,i]=i13[*,*,i]/flat13
tvim,(i1[*,*,6]/temp1[6]-i13[*,*,6]/temp13[6])/(i13[*,*,6]/temp13[6]),/scale,range=[-0.011,0.011],barwidth=0.5,tit='CF 1-13, position 6',stit='Re lative difference'

READ,pause

END
