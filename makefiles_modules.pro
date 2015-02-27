PRO makefiles_modules

d=fltarr(256,256,4)
openr,1,'non_tunable_phases_710660_June09_cal_256_2.bin'
readu,1,d
close,1
d2=fltarr(64,64,4)
for i=0,3 do d2[*,*,i]=REBIN(REFORM(d[*,*,i]),64,64)
openw,1,'non_tunable_phases_710660_June09_cal_64_2.bin'
writeu,1,d2
close,1 

openr,1,'non_tunable_contrasts_710660_June09_cal_256_2.bin'
d=fltarr(256,256,4)
readu,1,d
close,1
d2=fltarr(64,64,4)
for i=0,3 do d2[*,*,i]=REBIN(REFORM(d[*,*,i]),64,64)
openw,1,'non_tunable_contrasts_710660_June09_cal_64_2.bin'
writeu,1,d2
close,1 

openr,1,'tunable_contrasts_710660_June09_cal_256.bin'
d=fltarr(256,256,3)
readu,1,d
close,1
d2=fltarr(64,64,3)
for i=0,2 do d2[*,*,i]=REBIN(REFORM(d[*,*,i]),64,64)
openw,1,'tunable_contrasts_710660_June09_cal_64.bin'
writeu,1,d2
close,1

END
