RESTORE,'FeIline.bin'
plot,asin(findgen(1880)/1880.)*180./!dpi,res[0:1880],xrange=[0,90],xst=1,yrange=[0,16000],yst=1,thick=3; continuum (in DN)
plot,asin(findgen(1880)/1880.)*180./!dpi,res2[0:1880],xrange=[0,90],xst=1,yrange=[100,140],yst=1,thick=3; linewidth (FWHM in mA)

plot,asin(findgen(1880)/1880.)*180./!dpi,res3[0:1880]/res[0:1880],xrange=[0,90],xst=1,yrange=[0.3,0.6],yst=1,thick=3; linedepth (in DN)
plot,asin(findgen(1880)/1880.)*180./!dpi,res4[0:1880],xrange=[0,90],xst=1,yrange=[-50,50],yst=1,thick=3; magnetic field (Gauss)
oplot,[0,90],[0,0]
plot,asin(findgen(1880)/1880.)*180./!dpi,res5[0:1880],xrange=[0,90],xst=1,yst=1,thick=3; Doppler velocity
