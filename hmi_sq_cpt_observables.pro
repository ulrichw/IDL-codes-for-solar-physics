PRO hmi_sq_cpt_observables

; to read the data tkaen on January 30, 2008 at Goddard by Jesper

d = INTARR(4096,4096,240)

For i=0,239 do d[*,*,i]=mrdfits('/tmp20/schou/hmi080130/00'+STRTRIM(STRING(i+1003509),1)+'.fits')

read,pause


END
