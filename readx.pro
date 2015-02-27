function readx,filename

nin=0l
spawn,'wc -l '+filename,help
reads,help,nin
ntot=0l
spawn,'wc -w '+filename,help
reads,help,ntot
if ((ntot mod nin) eq 0) then begin
  ncol=ntot/nin
  q=dblarr(ncol,nin)
  openr,37,filename
  readf,37,q
  close,37
end

return,q

end

