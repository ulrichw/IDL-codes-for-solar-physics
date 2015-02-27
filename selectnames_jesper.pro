pro selectnames_jesper,names,nfiles
dir='qqq'
read,dir,prompt='Directory name: '
if (dir eq '') then dir='.'
dir=dir+'/'
spawn,'ls -1 '+dir,list0
repeat begin
  name='qqq'
  read,name,prompt='First file: '
  if (name eq '') then name='i_.*\.fit'
  whit=where(stregex(list0,name,/boolean) ne 0,nhit)
  if (nhit eq 0) then begin
    print,'No matches'
  endif else begin
    list=list0(whit)
    if (nhit eq 1) then print,'Found file: '+list(0)
    if (nhit ge 2) then begin
      print,'Found multiple files:'
      print,list
    end
  end
endrep until (nhit eq 1)
ifirst=whit(0)

repeat begin
  name='qqq'
  read,name,prompt='Last file: '
  if (name eq '') then name='i_.*\.fit'
  whit=where(stregex(list0,name,/boolean) ne 0,nhit)
  if (nhit eq 0) then begin
    print,'No matches'
  endif else begin
    list=list0(whit)
    if (nhit eq 1) then print,'Found file: '+list(0)
    if (nhit ge 2) then begin
      print,'Found multiple files:'
      print,list
    end
  end
endrep until (nhit eq 1)
ilast=whit(0)

names=dir+list0(ifirst:ilast)
nfiles=ilast-ifirst+1
print,'Found',nfiles,' files'

end

