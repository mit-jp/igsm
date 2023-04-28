pro makeeppamap_write, scale

; read in country numbers
dir='/users/noelle/documents/mit07/aerosol/idlcode/countries_eppa/'
restore,  'countrytemplate'
afr=read_ascii(dir+'afr.csv', template=countrytemplate,  /verbose)
anz=read_ascii(dir+'anz.csv', template=countrytemplate,  /verbose)
asi=read_ascii(dir+'asi.csv', template=countrytemplate,  /verbose)
can=read_ascii(dir+'can.csv', template=countrytemplate,  /verbose)
chn=read_ascii(dir+'chn.csv', template=countrytemplate,  /verbose)
eet=read_ascii(dir+'eet.csv', template=countrytemplate,  /verbose)
eur=read_ascii(dir+'eur.csv', template=countrytemplate,  /verbose)
fsu=read_ascii(dir+'fsu.csv', template=countrytemplate,  /verbose)
idz=read_ascii(dir+'idz.csv', template=countrytemplate,  /verbose)
ind=read_ascii(dir+'ind.csv', template=countrytemplate,  /verbose)
jpn=read_ascii(dir+'jpn.csv', template=countrytemplate,  /verbose)
lam=read_ascii(dir+'lam.csv', template=countrytemplate,  /verbose)
mes=read_ascii(dir+'mes.csv', template=countrytemplate,  /verbose)
mex=read_ascii(dir+'mex.csv', template=countrytemplate,  /verbose)
row=read_ascii(dir+'row.csv', template=countrytemplate,  /verbose)
usa=read_ascii(dir+'usa.csv', template=countrytemplate,  /verbose)

;set grid type
resol=scale
onebyone=ctm_grid(ctm_type('GENERIC',  resol=resol))
if (scale eq 1) then begin
newmap=fltarr(360, 180)     
scale='1x1'
endif

if (scale eq 2) then begin
newmap=fltarr(144, 91)
scale='2x25'
endif

if (scale eq 4) then begin
newmap=fltarr(72, 43)
scale='4x5'
endif

for ii=0,  N_Elements(afr.numbers)-1L do begin
afrcells=find_cells_by_country(afr.numbers[ii], onebyone)
newmap=newmap+afrcells
endfor

openw, lun, 'eppa.afr.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(anz.numbers)-1L do begin
anzcells=find_cells_by_country(anz.numbers[ii], onebyone)
newmap=newmap+anzcells
endfor

openw, lun, 'eppa.anz.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(asi.numbers)-1L do begin
asicells=find_cells_by_country(asi.numbers[ii], onebyone)
newmap=newmap+asicells
endfor
openw, lun, 'eppa.asi.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(can.numbers)-1L do begin
cancells=find_cells_by_country(can.numbers[ii], onebyone)
newmap=newmap+cancells
endfor
openw, lun, 'eppa.can.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(chn.numbers)-1L do begin
chncells=find_cells_by_country(chn.numbers[ii], onebyone)
newmap=newmap+chncells
endfor

openw, lun, 'eppa.chn.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(eet.numbers)-1L do begin
eetcells=find_cells_by_country(eet.numbers[ii], onebyone)
newmap=newmap+eetcells
endfor
openw, lun, 'eppa.eet.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0


for ii=0,  N_Elements(eur.numbers)-1L do begin
eurcells=find_cells_by_country(eur.numbers[ii], onebyone)
newmap=newmap+eurcells
endfor
openw, lun, 'eppa.eur.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(fsu.numbers)-1L do begin
fsucells=find_cells_by_country(fsu.numbers[ii], onebyone)
newmap=newmap+fsucells
endfor
openw, lun, 'eppa.fsu.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0


for ii=0,  N_Elements(idz.numbers)-1L do begin
idzcells=find_cells_by_country(idz.numbers[ii], onebyone)
newmap=newmap+idzcells
endfor
openw, lun, 'eppa.idz.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(ind.numbers)-1L do begin
indcells=find_cells_by_country(ind.numbers[ii], onebyone)
newmap=newmap+indcells
endfor
openw, lun, 'eppa.ind.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(jpn.numbers)-1L do begin
jpncells=find_cells_by_country(jpn.numbers[ii], onebyone)
newmap=newmap+jpncells
endfor
openw, lun, 'eppa.jpn.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(lam.numbers)-1L do begin
lamcells=find_cells_by_country(lam.numbers[ii], onebyone)
newmap=newmap+lamcells
endfor
openw, lun, 'eppa.lam.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(mes.numbers)-1L do begin
mescells=find_cells_by_country(mes.numbers[ii], onebyone)
newmap=newmap+mescells
endfor
openw, lun, 'eppa.mes.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(mex.numbers)-1L do begin
mexcells=find_cells_by_country(mex.numbers[ii], onebyone)
newmap=newmap+mexcells
endfor
openw, lun, 'eppa.mex.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(row.numbers)-1L do begin
rowcells=find_cells_by_country(row.numbers[ii], onebyone)
newmap=newmap+rowcells
endfor
openw, lun, 'eppa.row.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

for ii=0,  N_Elements(usa.numbers)-1L do begin
usacells=find_cells_by_country(usa.numbers[ii], onebyone)
newmap=newmap+usacells
endfor
openw, lun, 'eppa.usa.'+scale, /get_lun,  /f77_unformatted
newmap=double(newmap)
writeu,  lun,  newmap
close,  lun
free_lun,  lun
newmap=0d0

end
