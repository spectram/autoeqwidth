;---------------------------------------------------------------------------------------------------
;----------------------------------------Autoeqwidth v1.2-------------------------------------------
;---------------------Written by Sriram Sankar and Jayadev Pradeep 2019-----------------------------
;-----------Some parts have been adopted from codes written by Dr Anand Narayanan-------------------
;---------------------------------------------------------------------------------------------------
;The code calculates equivalent-width and Apparent column density for a list of ions 
;For non-detections the linear relation of CoG is used to estimate the upper limits. 
;The code does not account for continuum fitting uncertainties.
;
;Input file format: 'Ion' 'wavelength' 'oscillator strength' 'gamma-value' 'Detection Flag'
;1 - Detection & 0 - Non-detection - [Optional input]
;
;Equations for equivalent width calculation and column density measurement for the case of detections 
;have been adopted from the thesis of Chris Churchill (1997) - The Low Ionization Gaseous Content in 
;Intermediate Redshift Galaxies: http://astronomy.nmsu.edu/cwc/Research/thesis/Pages/ch2.html. 
;Labelled as 'CC Eq: 2.xx'

;Curve of Growth Linear regime relation for the case of non-detections adopted from: 
;Meiksin, Avery A. “The Physics of the Intergalactic Medium.” Reviews of Modern Physics 81, no. 4 
;(October 5, 2009): 1405–69. https://doi.org/10.1103/RevModPhys.81.1405. 
;Labelled as 'MA Eq:xx'
;Note that the order of magnitude of the constant of proportionality has been changed to account for 
;values input in cm and Angstrom.
;---------------------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------------------
;-----------If you wish to thank me or report bugs mail me at sriram10sankar@gmail.com--------------
;---------------------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------------------
;--------------------------------------------Features-----------------------------------------------
;Scale the continuum, apply velocity offsets, generate .csv table and Latex table of measurements,
;choose limits by clicking on the plot, and generate input file for system plot or for later runs...
;---------------------------------------------------------------------------------------------------
;------------------------------------------Aug 01 2019----------------------------------------------
;---------------------------------------------------------------------------------------------------
;Changes after release
;saving v1 and v2 as long in input file
;Significance level, detection prediction and option to flip detection flag
;FIxed lname bug for wavelength outside coverage
;Bug fix in formula to accurately measure equivalent width of saturated features - 13-Dec-2019
;/detflag added, automatic run from previous run file added and other minor bug fixes - 5-Jan-2020
;chknew1 declared to better differentiate between new runs and atmsip runs. {vl not declared error} - 26-Apr-2021
;
;
;
;


PRO autoeqwidth,ipfile=ipfile,ionlist=ionlist,z=z,detflag=dflag

if NOT KEYWORD_SET(ipfile) then begin
	 print,'Syntax - ' + $
        'autoeqwidth, ipfile=xyz.asc (within quotes), ionlist=xyz.dat (within quotes), z=00.00,/detflag'
      return
endif 

if NOT KEYWORD_SET(ionlist) then ionlist='atomsnew.dat'

print,'-------- Welcome to autoeqwidth v1.2 2019 --------'

close,1,2,3
c=299792.458 ;speed of light in vacuum in km/s

;Default range
X1d=0
X2d=0
;Previous limit set to default for first time
X1=X1d
X2=X2d

j=0L
k=0
vmax=600 ;km/s
v_xmax=vmax
v_xmin=-1.0*vmax

cscale=0
ni=0 ; Number of iterations
chkc=0 ;Continuum scaling chk
chkv=0 ;Velocity shift chk
chkdf=0 ;Flip detection flag chk

readcol,ipfile,wav,flux,err,cont,/SILENT
i1=n_elements(wav)
contnew=cont	

ion=strarr(file_lines(ionlist))
xyz=fltarr(file_lines(ionlist))
t3=strarr(file_lines(ionlist))
lamdao=fltarr(file_lines(ionlist))
f=fltarr(file_lines(ionlist))
lamda=fltarr(file_lines(ionlist))
wl=fltarr(file_lines(ionlist))
wl_e=fltarr(file_lines(ionlist))
nf=fltarr(file_lines(ionlist))
nf_e=fltarr(file_lines(ionlist))
range=strarr(file_lines(ionlist))
v1=fltarr(file_lines(ionlist))
v2=fltarr(file_lines(ionlist))
DUM_ar=fltarr(1)
adet=strarr(file_lines(ionlist))

dum=0
k=0

;------------------------------------------------------------------------------------------------
chkatm=0
chknew=0  ;re-measure
chknew1=0 ;new measurements
atmsip='atomsip'+strtrim(string(z),2)+'.dat'
atmip=file_search(atmsip)
if atmip eq atmsip then begin
print,'Input file from previous run found: ',atmsip

ysno=''
READ,ysno,PROMPT='Do you wish to re-generate output files directly? (y/any key): '
chkatm=STRCMP(ysno,'y',/FOLD_CASE)
endif

if chkatm eq 1 then begin

;------------------------------Reading autoeq i/p file---------------------------------
j=0L

vl=fltarr(file_lines(atmsip))
vr=fltarr(file_lines(atmsip))

openr,2,atmsip
 while NOT EOF(2) do begin
       var1=' '
	var2=''
	var3=''
	var4=''
	var5=''
	var6=''
	var7=''
	readf,2,var1,var2,var3,var4,var5,var7,var6,format='(A14,A13,A15,A10,A12,A4,A4)' 
	;If the previous one gives an error try - '(A14,A12,A19,A4,A12,A4,A4)'
    	ion[j]=var1
    	lamdao[j]=var2
    	f[j]=var3
    	vl[j]=long(var4)
    	vr[j]=long(var5)
	t3[j]=var6
	print,ion[j],' | ',lamdao[j],' | ',f[j],' | ',vl[j],' | ',vr[j],' | ',t3[j]
   j=j+1 ;j = Number of ions
 endwhile
close,2

det=fix(t3)
sm=total(det)


while (k lt j) do begin
	lamda[k]=lamdao[k]*(1+z)
	lname=ion[k]
	if (lamda[k] lt wav[0]) or (lamda[k] gt wav[i1-5]) then begin
	print,strtrim(ion[k],2)+' '+strtrim(string(nint(lamdao[k])),2)+' is outside wavelength coverage.'
	dum=dum+1
	DUM_ar=[DUM_ar,k]
  	endif else begin

v_wav=((wav-lamda[k])/lamda[k])*c

;--------------------------------Equivalent width calculation------------------------------------
ew=0
ew_err_sqsum=0
ew_err=0

w_X1=lamda[k]*(1+vl[k]/c)
p1=min(abs(wav-w_X1), index1)

w_X2=lamda[k]*(1+vr[k]/c)
p2=min(abs(wav-w_X2), index2)

if index2 lt index1 then begin
temp=index1
tempX=X1
index1=index2
X1=X2
index2=temp
X2=tempX
endif

for q=index1,index2 do begin

if flux[q] gt 0.0000 then begin
val=1-(flux[q]/cont[q]) 
endif else begin 
val=1 
endelse
w_diff=wav[q+1]-wav[q]
if q eq index2 then w_diff=wav[q]-wav[q-1]
ew=ew+val*w_diff 
tx=((err[q]^2)*(-1/(cont[q]))^2)						;CC Eq: 2.36 - Dimenssion corrected
ew_err_sqsum = ew_err_sqsum + tx*w_diff^2
endfor 

ew_err=sqrt(ew_err_sqsum)/(1+z)
ew=ew/(1+z)


;---------------------------------------Significance---------------------------------------------
;------------------------------------------------------------------------------------------------

if KEYWORD_SET(dflag) then begin
sign=ew/ew_err
sign=float(round(sign*100)/100.0d)
significance=strmid(strtrim(string(long(sign)),2),0,5)+' sigma '
if sign ge 3 then adet[k] = 1 else adet[k] = 0
if det[k] ne adet[k] and sm ne 0 then begin
print,''
print,'**Alert**'
print,'Mismatch noted between derived detection flag and input detection flag for ',strtrim(ion[k],2)
Print,'In certain cases of contamination, it is better to measure the upper limit using AOD method.'
print,'But be careful to modify the output for the ion and make it an upper limit'  
print,''
Print,'Please check the feature and enter the detection flag manually' 
print,'0 - Non-detection'
print,'1 - Detection'
print,''
if adet[k] eq 1 then print,'It could be a ',significance,'detection for the selected range'
if adet[k] eq 0 then print,'It could be a non-detection. Significance (< 3 sigma) = ',significance
print,''
detf=0
READ,detf,PROMPT='What will it be?: '
print,''
det[k]=detf
endif else begin
det[k]=adet[k]
sm=1
endelse
endif

;--------------------------------------------AOD-------------------------------------------------
tauv=alog(cont/flux)
nav=3.768e14*tauv/(f[k]*lamdao[k]) 
nav_12=nav/1e12
;-----------------------------------------Detection----------------------------------------------
if(det[k] ne 0.00) then begin
na=0
tau=0
sigma_tau_sumsq=0

for q=index1,index2 do begin
	if flux[q] gt 0 then begin
		v_diff=v_wav[q+1]-v_wav[q]
		if q eq index2 then v_diff = v_wav[q]-v_wav[q-1]
		na=na+(nav[q]*v_diff)
		tau=tau+(tauv[q]*v_diff)
		current = (err[q]^2)*((-1.0/flux[q])^2)*v_diff 		
		sigma_tau_sumsq= sigma_tau_sumsq + current                			
	endif 
endfor 

logna=alog10(na)
sigma_tau=sqrt(sigma_tau_sumsq)
sigma_na = na * (sigma_tau^2) / tau
sigma_logna=0.4343*sigma_na/na

wl[k]=ew*1e3
wl_e[k]=ew_err*1e3

nf[k]=logna
nf_e[k]=sigma_logna 
range[k]='['+strtrim(string(Long(vl[k])),2)+','+strtrim(string(Long(vr[k])),2)+']'  ;Just trimming white spaces and decimals
v1[k]=vl[k]
v2[k]=vr[k]

endif else begin
;----------------------------------------Non-Detection-------------------------------------------
wlamda= 3*ew_err
nu=wlamda/((8.855e-21)*f(k)*((lamdao[k])^2))
nuf=alog10(nu)

wl[k]=0.000
wl_e[k]= wlamda*1e3
nf[k]=0.000
nf_e[k]= nuf
range[k]='['+strtrim(string(Long(vl[k])),2)+','+strtrim(string(Long(vr[k])),2)+']'  ;Just trimming white spaces and decimals
v1[k]=vl[k]
v2[k]=vr[k]
endelse

endelse
k=k+1
ni=ni+1
endwhile

if ni ne 0 and dum ne 0 then DUM_ar=DUM_ar[1:*]

print,'---------------------------------------------------'
print,ni-dum,'  measurements made' ;'Number of iterations: ',ni
print,'---------------------------------------------------'
print,''

endif ;chkatm


;------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------
ysno=''

if chkatm eq 1 then begin
READ,ysno,PROMPT='Do you wish to re-measure specific ions manually? (if y, make sure that the ionlist and the input file from the previous run have the same index values for the ions): '
chknew=STRCMP(ysno,'y',/FOLD_CASE)
endif else begin
READ,ysno,PROMPT='Do you wish to measure specific ions manually?: '
chknew1=STRCMP(ysno,'y',/FOLD_CASE)
endelse

if chknew1 eq 1 or chknew eq 1 then begin
ni=0
;------------------------------------------------------------------------------------------------
j=0L
openr,1,ionlist
  while NOT EOF(1) do begin
        var1=' '
	var2=1.0D
	var3=1.0D
	var4=1.0D
	var5=''
	readf,1,var1,var2,var3,var4,var5,format='(A11,1x,F9.4,1x,F9.7,1x,E8.5,1x,A2)'
    	ion[j]=var1
    	lamdao[j]=var2
    	f[j]=var3
    	xyz[j]=var4
    	t3[j]=var5
    j=j+1 ;j = Number of ions
  endwhile
close,1

det=fix(t3)

sm=total(det)
yn=''
nil='-'
print,''
print,'If the detection flag has not been added to the ionlist, you will encounter type conversion errors.'
print,'No need to be alarmed you can choose the detection flag while measuring each ion.'
print,''

READ,yn,PROMPT='Do you wish to verify if the ionlist has been read correctly? (y/any key): '
chk=STRCMP(yn,'y',/FOLD_CASE)
;if chkatm eq 1 then det=dettemp

if chk eq 1 then begin
print,'Ion | Rest-wavelength | Oscillator strength | detection flag'
k=0
while k lt j do begin
if sm eq 0 then print,ion[k],' | ',lamdao[k],' | ',f[k],' | ',nil
if sm ne 0 then print,ion[k],' | ',lamdao[k],' | ',f[k],' | ',det[k]
k=k+1
endwhile
endif
;------------------------------------------------------------------------------------------------
dum=0
k=0
DUM_ar=fltarr(1)

print,''
print,'*********** PRO TIP ***********'
print,'For best flow of work, align the terminal to the left (58x37) and set it to be always visible on workspace.'
print,''
print,'Hit enter when you are ready.'
enter=get_kbrd()

if chknew eq 1 then begin
kt=0
while (kt lt j) do begin
print,kt+1,' --- ',ion[kt]
kt=kt+1
endwhile
Print,''
READ,knew,PROMPT='Enter the index value of the ion to be measured: '
k=knew-1
if chknew eq 1 then begin
X1d=vl[k]
X2d=vr[k]
endif 
endif

endif

while (k lt j) do begin
	lamda[k]=lamdao[k]*(1+z)
	lname=ion[k]
	outsidecoverage=0
	if (lamda[k] lt wav[0]) or (lamda[k] gt wav[i1-5]) then begin
	print,strtrim(ion[k],2)+' '+strtrim(string(nint(lamdao[k])),2)+' is outside wavelength coverage.'
	dum=dum+1
	DUM_ar=[DUM_ar,k]
	outsidecoverage=1
  	endif else begin
dv=0
dx=0
v_wav=0
cscale=0
cont=contnew
v_wav=((wav-lamda[k])/lamda[k])*c

;------------------------------------Operation on each ion---------------------------------------
;------------------------------------------------------------------------------------------------

;--------------------------------------Continuum scaling-----------------------------------------
if ni ne 0 AND chkc eq 1 then begin
cscale=1.0
print,''
print,'Continuum of the order of ',strtrim(string(cont[25]),2)
print,''
READ,cscale,PROMPT='Enter value to scale the continuum (value will be multiplied by 1e-15): '
print,'Scaling the continuum by: ',cscale
contnew=cont
cont=cont+(cscale*1e-15)
endif
;---------------------------------------Shift Spectrum--------------------------------------------
if ni ne 0 AND chkv eq 1 then begin
dv=1.0
READ,dv,PROMPT='Enter the shift to be applied in km/s: '
print,'Shifting the spectrum by: ',strtrim(string(dv),2)
dx=(dv/c)*lamda[k]
wavold=wav
wav=wavold+dx
v_wav=((wav-lamda[k])/lamda[k])*c
endif
;------------------------------------------------------------------------------------------------
w_xmin=lamda[k]*(1+v_xmin/c)
w_xmax=lamda[k]*(1+v_xmax/c)

;------------------------------------------------------------------------------------------------
;lname=ion[k]
print,''
print,'*****************************'
print,'Analysing now: ',lname
print,'*****************************'
print,''

;------------------------------------------------------------------------------------------------
norm_flux=flux/cont
norm_cont=cont/cont
norm_err=err/cont
;----------------------------------------------Plot----------------------------------------------
!MOUSE.button = 1

cgdisplay, 1000, 800 
ERASE

cgplot,v_wav, norm_flux, xtitle=textoidl('!6 Velocity (km/s)'), ytitle=textoidl('!6 Normalized Flux'),$
xstyle=1,charsize=1.5,ycharsize=1.2,xcharsize=1.2,$
yrange=[0.0,1.5],xrange=[v_xmin,v_xmax],/ystyle,PSym=10,xthick=5,ythick=5,thick=2

cgoplot, v_wav, norm_err,PSym=10, color='blue', thick=1
cgoplot, v_wav, norm_cont,PSym=10,linestyle=2, color='green', thick=3
vline, [0.0],/NOERASE,linestyle=2, thick=3
cgtext, v_xmin+0.13*(v_xmax-v_xmin),0.2, textoidl(lname), charthick=2.0, color='red', charsize=2 
obs_wav_text = '!6 \lambda ='+strtrim(string(lamda[k]),2)+'A^{o}'
cgtext, v_xmax-0.35*(v_xmax-v_xmin),0.2, textoidl(obs_wav_text), charthick=2.0, color='green', charsize=2

vline, [X1d], thick=2, linestyle=2
vline, [X2d], thick=2, linestyle=2

;------------------------------------------------------------------------------------------------
;---------------------------------------Default limits-------------------------------------------

input=''

print,'---------------------------------------------------'
Print,'d  - default  limits ',strtrim(string(Long(X1d)),2),' - ',strtrim(string(Long(X2d)),2)
Print,'o  - previous limits ',strtrim(string(Long(X1)),2),' - ',strtrim(string(Long(X2)),2)
Print,'e  - to enter range ' 
Print,'ds - to set new default limits and measure ' 
Print,'s  - to click and select '
Print,'Enter any other key to exit '
print,'---------------------------------------------------'
Print,''
READ,input,PROMPT='What will it be?: '
Print,''
chk=STRCMP(input,'d',/FOLD_CASE)
chk2=STRCMP(input,'o',/FOLD_CASE)
chk3=STRCMP(input,'e',/FOLD_CASE)
chk4=STRCMP(input,'s',/FOLD_CASE)
chk5=STRCMP(input,'ds',/FOLD_CASE)

if chk ne 1 AND chk2 ne 1 AND chk3 ne 1 AND chk4 ne 1 and chk5 ne 1 THEN begin
READ,yn,PROMPT='**ALERT** Do you really want to exit? (y/n): '
chke=STRCMP(yn,'y',/FOLD_CASE)
if chke eq 1 THEN begin
break 
endif else begin 
k=k
input=''
print,''
print,'---------------------------------------------------'
Print,'d  - default  limits ',strtrim(string(Long(X1d)),2),' - ',strtrim(string(Long(X2d)),2)
Print,'o  - previous limits ',strtrim(string(Long(X1)),2),' - ',strtrim(string(Long(X2)),2)
Print,'e  - to enter range ' 
Print,'ds - to set new default limits and measure ' 
Print,'s  - to click and select '
Print,'Enter any other key to exit '
print,'---------------------------------------------------'
Print,''
READ,input,PROMPT='What will it be?: '
Print,''
chk=STRCMP(input,'d',/FOLD_CASE)
chk2=STRCMP(input,'o',/FOLD_CASE)
chk3=STRCMP(input,'e',/FOLD_CASE)
chk4=STRCMP(input,'s',/FOLD_CASE)
chk5=STRCMP(input,'ds',/FOLD_CASE)
if chk ne 1 AND chk2 ne 1 AND chk3 ne 1 AND chk4 ne 1 and chk5 ne 1 THEN break
endelse
endif

if chk eq 1 THEN BEGIN
X1=X1d
X2=X2d
endif 

if chk5 eq 1 THEN BEGIN
Print,''
read, X1d, prompt= 'Enter new default lower limit of Velocity in km/s: '
vline, [X1d],/NOERASE, thick=2, linestyle=2
read, X2d, prompt= 'Enter new default upper limit of Velocity in km/s: '
vline, [X2d],/NOERASE, thick=2, linestyle=2
Print,''
X1=X1d
X2=X2d
endif

if chk3 eq 1 then begin
Print,''
read, X1, prompt= 'Enter the lower limit of Velocity in km/s: '
vline, [X1],/NOERASE, thick=2, linestyle=2
read, X2, prompt= 'Enter the upper limit of Velocity in km/s: '
vline, [X2],/NOERASE, thick=2, linestyle=2
Print,''
endif 
;------------------------------------------------------------------------------------------------
;---------------------------------------Manual limits--------------------------------------------

if chk4 eq 1 then begin
Print,''
print,'Click to choose the index manually'
Print,''
print,'---------------------------------------------------'
Print,'Right click to exit'
Print,'Left Click to choose the range'
Print,'Middle Mouse to enter the final range'
print,'---------------------------------------------------'
Print,''
cont_tester=0
ERASE
WHILE (!MOUSE.button NE 4) DO BEGIN 

cgdisplay, 1000, 800 

cgplot,v_wav, norm_flux, xtitle=textoidl('!6 Velocity (km/s)'), ytitle=textoidl('!6 Normalized Flux'),$
xstyle=1,charsize=1.5,ycharsize=1.2,xcharsize=1.2,$
yrange=[0.0,1.5],xrange=[v_xmin,v_xmax],/ystyle,PSym=10,xthick=5,ythick=5,thick=2

cgoplot, v_wav, norm_err,PSym=10, color='blue', thick=1
cgoplot, v_wav, norm_cont,PSym=10,linestyle=2, color='green', thick=3
vline, [0.0],/NOERASE,linestyle=2, thick=3
cgtext, v_xmin+0.13*(v_xmax-v_xmin),0.2, textoidl(lname), charthick=2.0, color='red', charsize=2 
obs_wav_text = '!6 \lambda ='+strtrim(string(lamda[k]),2)+'A^{o}'
cgtext, v_xmax-0.35*(v_xmax-v_xmin),0.2, textoidl(obs_wav_text), charthick=2.0, color='green', charsize=2

if cont_tester eq 0 then begin
cont_tester=1
;ERASE
CONTINUE
endif

if(!MOUSE.button EQ 4) then break

CURSOR, X1, Y1, /DOWN, /DATA

if(!MOUSE.button EQ 2) then begin

read, X1, prompt= 'Enter the lower limit of Velocity in km/s: '
vline, [X1],/NOERASE, thick=2, linestyle=2
read, X2, prompt= 'Enter the upper limit of Velocity in km/s: '
vline, [X2],/NOERASE, thick=2, linestyle=2

break
endif else begin

CURSOR, X1, Y1, /DOWN, /DATA
print,'Lower limit: ', X1
vline, [X1],/NOERASE, thick=1, linestyle=2,color=106


CURSOR, X2, Y2, /DOWN, /DATA

if(!MOUSE.button EQ 4) then break

print,'Upper limit: ', X2
vline, [X2],/NOERASE, thick=1,linestyle=2,color=106

CURSOR, X3, Y3, /DOWN, /DATA
if(!MOUSE.button EQ 4) then break

endelse

endwhile
endif

;------------------------------------------------------------------------------------------------
;--------------------------------Equivalent width calculation------------------------------------

ew=0
ew_err_sqsum=0
ew_err=0
w_diff=wav[2]-wav[1]

w_X1=lamda[k]*(1+X1/c)
p1=min(abs(wav-w_X1), index1)

w_X2=lamda[k]*(1+X2/c)
p2=min(abs(wav-w_X2), index2)

if index2 lt index1 then begin
temp=index1
tempX=X1
index1=index2
X1=X2
index2=temp
X2=tempX
endif

print,'Continuum scaled by: ',cscale
print,'Spectrum shifted by: ',dv,' km/s'

for q=index1,index2 do begin

if flux[q] gt 0.0000 then begin
val=1-(flux[q]/cont[q]) 
endif else begin 
val=1 
endelse
w_diff=wav[q+1]-wav[q]
if q eq index2 then w_diff=wav[q]-wav[q-1]

ew=ew+val*w_diff 
tx=((err[q]^2)*(-1/(cont[q]))^2)						;CC Eq: 2.36 - Dimenssion corrected
ew_err_sqsum = ew_err_sqsum + tx*w_diff^2
endfor 
ew_err=sqrt(ew_err_sqsum)/(1+z)
ew=ew/(1+z)

;ERASE
;Cgdisplay,1000,800

for q=index1,index2 do begin
xpoly=[v_wav[q]-w_diff/2,v_wav[q]+w_diff/2, v_wav[q]-w_diff/2,v_wav[q]+w_diff/2]
ypoly=[norm_flux[q],norm_flux[q], 1, 1]
cgplots,xpoly,ypoly, color='dodger blue', thick=5
endfor

;If you get multiple plots one over the other comment the lines below
cgoplot,v_wav, norm_flux,xtitle=textoidl('!6 Velocity (km/s)'), ytitle=textoidl('!6 Normalized Flux'),$
xstyle=1,charsize=1.5,ycharsize=1.2,xcharsize=1.2,$
yrange=[0.0,1.5],xrange=[v_xmin,v_xmax],/ystyle,PSym=10,xthick=5,ythick=5,thick=2

cgoplot, v_wav, norm_cont,PSym=10,linestyle=2, color='green', thick=3
cgoplot, v_wav, norm_err,PSym=10, color='blue', thick=1
vline, [0.0],/NOERASE,linestyle=3, thick=2
cgtext, v_xmin+0.13*(v_xmax-v_xmin),0.2, textoidl(lname), charthick=2.0, color='red', charsize=2 
obs_wav_text = '!6 \lambda ='+strtrim(string(lamda[k]),2)+'A^{o}'
cgtext, v_xmax-0.35*(v_xmax-v_xmin),0.2, textoidl(obs_wav_text), charthick=2.0, color='green', charsize=2
;till here

vline, [X1],/NOERASE, thick=2, linestyle=2
vline, [X2],/NOERASE, thick=2,linestyle=2

ew1 = ew * 1e3
ew_err1 = ew_err * 1e3

;------------------------------------------------------------------------------------------------
;---------------------------------------Significance---------------------------------------------
;------------------------------------------------------------------------------------------------
sign=ew/ew_err
sign=float(round(sign*100)/100.0d)
significance=strmid(strtrim(string(long(sign)),2),0,5)+' sigma '
if sign ge 3 then adet[k] = 1 else adet[k] = 0									;Significance

output1 = 'Equivalent Width = '+strtrim(string(ew1),2)+' +/- '+strtrim(string(ew_err1),2)+'mA^{o}'
cgtext, v_xmin+0.01*(v_xmax-v_xmin),1.5+0.03, textoidl(output1), charthick=2.0, color=109, charsize=2.0

;------------------------------------------------------------------------------------------------
if total(det) eq 0 AND chkdf ne 1 then sm = 0

if KEYWORD_SET(dflag) then begin
if det[k] ne adet[k] and sm ne 0 then begin
print,''
print,'**Alert**'
print,'Mismatch noted between derived detection flag and input detection flag for ',strtrim(ion[k],2)
Print,'In certain cases of contamination, it is better to measure the upper limit using AOD method.'
print,'But be careful to modify the output for the ion and make it an upper limit'  
print,''
Print,'Please check the feature and enter the detection flag manually' 
print,'0 - Non-detection'
print,'1 - Detection'
print,''
if adet[k] eq 1 then print,'It could be a ',significance,'detection for the selected range'
if adet[k] eq 0 then print,'It could be a non-detection. Significance (< 3 sigma) = ',significance
print,''
detf=0
READ,detf,PROMPT='What will it be?: '
print,''
det[k]=detf
endif else begin
det[k]=adet[k]
sm=1
endelse
endif

if det[k] ne adet[k] AND sm ne 0 then begin
print,''
print,'**Is the detection flag for ',lname,'accurate?**'
if adet[k] eq 1 then print,'It could be a ',significance,'detection for the selected range'
if adet[k] eq 0 then print,'It could be a non-detection. Significance (< 3 sigma) = ',significance
endif

;------------------------------------------------------------------------------------------------
;--------------------------------------------AOD-------------------------------------------------

tauv=alog(cont/flux)
nav=3.768e14*tauv/(f[k]*lamdao[k]) 
nav_12=nav/1e12
print,''
detf=0
if sm eq 0 then begin
Print,'Please enter the detection flag manually' 
print,'0 - Non-detection'
print,'1 - Detection'
print,''
if adet[k] eq 1 then print,'It could be a ',significance,'detection for the selected range'
if adet[k] eq 0 then print,'It could be a non-detection. Significance (< 3 sigma) = ',significance
print,''
READ,detf,PROMPT='What will it be?: '
print,''
det[k]=detf
endif
;-----------------------------------------Detection----------------------------------------------
if(det[k] ne 0.00) then begin
na=0
v_diff=v_wav[2]-v_wav[1]
tau=0
sigma_tau_sumsq=0

for q=index1,index2 do begin
	if flux[q] gt 0 then begin
		v_diff=v_wav[q+1]-v_wav[q]
		if q eq index2 then v_diff = v_wav[q]-v_wav[q-1]
		na=na+(nav[q]*v_diff)
		tau=tau+(tauv[q]*v_diff)
		current = (err[q]^2)*((-1.0/flux[q])^2)*v_diff 		
		sigma_tau_sumsq= sigma_tau_sumsq + current                			
	endif 
endfor 

logna=alog10(na)
sigma_tau=sqrt(sigma_tau_sumsq)
sigma_na = na * (sigma_tau^2) / tau
sigma_logna=0.4343*sigma_na/na

print,significance,'detection: Proceeding to AOD'
Print,''
print,'*****************************'
print,lname
print,'*****************************'
print,'Lower Velocity Limit: ', strtrim(string(Long(X1)),2)
print,'Upper Velocity Limit:  ', strtrim(string(Long(X2)),2)
print,'Rest-frame Equivalent Width = ',strtrim(string(ew*1e3),2),' +/- ',strtrim(string(ew_err*1e3),2),' mA'
print,'log10(Integrated Apparent Column Density) = ',strtrim(string(logna),2),' +/- ',strtrim(string(sigma_logna),2)
print,'*****************************'
wl[k]=ew*1e3
wl_e[k]=ew_err*1e3

nf[k]=logna
nf_e[k]=sigma_logna 
range[k]='['+strtrim(string(Long(X1)),2)+','+strtrim(string(Long(X2)),2)+']'  ;Just trimming white spaces and decimals
v1[k]=long(X1)
v2[k]=long(X2)
endif else begin

;------------------------------------------------------------------------------------------------
;----------------------------------------Non-Detection-------------------------------------------

Print,''
print,'Non-detection: Proceeding to curve of growth'
print,'Significance: ',significance
Print,''
print,'Rest-frame Equivalent Width = ',strtrim(string(ew*1e3),2),' +/- ',strtrim(string(ew_err*1e3),2),' mA'
wlamda= 3*ew_err
print,'3 sigma upper limit: ',strtrim(string(wlamda*1e3),2),' A'
print,'Error underpredicted as continuum fitting uncertainty not considered'
print,'Oscillator strength: ',f[k]

nu=wlamda/((8.855e-21)*f(k)*((lamdao[k])^2))
nuf=alog10(nu)
Print,''
print,'*****************************'
print,lname
print,'*****************************'
print,'Lower Velocity Limit: ', strtrim(string(Long(X1)),2)
print,'Upper Velocity Limit:  ', strtrim(string(Long(X2)),2)
print,'Equivalent width upper limit < ',strtrim(string(wlamda*1e3),2),' mA'
print,'log10(Apparent Column Density limit) < ',strtrim(string(nuf),2)
print,'*****************************'

wl[k]=0.000
wl_e[k]= wlamda*1e3
nf[k]=0.000
nf_e[k]= nuf
range[k]='['+strtrim(string(Long(X1)),2)+','+strtrim(string(Long(X2)),2)+']'  ;Just trimming white spaces and decimals
v1[k]=long(X1)
v2[k]=long(X2)
endelse

endelse ;Main else

print,''
if outsidecoverage eq 1 then begin
k=k+1
if k ge j then break
print,'Proceeding to the next ion ',ion[k]
endif else begin
ni=ni+1
yn=''

if k eq j-2 then print,'*** Attention: Next ion is the last in the list ***'
if k eq j-1 then print,'*** Attention: Last ion in the list ***'

Print,''
print,'---------------------------------------------------'
if k ne j-1 then Print,'y  - to continue to the next ion ',ion[k+1] 
Print,'r  - to remeasure ',lname
Print,'c  - to scale the continuum and remeasure ',lname 
Print,'v  - to shift the spectrum and remeasure ',lname 
Print,'df - to flip the detection flag and remeasure ',lname 
Print,'k  - to measure a specific ion '
Print,'Enter any other key to exit '
print,'---------------------------------------------------'
Print,''
READ,yn,PROMPT='What will it be?: '
Print,''

chk=STRCMP(yn,'y',/FOLD_CASE)
chk2=STRCMP(yn,'r',/FOLD_CASE)
chkc=STRCMP(yn,'c',/FOLD_CASE)
chk4=STRCMP(yn,'k',/FOLD_CASE)
chkv=STRCMP(yn,'v',/FOLD_CASE)
chkdf=STRCMP(yn,'df',/FOLD_CASE)

if chk ne 1 AND chk2 ne 1 AND chkc ne 1 AND chk4 ne 1 AND chkv ne 1 AND chkdf ne 1 THEN begin
READ,yn,PROMPT='**ALERT** Do you really want to exit? (y/any key): ' 
chke=STRCMP(yn,'y',/FOLD_CASE)
if chke eq 1 THEN begin 
kn=k+1  ;Note about Kn - It is not used anywhere. It was intended to serve as the number of ions measured
break 
endif else begin
k=k+1 
cscale=0
endelse
endif 

if chk4 eq 1 THEN begin
kt=0
while (kt lt j) do begin
print,kt+1,' --- ',ion[kt]
kt=kt+1
endwhile
Print,''
READ,knew,PROMPT='Enter the index value of the ion to be measured: '
k=knew-1
if chkatm eq 1 then begin
X1d=vl[k]
X2d=vr[k]
endif
endif

if chk2 eq 1 or chkc eq 1 or chkv eq 1 THEN begin
k=k 
kn=k
endif

if chkdf eq 1 THEN begin
k=k 
kn=k
if det[k] eq 0 then det[k] = 1 else det[k] = 0
if sm eq 0 then sm=1
endif

if chk eq 1 THEN begin 
k=k+1 
kn=k 
if chkatm eq 1 then begin
X1d=vl[k]
X2d=vr[k]
endif
endif

endelse ;outside coverage check

endwhile ;Main While

;----------------------------------------End of loop---------------------------------------------

if ni ne 0 and dum ne 0 then DUM_ar=DUM_ar[1:*]

print,'---------------------------------------------------'
print,ni-dum,'  measurements made' ;'Number of iterations: ',ni
print,'---------------------------------------------------'
print,''
;------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------
;------------------------------------Writing to CSV file-----------------------------------------

yn=''
READ,yn,PROMPT='Do you wish to generate a .csv table? (y/any key): '
chk=STRCMP(yn,'y',/FOLD_CASE)
head=strarr(6)
head(0)='Ion'
head(1)='w (mA)'
head(2)='sigma_w (mA)'
head(3)='Na (cm^-2)'
head(4)='sigma_Na (cm^-2)'
head(5)='v (km/s)'
if chk eq 1 THEN begin
opfile=''
READ,opfile,PROMPT='Enter the output table filename: '
output=opfile+'.csv'
write_csv, output,ion,wl,wl_e,nf,nf_e,range,HEADER=head
print,'File saved: ',output
print,'---------------------------------------------------'
if dum ne 0 then begin
print,''
print,'It seems like your ionlist contains lines that are outside the wavelength coverage.'
print,ion[DUM_ar] 
print,''
print,'Please remove them from the table.'
endif
endif
;------------------------------------------------------------------------------------------------
;------------------------------------Writing LaTex output----------------------------------------
print,''
yn=''
READ,yn,PROMPT='Do you wish to generate a LaTex table? (y/any key): '
chk=STRCMP(yn,'y',/FOLD_CASE)
if chk eq 1 THEN begin
wl_r=strarr(file_lines(ionlist))
wl_e_r=strarr(file_lines(ionlist))
nf_r=strarr(file_lines(ionlist))
nf_e_r=strarr(file_lines(ionlist))
col1=strarr(file_lines(ionlist))

k=0
while (k lt j) do begin

if wl[k] eq 0 then begin 							;Non-Detection
wl_e_r[k]=strmid(strtrim(string(float(Round(wl_e[k]))),2),0,5)
nf_e_r[k]=strmid(strtrim(string(float(Round(nf_e[k]*100)/100.0d)),2),0,4)
wl_r[k]='$<$ '+ wl_e_r[k]
nf_r[k]='$<$ '+ nf_e_r[k]
endif else begin								;Detection
if wl_e[k] ge 10 then wl_e_r[k]=strmid(strtrim(string(Round(wl_e[k])),2),0,5)
;wl_r[k]=strmid(strtrim(string(float(Round(wl[k]*100)/100.0d)),2),0,5)
;wl_e_r[k]=strmid(strtrim(string(float(Round(wl_e[k]*100)/100.0d)),2),0,4)
wl_r[k]=strmid(strtrim(string(Round(wl[k])),2),0,5)
wl_e_r[k]=strmid(strtrim(string(Round(wl_e[k])),2),0,4)
nf_r[k]=strmid(strtrim(string(float(Round(nf[k]*100)/100.0d)),2),0,5)
nf_e_r[k]=strmid(strtrim(string(float(Round(nf_e[k]*100)/100.0d)),2),0,4)

if nf_e[k] gt 1.0 and det[k] eq 1 then begin				;Saturated
sat=''
READ,sat,PROMPT='Is '+ion[k]+' saturated? (y/n): '
chk=STRCMP(sat,'y',/FOLD_CASE)
if chk eq 1 then wl_r[k]= '$>$ '+ wl_r[k] & nf_r[k]= '$>$ '+ nf_r[k]
endif else begin

wl_r[k]= wl_r[k]+' $\pm$ '+wl_e_r[k]
nf_r[k]= nf_r[k]+' $\pm$ '+nf_e_r[k]
endelse

endelse
range[k]='$'+range[k]+'$'
col1[k]=ion[k]+' & '+wl_r[k]+' & '+nf_r[k]+' &   & '+range[k]+' \\'

k=k+1
endwhile

lopfile=''
READ,lopfile,PROMPT='Enter the LaTex output filename: '
loutput=lopfile+'.txt'
writecol,loutput,col1

print,'File saved: ',loutput
Print,'Values rounded to 2 decimal places. Please check the values and make neccesary approximations.'
Print,''
print,'Ion & Equivalent width (mA) & Column density (cm^-2) & Doppler parameter (blank) & Velocity (km/s)'
print,'---------------------------------------------------'
if dum ne 0 then begin
print,''
print,'I can see that your ionlist contains lines that are outside the wavelength coverage.'
print,ion[DUM_ar] 
print,''
print,'Please remove them from the table.'
endif
endif

;------------------------------------------------------------------------------------------------
;----------------------------Writing Input file for systemplot-----------------------------------
print,''
READ,yn,PROMPT='Do you wish to generate an input file? (y/any key): '
chk=STRCMP(yn,'y',/FOLD_CASE)
if chk eq 1 then begin

sipfile='atomsip'+strtrim(string(z),2)+'.dat'
writecol,sipfile,ion,lamdao,f,v1,v2,det

print,'File saved: ',sipfile
print,''
print,'Ion | Rest-wavelength | Oscillator strength | Velocity lower limit | Velocity upper limit'
print,'---------------------------------------------------'
if dum ne 0 then begin
print,''
print,'I can see that your ionlist contains lines that are outside the wavelength coverage.'
print,ion[DUM_ar] 
print,''
print,'Please remove them from the table.'
endif
endif
;------------------------------------------------------------------------------------------------

END
