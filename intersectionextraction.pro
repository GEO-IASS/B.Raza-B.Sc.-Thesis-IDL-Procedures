;+
; NAME:
;    INTERSECTIONEXTRACTION
;
; PURPOSE:
;    Active Regions defined by thresholding at a certain value are smaller than the real ones. 
;    To find the real active regions one has to match patches of same polarity with the are defined by thresholding.
;    This mathematical operation is equivalent to the intersection of two sets. 
;    The set of pixels S={s_i| where |y_ij| gt thr} is a superset of the set of patches with same polarity M:
;    S <= M. If S_i = M_i then M_i can be defines as an active region.
;    This function extracts the array indices of the pixel which are corresponding to active regions. 
;    These indices are then used to create a logical mask for active region pixels
;
; CATEGORY:
;    Image Processing
;
; CALLING SEQUENCE:
;    result = intersectionextraction(chart, thr)
;
; INPUTS:
;    CHART:  A synoptic chart magnetogram
;    THR:  This function uses the extrema of the array as reference points for thresolding
;    THR=0.5 correspons to threshold=0.5*max/min(chart). Other methods need changes in the code.
;
;
; OUTPUTS:
;    RESULT:  A logical array with zeroes and ones, which can be used to find active regions.
;
; USES:
;    AR/QS Separation
;

; EXAMPLE:
;    data = randomn(seed, 10000, 10000 )
;    thr = 0.5
;    result = intersection(data,thr)
;
; MODIFICATION HISTORY:
;    Written by:  Belal-Kevin Raza(belal-kevin.raza@stud.uni-goettingen.de); 19.12.2016
;-



FUNCTION INTERSECTIONEXTRACTION, chart,thr 


sz=size(chart)

;begin with defining the thresholds. Because there exists an asymmetry between positive and negative magnetic flux
;on my data, I will process them independently.

thrp = thr*max(chart)
thrm=thr*min(chart)

; In the next step we will label all independent regions with the LABEL_REGION function. 
;for further information check https://www.cis.rit.edu/class/simg782/lectures/lecture_01/lec782_05_01a.pdf
chartthrplus=chart
chartthrminus=chart


chartthrplus[where(chartthrplus lt thrp)] =0
chartthrminus[where(chartthrminus gt thrm)] =0
regions_m=LABEL_REGION(chartthrminus,/ALL_NEIGHBORS,/ULONG)
regions_p=LABEL_REGION(chartthrplus,/ALL_NEIGHBORS,/ULONG)
maxregp=max(regions_p)
maxregm=max(regions_m)

mtable=lindgen(maxregm)+1
mlist=list()
for j=1l,maxregm do begin
  mlist.add, where(regions_m eq j)
endfor
m=orderedhash(mtable, mlist)

plist=list()
ptable=ulindgen(maxregp)+1
for j=1l,maxregp do begin

  plist.add,where(regions_p eq j)
endfor

p=orderedhash(ptable, plist)

undefine, plist,mlist
undefine, chartthrminus,chartthrplus


;now do the same thing for the polarity patches
;polarity +-
;
plus=intarr(sz[1],sz[2])
plus[where(chart gt 0)] = 1




minus=intarr(sz[1],sz[2])
minus[where(chart lt 0)] = 1


;Use erode to remove all negligble regions
s=make_array(5,5,value=1)
plus=erode(plus,s)
plus=dilate(plus,s)
minus=erode(minus,s)
minus=dilate(minus,s)


regions_pp=LABEL_REGION(plus,/ALL_NEIGHBORS,/ULONG)
regions_mm=LABEL_REGION(minus,/ALL_NEIGHBORS,/ULONG)

maxpp=max(regions_pp)
maxmm=max(regions_mm)



mmlist=list()
mmtable=lindgen(maxmm)+1
for j=1l,maxmm do begin

  mmlist.add, where(regions_mm eq j)
endfor
mm=orderedhash(mmtable,mmlist)




pplist=list()
pptable=lindgen(maxpp)+1
for j=1l,maxpp do begin
  pplist.add, where(regions_pp eq j)
endfor
pp=orderedhash(pptable,pplist)


undefine, pplist,mmlist

;now do the set intersection. 2 for + 1 if is not very elegant but gets the job done
;Maybe someone else can make use of IDL's matrix notation to decrease computation time.


maskadd=list()

for j=1l,maxregp do begin
  for i=1l,maxpp do begin
    if (SetIntersection(p[j],pp[i]) ne [-1]) then begin
      maskadd.add,pp[i],/extract
    endif
  endfor
endfor
undefine, pp,p

for j=1l,maxregm do begin
  for i=1l,maxmm do begin
    if (SetIntersection(m[j],mm[i]) ne [-1]) then begin
      maskadd.add,mm[i],/extract
    endif
  endfor
endfor
undefine, mm,m
;build the mask
maske=dblarr(sz[1],sz[2])
maskindices=maskadd.toarray(type='ulong')
maske[maskindices]=1
undefine, maskindices,maskadd

;DONE
return, maske
end

