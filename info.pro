;BEGIN wrapper functions for InfoLib

@info_sigma
@info_sigma_multi
@info_surface
@infocalc
@info_sigma_fish_multi


;;BEGIN FUNCTION
;
;Name: info_multi
;
;purpose: calculates FT information for the N+1 version of the passed image
;
;inputs: Image, The image array used to calculate sigma
;
;keywords:
;
;output: returns a structure {Iks:Iks, Hks:Hks, KSI:KSI}
;
;
;calling example:
;                  result = Info_multi(Image)
;
;

function info_multi, Image,ORDERED=ordered,SIG=sig,_EXTRA=xtra,MANSIG=mansig,NOMULTI=nomulti

  ;;use manual sigma if asked added 02-26-2010
  if keyword_set(mansig) then begin
    sig = mansig
  endif else begin
    sig = info_sigma_multi(image)
  endelse
  
  ;;removed 08-28-12
  ;img = info_surface(image)
  ;
  ;;skip the multi if asked -- added 08-28-12
  if keyword_Set(nomulti) then begin
    img = image
  endif else begin
    img = info_surface(image)
  endelse
  
  
  if keyword_set(ordered) then img=info_order(img)
  
  return, infocalc(Img, sig,_EXTRA=xtra)
  
end

;;this is a version for the zebrafish data
function info_fish_multi,embryo_dex,nuclei,dims,ORDERED=ordered,SIG=sig,_EXTRA=xtra,MANSIG=mansig

  if keyword_set(mansig) then begin
    sig = mansig
  endif else begin
    sig = info_sigma_fish_multi(nuclei,dims)
  endelse
  
  img = dblarr(dims)
  img[embryo_dex] = 1.0
  
  if keyword_set(ordered) then begin
    img=info_order(img)
    
  endif
  
  return, infocalc(Img, sig,_EXTRA=xtra)
  
end

;;BEGIN FUNCTION
;
;Name: info
;
;purpose: calculates FT information of the passed image
;
;inputs: Image, The image array used to calculate sigma
;
;keywords:
;
;output: returns a structure {Iks:Iks, Hks:Hks, KSI:KSI}
;
;
;calling example:
;                  result = Info_multi(Image)
;
;

function info, image,ORDERED=ordered,SIG=sig,_EXTRA=xtra

  sig = info_sigma(image)
  
  if keyword_set(ordered) then return,info_calc(info_order(image), sig)
  
  return, infocalc(image, sig,_EXTRA=xtra)
  
end
