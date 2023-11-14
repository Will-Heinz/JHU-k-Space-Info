;;BEGIN FUNCTION
;
;Name: Info_surface
;
;purpose: takes an array and creates an new array with N+1 dimensions 
;          eg Im[x,y]=i -> Im_M[x,y,z]=1 for z=i and = 0 for z!=i 
;
;inputs: image, An image array to be processed, will be bytscaled 
;        
;
;keywords: Top, sets the size for the new dimension autoset to 256
;
;output: returns a N+1 dimensional surface of the image 
;
;calling example: 
;                  Image = Info_Surface(Image)
;                  
;



function info_surface,image,TOP=top
  COMPILE_OPT IDL2
  
  ;;size of the array
  n = n_elements(image)
  
  ;;defasult size of new dimension
  if not(keyword_set(top)) then top = 256
  
  ;;make sure the image is some kind of integer type
  t = size(image,/TYPE)
  if NOT( (t gt 0) AND ((t lt 4) OR (t gt 11) )) then begin
    image = bytscl(image,TOP=top)
  endif
  
  ;;make sure all values are positve
  if (total(image ge 0) ne n) then return,-1L ;;not all positive values
  
  ;;make an array to index pixels of value 1 in volume
  ones_index = lindgen(n) + n*image
  
  ;;image size
  s = size(image,/DIMENSIONS)
  
  ;;build it
  image3d = dblarr([s,top])
  image3d[ones_index] = 1D
  
  return, image3d
end

