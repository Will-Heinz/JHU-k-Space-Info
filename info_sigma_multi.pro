;;BEGIN FUNCTION
;
;Name: info_sigma_multi
;
;purpose: quick calculates a parsevals sigma for the N+1 version of the passed image
;
;inputs: data, The image array used to calculate sigma
;        
;keywords: 
;
;output: returns the sigma of the normal distribution of the FT coeffs from parceval's theorm
;        for a image array that has been processed by info_surface
;
;calling example: 
;                  sigma = Info_sigma_multi(Image)
;                  
;



function info_sigma_multi,data
  COMPILE_OPT IDL2
  
  ;;set the maximum value of the bytscl version of data
  top = 256L
  
  ;data = bytscl(data,TOP=top)
  
  ;;get some statistics of the data
  ;;size of data
  n=double(n_elements(data))
  
  ;;total size of multi-dimensional volume
  n_multi = n * double(top)
  
  ;;size of data array
  sz = double(size(data,/DIMENSIONS))
  
  ;;find the sigma
    ;;first find the mean.  The mean of a multi-d is the number of elements
    ;  divided by the size of the multiD array
    dc = n/(product(sz)*double(top))
    
    ;;get the sigma without the dc component
    sigma = sqrt( (n*abs(1D - dc)^2.0D + (n_multi-n)*(-dc)^2.0D))/n_multi
    
  ;;split it for real and imaginary parts
  
    result = sigma/sqrt(2.0D)
  
  ;;return result
  return,result
end