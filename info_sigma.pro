;;BEGIN FUNCTION
;
;Name: info_sigma
;
;purpose: quick calculates a parsevals sigma for the passed image
;
;inputs: data, The image array used to calculate sigma
;        
;keywords: 
;
;output: returns the sigma of the normal distribution of the FT coeffs from parceval's theorm
;
;calling example: 
;                  sigma = Info_sigma(Image)
;                  
;


function info_sigma,data
  COMPILE_OPT IDL2
  
  ;size of the array
   n=n_elements(data)
   
   ;;find the mean
    dc = mean(data,/NAN,/DOUBLE)
    
    ;;get the sigma without the dc component
    result = sqrt(total(abs(data-dc)^2.0D,/NAN,/DOUBLE))/n
    
    ;;split it
    result/=sqrt(2.0D)
  
    return,result
 

end
