;;BEGIN FUNCTION
;
;Name: InfoCalc
;
;purpose: calculates standard FT information from an image and a sigma
;
;inputs: image, The image array to be processed
;        sig,  The parcevals sigma calculated for the image
;
;keywords:
;
;output: returns a structure {Iks:Iks, Hks:Hks, KSI:KSI}
;
;calling example:
;                  result = InfoCalc(Image,Sigma)
;                  result = InfoCalc(Image, ParSig(Image))
;
;MODIFiCATIONS:
;   08-19-2009 WH Added SAVE_FILENAME keyword.
;                 Set this to a string contianing full path to a file to save
;                 the real and imaginary parts of fft for future use.
;   08-24-2009 WH Removed SAVE_FILENAME keyword and related code
;   08-24-2009 WH Added MBAND keyword.
;                 Set this to a structure to generate a band analysis of the current information.
;                  mband = {$
;                            nbands:n,$          ; This is the number of bands.  Also the number of bins in rev_ptr
;                           hist_ptr:ptr_new(hist),$  ;;this holds the histogram
;                           locs_ptr:ptr_new(locs),$ ;; this holds the locations of the bins
;                           rev_ptr:ptr_new(r),$   ;This holds the output of REVERSE_INDICES of HISTOGRAM
;                           result_ptr:ptr_new()}; this will hold the result of the analysis
function InfoCalc, Image, Sig,MBAND=mband
  divide = 100.0d
  binsize = Sig/divide
  s = binsize/2d
  Nsig = 10.0d
  
  F = FFT(Image,/DOUBLE)
  Fr = REAL_PART(F)
  Fi = IMAGINARY(F)
  
  ;clear memory of F
  F=0b
  
  ;Sets the dc components to 0.0
  Fr[0] = 0.0D
  Fi[0] = 0.0D
  
  ;produces the probabilities of Fr and Fi
  ;This is done by integrating the normal dist (mean=0,sig) about x over a window of +-binsize/2
  ;This is done by subtracting two cumulative integrals from each other using the gauss int function
  
  ;;check if sig is an array (if we are doing manual entry of sigmas, this will occur)
  ;;if it is then we must do things for real and imaginary separately
  if (n_elements(sig) gt 1) then begin
    pr = gaussint((Fr+s[0])/(sig[0]))-gaussint((Fr-s[0])/(sig[0]))
    pi = gaussint((Fi+s[1])/(sig[1]))-gaussint((Fi-s[1])/(sig[1]))
  endif else begin
    pr = gaussint((Fr+s)/(sig))-gaussint((Fr-s)/(sig))
    pi = gaussint((Fi+s)/(sig))-gaussint((Fi-s)/(sig))
  endelse
  
  ;;free Fr Fi memory
  Fr = 0b
  Fi = 0b
  
  ;Calculates information with units of bits
  ;by summing -log2(p) over the whole image
  Ir = -1D*ALOG(temporary(pr))/ALOG(2D)
  Ii = -1D*ALOG(temporary(pi))/ALOG(2D)
  
  ;;do marching band here if asked
  if keyword_set(mband) then begin
    ;;set up a results array for the banding
    bresults = dblarr((*mband).nbands)
    
    ;;loop through the bands
    for i = 0L, (*mband).nbands-1L do begin
      ;;if there are elements in the bin, then sum them
      IF ((*(*mband).rev_ptr)[i] NE (*(*mband).rev_ptr)[i+1]) THEN begin
        the_binR = Ir[(*(*mband).rev_ptr)[(*(*mband).rev_ptr)[i] : (*(*mband).rev_ptr)[i+1]-1]]
        the_binI = Ii[(*(*mband).rev_ptr)[(*(*mband).rev_ptr)[i] : (*(*mband).rev_ptr)[i+1]-1]]
        
        ;;the result is the average info in that band
        bresults[i] = (TOTAL(the_binR,/NAN)+TOTAL(the_binI,/NAN))/(1.0*(*(*mband).hist_ptr)[i])
      endif else begin ;; no elements in that band
        bresults[i] = 0D
      endelse
    endfor
    
    ;;free up results pointer
    if ptr_valid((*mband).result_ptr) then (*(*mband).result_ptr)=0b
    ptr_free,(*mband).result_ptr
    ;;load up the pointer for results
    (*mband).result_ptr = ptr_new(bresults)
  endif
  
  ;;total infos
  Irtot = total(temporary(Ir),/NAN)
  Iitot = total(temporary(Ii),/NAN)
  
  ;;total info
  Iks = Irtot + Iitot
  
  ;  ;; free Ir & Ii memory
  ;  Ir = 0b
  ;  Ii = 0b
  
  
  ;This section calculates the Entropy of the Image
  ;This is done by suming -plog(p) over the range -Nsig*sig to Nsig*sig
  ;this is done in discrete bin sized chunks
  
  ;;do it for real and imaginary if manual sigmas passed
  if (n_elements(sig) gt 1) then begin
    xr = (dindgen(2D*Nsig*divide+1D)-Nsig*divide)*binsize[0] ;Produces the range for calculating Hks
    phr = gaussint(shift(xr,-1)/sig[0])-gaussint(xr/sig[0])
    xi = (dindgen(2D*Nsig*divide+1D)-Nsig*divide*binsize[1]);Produces the range for calculating Hks
    phi = gaussint(shift(xi,-1)/sig[01])-gaussint(xi/sig[1])
    
    ;removes last element because it is not a valid probability
    ph=[ phr[0:2D*Nsig*divide-1D], phi[0:2D*Nsig*divide-1D] ]
    Hc_scale = 1D ;; inserted 07-15-11 to account for real and imaginary contributions to entropy
  endif else begin
    x = (dindgen(2D*Nsig*divide+1D)-Nsig*divide)*binsize ;Produces the range for calculating Hks
    
    ph = gaussint(shift(x,-1)/sig)-gaussint(x/sig)
    
    ;removes last element because it is not a valid probability
    ph=ph[0:2D*Nsig*divide-1D]
    Hc_scale = 2D ;; inserted 07-15-11 to account for real and imaginary contributions to entropy
  endelse
  
  
  ;Caculates Hks by multiplying H for an individual coef value: TOTAL(-ph*ALOG(ph)/ALOG(2), /NaN)
  ;by the number of elements in the FFT 2*N_ELEMENTS(image) (takes into account both real and imaginary parts)
  Hc = TOTAL(-ph*ALOG(ph)/ALOG(2D), /NaN)*Hc_scale ;;entropy per coefficient
   ;;Hc_scale mutiplier inserted 07-15-11 to account for real and imaginary contributions to entropy
   
  Hks = 2D*N_ELEMENTS(image)*Hc/Hc_scale ;; Hks already has the scaling, so we need to drop Hc_scale here.
  
  
  return, {Iks:Iks, Hks:Hks, KSI:Hks-Iks,Hc:Hc}
  
END