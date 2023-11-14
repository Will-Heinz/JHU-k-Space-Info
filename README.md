# JHU-k-Space-Info
IDL routines to calculate spatial information content of images using k-space information.

The k-space information (kSI) quantifies the spatial information of an image (or other data of arbitrary dimensions) using an approach based on Fourier coefficient distributions.

Details of the method can be found here: 
Heinz, W.F., Werbin, J.L., Lattman, E. et al. Computing Spatial Information from Fourier Coefficient Distributions. J Membrane Biol 241, 59–68 (2011). https://doi.org/10.1007/s00232-011-9362-x

## Basic usage
Info.pro is a wrapper function that calculates the kSI, Hks, and Iks. These terms are defined in the J.Membrane Biol. article.

Call info.pro as result = info(image). 

Input: image - an image (more generally, an array of up to 8 dimensions).

Output: result - a structure {Iks:Iks, Hks:Hks, KSI:KSI}.

See the code for details and keyword options.

## Disclaimer
These routines were written between 2009 and 2011, using IDL v 8.1 an an Intel-based Apple Macintosh computer.  Much has changed since then, and the routines have not been updated. Due to changes in computer architectures and IDL itself in the last 10+ years, I can not guarantee that they will work in other versions/installations of IDL. 

