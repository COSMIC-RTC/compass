plug_in,"yoga_template";

require,"yoga.i";

extern yoga_magma_syevd;
/* DOCUMENT yoga_magma_syevd
   ev = yoga_magma_syevd( jobz, uplo, h_A, d_R)

   input:
     - jobz: 'N' ou 'V'
     - uplo: 'U' ou 'L'
     - h_A : host matrix to decompose
     - d_R : CArMA obj containing eigen vectors

   return:
     - ev : eigen value
 
   SEE ALSO: 
 */

extern yoga_magma_getri;
/* DOCUMENT yoga_magma_getri
   yoga_magma_getri, h_A, d_R
   
   equivalent of LUsolve in Yorick : compute the invert of matrix h_A
   
   input:
     - h_A : host matrix to invert
     - d_R : CArMA obj containing result

   SEE ALSO: 
 */
