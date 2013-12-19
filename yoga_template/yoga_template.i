plug_in,"yoga_template";

require,"yoga.i";

extern yoga_immult;
/* DOCUMENT yoga_immult
   yoga_immult,dest,src

   This routine multiplies image dest by src

   can only be called as a subroutine

   SEE ALSO: 
 */

extern yoga_magma_evd;
/* DOCUMENT yoga_magma_evd
   ev = yoga_magma_evd( jobz, uplo, h_A, d_R)

   input:
     - jobz: 'N' ou 'V'
     - uplo: 'U' ou 'L'
     - h_A : host matrix to decompose
     - d_R : KARMA obj containing eingen vectors

   return:
     - ev : eigen value
 
   SEE ALSO: 
 */
