//kp_flopscalc.h

#ifndef __SEGER__KP_FLOPSCALC_H__
#define __SEGER__KP_FLOPSCALC_H___

class kp_fpc_data
{
 public:
   enum accumulators
     {
	A_GEMM = 0, A_GEMV, A_EIG, A_INVERSE, A_SPARCE, A_OTHER, A_SIZE
     };
   kp_fpc_data() {reset();};
   void add_gemm(long n)    {check(); a[A_GEMM]    += n;};
   void add_gemv(long n)    {check(); a[A_GEMV]    += n;};
   void add_eig(long n)     {check(); a[A_EIG]     += n;};
   void add_inverse(long n) {check(); a[A_INVERSE] += n;};
   void add_sparce(long n)  {check(); a[A_SPARCE]  += n;};
   void add_other(long n)   {check(); a[A_OTHER]   += n;};
   
   long get_gemm()    {return a[A_GEMM];};
   long get_gemv()    {return a[A_GEMV];};
   long get_eig()     {return a[A_EIG];};
   long get_inverse() {return a[A_INVERSE];};
   long get_sparce()  {return a[A_SPARCE];};
   long get_other()   {return a[A_OTHER];};
   
   long sum();
   void reset();
   void check();
   void operator+=(const kp_fpc_data&);
 private:
   long a[A_SIZE];
};

extern kp_fpc_data kp_fpc_global;

#endif
