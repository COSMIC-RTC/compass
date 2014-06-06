//kp_kalman_core_full_GPU.cpp

#include "kp_kalman_core_full_GPU.h"
#include <fstream>
#include <iomanip>

#define __SP setprecision(20)<<


kp_kalman_core_full_GPU::kp_kalman_core_full_GPU(const kp_matrix& D_Mo_,
	const kp_matrix& N_Act_, const kp_matrix& PROJ_, bool isZonal_, cublasHandle_t cublasHandle_) : kp_kalman_core_full(D_Mo_, N_Act_, PROJ_, isZonal_), cublasHandle(cublasHandle_)
{

	cu_D_Mo = D_Mo_;
	cu_H_inf.resize(0,0);
	
	cu_N_Act = N_Act_;
	cu_PROJ = PROJ_;
	cu_U_km2.resize(nb_act);
	cu_U_km2.zeros();
	cu_U_km1.resize(nb_act);
	cu_U_km1.zeros();
	cu_X_kskm1.resize(nb_n);
	cu_X_kskm1.zeros();
		
	//variables temporaires de next_step
	cu_Y_k.resize(nb_p);
	cu_Nact_Ukm2.resize(nb_az);
	cu_tmp_vec1.resize(nb_az);
	cu_innovation.resize(nb_p);
	cu_X_kskm1_tmp.resize(nb_n);	
	cu_X_kp1sk.resize(nb_n);
	cu_X_kp1sk_tmp.resize(nb_az);
	cu_Y_kskm1.resize(nb_p);
	cu_A1_00_Xkdebut.resize(nb_az);
	cu_A1_01_Xkfin.resize(nb_az);
	cu_X_kp1sk_debut.resize(nb_az);
	cu_U_k.resize(nb_act);
		
}



kp_kalman_core_full_GPU::~kp_kalman_core_full_GPU()
{
}

void kp_kalman_core_full_GPU::calculate_gain(real bruit_pix,
		real k_W, 
		const  kp_matrix&  SigmaV,
		const  kp_vector& atur_,
		const  kp_vector& btur_)
{
	//atur = atur_;
	//btur = btur_;
	cu_atur = atur_;
	cu_btur = btur_;

	
	

	if (cu_atur.size() != nb_az)
	{
		cerr<<"Error | kp_kalman_core_full_GPU::kp_kalman_core_full_GPU | size problem atur"<<endl;
		exit(EXIT_FAILURE);
	}
	if (cu_btur.size() != nb_az)
	{
		cerr<<"Error | kp_kalman_core_full_GPU::kp_kalman_core_full_GPU | size problem btur"<<endl;
		exit(EXIT_FAILURE);
	}
	int i;

	bool AR1 = true;
	for(i = 0 ; i < btur_.size() ; i++) AR1 &= (btur_.d[i]==0);
	if (AR1) ordreAR = 1 ; else ordreAR = 2;

int expression = 2;
	real seuil = 1/pow(10,14);
	int boucle = 0;
	real ecart = 1.0;
	const int boucle_max = 50;
	const real SigmaW = k_W*bruit_pix;
	//real Trac_T[boucle_max];
	real trac1_tmp, trac2_tmp;

//kp_cu_timer temps_alphak, temps_betak, temps_Tk, temps_inversion ;

	cu_H_inf.resize(0,0);


	vector<kp_cu_matrix*> cu_ms ;

	kp_cu_matrix cu_alpha_kp1(0,0);
	if (ordreAR == 1) cu_alpha_kp1.resize(nb_az,nb_az);
	else cu_alpha_kp1.resize(nb_n,nb_n);
	cu_alpha_kp1.zeros();
	
	// pour ordreAR 1 et 2
	kernel_memcpy_diag(cu_alpha_kp1.d_cu, cu_atur.d_cu, 0, 0, nb_az, cu_alpha_kp1.dim1);
	//pour ordreAR 2 uniquement
	if (ordreAR == 2)
	{
		kernel_memset_diag(cu_alpha_kp1.d_cu, 1.0, 0, nb_az, nb_az, cu_alpha_kp1.dim1);
		kernel_memcpy_diag(cu_alpha_kp1.d_cu, cu_btur.d_cu, nb_az, 0, nb_az, cu_alpha_kp1.dim1);
	}
	kp_cu_matrix cu_beta_kp1(0,0);
	if (ordreAR == 1) cu_beta_kp1.resize(nb_az, nb_az);
	else cu_beta_kp1.resize(nb_n, nb_n);
	cu_beta_kp1.zeros();



	// zeros_Dmo = [0 D_Mo]
	kp_cu_matrix cu_zeros_Dmo(0,0);
	kp_cu_matrix cu_zeros_nbp_nbaz(nb_p,nb_az);
	cu_zeros_nbp_nbaz.zeros();
	cu_ms.push_back(&cu_zeros_nbp_nbaz);
	cu_ms.push_back(&cu_D_Mo);
	kp_cu_horizcat(cu_zeros_Dmo,cu_ms);
	cu_ms.clear();
	cu_zeros_nbp_nbaz.resize(0,0) ; 

	kp_cu_matrix cu_C1;	
	// si ordreAR 1 C1 = D_Mo
	if (ordreAR == 1) cu_C1 = cu_D_Mo;
	//si ordreAR 2 C1 = [zeros(nb_p,nb_az) D_Mo]
	else cu_C1 = cu_zeros_Dmo;

	if (expression==2) cu_zeros_Dmo.resize(0,0);


	// VRAI SEULEMENT SI SigmaW = k*Id (k reel)	
	// beta_kp1 = (C1)T * C1
	kp_cu_gemm (cublasHandle,'T', 'N', 1, cu_C1, cu_C1, 0, cu_beta_kp1);

	// beta_kp1 = (C1)T * C1 / SigmaW (avec SigmaW reel)
	cu_beta_kp1 *= 1/SigmaW ;

		



	kp_cu_matrix cu_SigmaV(SigmaV); 
	kp_cu_matrix cu_T_kp1sk(0,0);
	if (ordreAR == 1) 
	{
		cu_T_kp1sk = cu_SigmaV;
	}
	else
	{
		// beta_kp1_tmp = [SigmaV ; zeros(nb_az, nb_az)];
		kp_cu_matrix cu_zeros_nbaz_nbaz(nb_az, nb_az);
		cu_zeros_nbaz_nbaz.zeros();
		kp_cu_matrix cu_T_kp1sk_tmp(cu_SigmaV);
		cu_ms.push_back(&cu_T_kp1sk_tmp); cu_ms.push_back(&cu_zeros_nbaz_nbaz);
		kp_cu_matrix cu_beta_kp1_tmp(nb_n, nb_az);
		kp_cu_vertcat(cu_beta_kp1_tmp, cu_ms);
		cu_zeros_nbaz_nbaz.resize(0,0) ;
		cu_T_kp1sk_tmp.resize(0,0);
		cu_ms.clear();

		// T_kp1sk = [beta_kp1_tmp zeros_nbn_nbaz];
		kp_cu_matrix cu_zeros_nbn_nbaz(nb_n, nb_az);
		cu_zeros_nbn_nbaz.zeros();
		cu_ms.push_back(&cu_beta_kp1_tmp); cu_ms.push_back(&cu_zeros_nbn_nbaz);
		kp_cu_horizcat(cu_T_kp1sk, cu_ms);
		cu_beta_kp1_tmp.resize(0,0);
		cu_zeros_nbn_nbaz.resize(0,0);
		cu_ms.clear();
	}

	if(ordreAR==2 || expression==2) cu_SigmaV.resize(0,0) ; 



	kp_cu_matrix cu_alpha_k(cu_alpha_kp1.dim1, cu_alpha_kp1.dim2);
	kp_cu_matrix cu_T_k(cu_T_kp1sk.dim1, cu_T_kp1sk.dim2);
	kp_cu_matrix cu_beta_k(cu_beta_kp1.dim1, cu_beta_kp1.dim2);
	
	kp_cu_matrix cu_IBG_1(cu_alpha_k.dim2, cu_beta_k.dim1);



	kp_cu_matrix cu_Tk_IBG1(cu_T_k.dim1, cu_IBG_1.dim2);
	kp_cu_matrix cu_betak_Tk(cu_beta_k.dim1, cu_T_k.dim2);
	kp_cu_matrix cu_alphak_IBG1(cu_alpha_k.dim1, cu_IBG_1.dim2);
	kp_cu_matrix cu_alphak_IBG1_betak(cu_alpha_k.dim1, cu_beta_k.dim2);
	kp_cu_matrix cu_Tk_IBG1_alphak(cu_T_k.dim1, cu_alpha_k.dim2);

	kp_cu_vector cu_diag_cu_Tkp1sk(cu_T_kp1sk.dim1);
	kp_cu_vector cu_diag_cu_Tk(cu_T_k.dim1);


	ofstream fichier;



	while( (ecart>seuil) && (boucle < boucle_max) )
	{
		// mise a jour
		cu_alpha_k = cu_alpha_kp1;
		cu_beta_k = cu_beta_kp1;
		cu_T_k = cu_T_kp1sk;


		// Calcul de IBG_1
		// betak_Tk = beta_k * T_k
/*cout<<"cu_beta_k="<<cu_beta_k.dim1<<"x"<<cu_beta_k.dim2<<endl;
cout<<"cu_T_k="<<cu_T_k.dim1<<"x"<<cu_T_k.dim2<<endl;
cout<<"cu_betak_Tk="<<cu_betak_Tk.dim1<<"x"<<cu_betak_Tk.dim2<<endl;
kp_matrix betak_test, Tk_test, betak_Tk_test;
kp_cu2kp_matrix(betak_test, cu_beta_k);
kp_cu2kp_matrix(Tk_test, cu_T_k);
kp_cu2kp_matrix(betak_Tk_test, cu_betak_Tk);
fichier.open("betak_test.dat",ios::out);
   for (int i=0 ; i<betak_test.dim1 ; i++)
   {
       for (int j=0 ; j<betak_test.dim2 ; j++)
       {
          fichier << betak_test(i,j)<< " ";
       }
      fichier <<endl;
   }
fichier.close();
fichier.open("Tk_test.dat",ios::out);
   for (int i=0 ; i<Tk_test.dim1 ; i++)
   {
      for (int j=0 ; j<Tk_test.dim2 ; j++)
      {
         fichier << Tk_test(i,j)<< " ";
      }
      fichier <<endl;
   }
fichier.close();
fichier.open("betak_Tk_test.dat",ios::out);
   for (int i=0 ; i<betak_Tk_test.dim1 ; i++)
   {
      for (int j=0 ; j<betak_Tk_test.dim2 ; j++)
      {
         fichier << betak_Tk_test(i,j)<< " ";
      }
      fichier <<endl;
   }
fichier.close();*/
		kp_cu_gemm (cublasHandle, 'N', 'N', 1 , cu_beta_k, cu_T_k, 0, cu_betak_Tk);

		/*if (ordreAR ==1)
			
			for (i=0 ; i<nb_az ; i++) (cu_betak_Tk.d_cu[i * cu_betak_Tk.dim1 + i]) += 1;
		else
		{
			for (i=0 ; i<nb_n ; i++) (cu_betak_Tk.d_cu[i * cu_betak_Tk.dim1 + i]) += 1;
		}*/
		kernel_add_diag_const(cu_betak_Tk.d_cu, 1, cu_betak_Tk.dim1);



		cu_IBG_1 = cu_betak_Tk;
//temps_inversion.start();
		cu_IBG_1.inverse();
//temps_inversion.pause();

		


		

	
//temps_alphak.start();
		// Calcul de alpha_kp1
		// alphak_IBG1 = alpha_k * IBG_1
		kp_cu_gemm (cublasHandle, 'N', 'N', 1 , cu_alpha_k, cu_IBG_1 , 0, cu_alphak_IBG1);

        	// alpha_kp1 = alphak_IBG1 * alphak (= alpha_k * IBG1 * alpha_k) 
		kp_cu_gemm (cublasHandle, 'N', 'N', 1 , cu_alphak_IBG1, cu_alpha_k , 0, cu_alpha_kp1);

//temps_alphak.pause();




//temps_betak.start();
		// Calcul de beta_kp1

		// alphak_IBG1_betak = alphak_IBG1 * beta_k (= alpha_k * IBG1 * beta_k)
		kp_cu_gemm (cublasHandle, 'N', 'N', 1 , cu_alphak_IBG1, cu_beta_k , 0, cu_alphak_IBG1_betak);


		//beta_kp1 = alphak_IBG1_betak * (alpha_k)T (= alpha_k * IBG1 * beta_k * (alpkha_k)T)
		kp_cu_gemm (cublasHandle, 'N', 'T', 1 , cu_alphak_IBG1_betak, cu_alpha_k , 1, cu_beta_kp1);
//temps_betak.pause();




//temps_Tk.start();
		//Calcul de T_kp1sk

		// Tk_IBG1 = T_k * IBG1
		kp_cu_gemm (cublasHandle, 'N', 'N', 1, cu_T_k ,cu_IBG_1, 0, cu_Tk_IBG1);
		// Tk_IBG1_alphak = Tk_IBG1 * alpha_k (= T_k * IBG1 * alpha_k)
		kp_cu_gemm (cublasHandle, 'N', 'N', 1 , cu_Tk_IBG1 , cu_alpha_k , 0, cu_Tk_IBG1_alphak);		

		// T_kp1sk = (alpha_k)T * Tk_IBG1_alphak (= (alpha_k)T * T_k * IBG1 * alpha_k)
		kp_cu_gemm (cublasHandle, 'T', 'N', 1 , cu_alpha_k, cu_Tk_IBG1_alphak , 1, cu_T_kp1sk);
		
//temps_Tk.pause();
	



		kernel_get_diag(cu_diag_cu_Tkp1sk.d_cu, cu_T_kp1sk.d_cu, cu_T_kp1sk.dim1);	
		kernel_get_diag(cu_diag_cu_Tk.d_cu, cu_T_k.d_cu, cu_T_k.dim1);	

	
		trac1_tmp = kp_cu_reduce(cu_diag_cu_Tkp1sk);
		trac2_tmp = kp_cu_reduce(cu_diag_cu_Tk);
		
		//Trac_T[boucle]=trac1_tmp;
		ecart =fabs(trac1_tmp/trac2_tmp-1.0);
		boucle ++;


	}


/*cout<<"nb_boucle = "<<boucle<<endl;
cout<< "temps inversion = "<<temps_inversion.rez()<<endl;
cout<< "temps alpha_k = "<<temps_alphak.rez()<<endl;
cout<< "temps beta_k = "<<temps_betak.rez()<<endl;
cout<< "temps T_k = "<<temps_Tk.rez()<<endl;*/



	cu_T_k.resize(0,0);
	cu_alpha_kp1.resize(0,0) ;
	cu_alpha_k.resize(0,0);
	cu_beta_k.resize(0,0);
	cu_beta_kp1.resize(0,0) ;
	cu_IBG_1.resize(0,0);
	cu_Tk_IBG1.resize(0,0);
	cu_Tk_IBG1_alphak.resize(0,0);
	cu_alphak_IBG1_betak.resize(0,0);
	cu_alphak_IBG1.resize(0,0);
	cu_betak_Tk.resize(0,0);





	if ( (expression == 1) && (ordreAR == 1))
	{
		//calcul de Sinf_0_0 (matrice superieure gauche de S_inf)
        	kp_cu_matrix cu_Atur_Tkp1sk(nb_az,nb_az);
		kp_cu_matrix cu_Sinf_0_0(nb_az,nb_az);
		kernel_diag_mult3(cu_Atur_Tkp1sk.d_cu, cu_Sinf_0_0.d_cu, cu_T_kp1sk.d_cu, cu_atur.d_cu, cu_atur.size(), cu_T_kp1sk.dim1*cu_T_kp1sk.dim2);


		// Sinf_0_0 = Sinf_0_0 + SigmaV
		cu_Sinf_0_0 += cu_SigmaV ;
		cu_SigmaV.resize(0,0) ; 
		
		
		// calcul de Sinf_0 (matrice superieure de S_inf)
		kp_cu_matrix cu_Sinf_0(nb_az,nb_n);
		cu_ms.push_back(&cu_Sinf_0_0) ; cu_ms.push_back(&cu_Atur_Tkp1sk) ;
		// Sinf = [ Sinf_0_0  A_tur*T_kp1sk ] 
		kp_cu_horizcat(cu_Sinf_0, cu_ms);
		cu_ms.clear();
		cu_Sinf_0_0.resize(0,0);
		cu_Atur_Tkp1sk.resize(0,0);


		// calcul de Sinf_1 (matrice inferieure de S_inf)
		kp_cu_matrix cu_Sinf_1_0(nb_az,nb_az);
		kernel_diag_mult(cu_Sinf_1_0.d_cu, cu_T_kp1sk.d_cu, cu_atur.d_cu, cu_atur.size(), cu_T_kp1sk.dim1*cu_T_kp1sk.dim2);


		kp_cu_matrix cu_Sinf_1(nb_az,nb_n);
		cu_ms.push_back(&cu_Sinf_1_0) ; cu_ms.push_back(&cu_T_kp1sk) ;
		// Sinf_1 = [ Sinf_1_0  T_kp1sk] 
		kp_cu_horizcat(cu_Sinf_1, cu_ms);
		cu_ms.clear();
		cu_Sinf_1_0.resize(0,0);
		cu_T_kp1sk.resize(0,0) ; 
	
		// S_inf = [ Sinf_0 ; Sinf_1 ]
		kp_cu_matrix cu_S_inf(nb_n, nb_n);
		cu_ms.push_back(&cu_Sinf_0) ; cu_ms.push_back(&cu_Sinf_1) ;	 
		kp_cu_vertcat(cu_S_inf, cu_ms);	
		cu_ms.clear();
		cu_Sinf_0.resize(0,0);
		cu_Sinf_1.resize(0,0);
	






		kp_cu_matrix cu_Sinf_zerosDmot(nb_n,nb_p);

		// Sinf_zerosDmot = S_inf * [0 D_Mo]T 
		kp_cu_gemm (cublasHandle, 'N', 'T', 1 , cu_S_inf, cu_zeros_Dmo , 0, cu_Sinf_zerosDmot);

		cu_S_inf.resize(0,0);

	/*kp_matrix Sinf_zerosDmot;
	kp_cu2kp_matrix(Sinf_zerosDmot, cu_Sinf_zerosDmot);
	fichier.open("Sinf_zerosDmot.dat",ios::out);
	for(i=0;i<Sinf_zerosDmot.dim1;i++)
	{
		for (j=0;j<Sinf_zerosDmot.dim2;j++)
		{
			fichier<< __SP Sinf_zerosDmot(i,j)<<" ";
		}
		fichier << endl;
	}	
	fichier.close();*/






		// calcul de Sigma_tot = [0 D_Mo] * S_inf * [0 D_Mo]' + SigmaW
		kp_cu_matrix cu_inv_Sigmatot(nb_p,nb_p);
		// Sigma_tot = zeros_Dmo * Sinf_zerosDmot  ( = [0 D_Mo] * Sinf * [0 ; (D_Mo)T] )
		kp_cu_gemm (cublasHandle, 'N', 'N', 1 , cu_zeros_Dmo, cu_Sinf_zerosDmot , 0, cu_inv_Sigmatot);

		// Sigma_tot = Sigma_tot + SigmaW*Id (avec SigmaW reel)  (= [0 D_Mo] * Sinf * [0 ; (D_Mo)T] + SigmaW*Id)
		//for(i = 0 ; i < cu_inv_Sigmatot.dim1 ; i++) cu_inv_Sigmatot.d_cu[i * cu_inv_Sigmatot.dim1 + i] += SigmaW;
		
		
		kernel_add_diag_const(cu_inv_Sigmatot.d_cu, SigmaW, cu_inv_Sigmatot.dim1);
		cu_zeros_Dmo.resize(0,0) ; 
		
		//inversion de Sigma_tot
		// inv_Sigmatot = inv(Sigma_tot)
		cu_inv_Sigmatot.inverse();
		// ATTENTION !!!! Le pointeur Sigma_tot n'existe plus
		// H_inf = Sinf_zerosDmot * inv_Sigmatot  ( = Sinf * [0 (D_Mo)T] * inv(Sigma_tot))
		cu_H_inf.resize(nb_n,nb_p);
		kp_cu_gemm (cublasHandle, 'N', 'N', 1 , cu_Sinf_zerosDmot, cu_inv_Sigmatot , 0, cu_H_inf);
		cu_Sinf_zerosDmot.resize(0,0);
		cu_inv_Sigmatot.resize(0,0);
		
		
	}
	else
	{
		//Calcul de T_kp1sk * (C1)T

        	kp_cu_matrix cu_Tkp1sk_C1t(0,0);
		if (ordreAR==1) cu_Tkp1sk_C1t.resize(nb_az,nb_p);
		else cu_Tkp1sk_C1t.resize(nb_n,nb_p);
        	kp_cu_gemm (cublasHandle, 'N', 'T', 1 , cu_T_kp1sk, cu_C1 , 0, cu_Tkp1sk_C1t);

		cu_T_kp1sk.resize(0,0);

		
		//Calcul de C1 * T_kp1sk * C1' + SigmaW
        	kp_cu_matrix cu_inv_Sigmatot(cu_C1.dim1, cu_Tkp1sk_C1t.dim2);
		// Sigma_tot = SigmaW * Id (avec SigmaW reel)
		cu_inv_Sigmatot.zeros();
		//for(i = 0 ; i < cu_inv_Sigmatot.dim1 ; i++) cu_inv_Sigmatot.d_cu[i * cu_inv_Sigmatot.dim1 + i] += SigmaW;
		kernel_add_diag_const(cu_inv_Sigmatot.d_cu, SigmaW, cu_inv_Sigmatot.dim1);

	
		// Sigma_tot = Sigma_tot + C1 * Tkp1sk_C1t  ( = (C1 * T_kp1sk * (C1)T) + SigmaW*Id )
		kp_cu_gemm (cublasHandle, 'N', 'N', 1 , cu_C1, cu_Tkp1sk_C1t , 1, cu_inv_Sigmatot);

		// inversion de (Sigmatot)
		//inv_Sigmatot = inv(inv_Sigmatot)
        	cu_inv_Sigmatot.inverse();
		// ATTENTION !!!! Le pointeur Sigma_tot n'existe plus
        
		
		if (ordreAR == 1)
		{
			kp_cu_matrix cu_H_opt(cu_Tkp1sk_C1t.dim1, cu_inv_Sigmatot.dim2);

			// H_opt = Tkp1sk_C1t * inv_Sigmatot
			kp_cu_gemm(cublasHandle, 'N', 'N', 1 , cu_Tkp1sk_C1t, cu_inv_Sigmatot , 0, cu_H_opt);
        	
			cu_Tkp1sk_C1t.resize(0,0);
			cu_inv_Sigmatot.resize(0,0);
			
			kp_cu_matrix cu_Atur_Hopt(nb_az, nb_p);
			// Atur_Hopt = Atur * H_opt
			/*for(i = 0 ; i < cu_atur.size() ; i++)
				for (j = 0 ; j < cu_Atur_Hopt.dim2 ; j++)
					cu_Atur_Hopt.d_cu[j * cu_Atur_Hopt.dim1 + i] = cu_atur.d_cu[i] * cu_H_opt.d_cu[j * cu_H_opt.dim1 + i];*/
			kernel_diag_mult(cu_Atur_Hopt.d_cu, cu_H_opt.d_cu, cu_atur.d_cu, cu_atur.size(), cu_H_opt.dim1*cu_H_opt.dim2);

			cu_H_inf.resize(nb_n,nb_p);
			cu_ms.push_back(&cu_Atur_Hopt) ; cu_ms.push_back(&cu_H_opt);
			// H_inf = [ Atur*H_opt ; H_opt]
			kp_cu_vertcat(cu_H_inf, cu_ms);
			cu_ms.clear();
			cu_Atur_Hopt.resize(0,0);
			cu_H_opt.resize(0,0);
		}
		
		else
		{
			cu_H_inf.resize(nb_n,nb_p);
			kp_cu_gemm(cublasHandle, 'N', 'N', 1 , cu_Tkp1sk_C1t, cu_inv_Sigmatot , 0, cu_H_inf);
			cu_Tkp1sk_C1t.resize(0,0);
			cu_inv_Sigmatot.resize(0,0);
		}
	
	
	}




	/*kp_matrix H_inf;
	kp_cu2kp_matrix(H_inf, cu_H_inf);
	fichier.open("H_inf.dat",ios::out);
	for(i=0;i<H_inf.dim1;i++)
	{
		for (j=0;j<H_inf.dim2;j++)
		{
			fichier<< __SP H_inf(i,j)<<" ";
		}
		fichier << endl;
	}	
	fichier.close();*/


	gainComputed = true;
	




}




void kp_kalman_core_full_GPU::next_step(const kp_vector& Y_k, kp_vector& U_k)
{

	if(!gainComputed)
	{
		cerr << "Error | kp_kalman_core_full_GPU::next_step | gain has not been initialized" << endl;
		exit(EXIT_FAILURE);
	}

	real mean_Xkp1skdebut;
	cu_Y_k = Y_k;
	cu_Nact_Ukm2.zeros();
	cu_tmp_vec1.zeros();
	cu_innovation.zeros();
	cu_X_kskm1_tmp.zeros();
	cu_X_kp1sk.zeros();
	cu_Y_kskm1.zeros();
	cu_A1_00_Xkdebut.zeros();
	cu_A1_01_Xkfin.zeros();



//temps_op1.start();
	// VECTEUR d'ESTIMATION de MESURE ( A l' INSTANT K )
	// Nact_Ukm2 = N_Act * U_km2 
	//kp_gemv (1, N_Act, *U_km2, 0, *Nact_Ukm2);
	kp_cu_gemv (cublasHandle, 'N', 1, cu_N_Act, cu_U_km2, 0, cu_Nact_Ukm2);
	
	// tmp_vec1 = X_kskm1 - Nact_Ukm2 (= X_kskm1 - N_Act * U_km2)
	//kp_cu_cudaMemcpy(tmp_vec1->d, cu_X_kskm1->d_cu+nb_az, nb_az*sizeof(real), cudaMemcpyDeviceToHost);
	kernel_memcpy_real(cu_tmp_vec1.d_cu, cu_X_kskm1.d_cu+nb_az, nb_az);
	cu_tmp_vec1 -= cu_Nact_Ukm2 ;

		
	// Y_kskm1 = D_Mo * tmp_vec1 (= D_Mo * (X_kskm1 - N_Act * U_km2))
	//kp_gemv (1,D_Mo, *tmp_vec1,0,*Y_kskm1); 
	kp_cu_gemv (cublasHandle, 'N', 1, cu_D_Mo, cu_tmp_vec1, 0, cu_Y_kskm1); 
//temps_op1.pause();



		
//temps_op2.start();			
	// VECTEUR D'ESTIMATION de PREDICTION ( A l' INSTANT K )
		
	// innovation = Y_k - Y_kskm1
	cu_innovation = cu_Y_k;
	//cu_Y_kskm1 = *Y_kskm1;
	cu_innovation -= cu_Y_kskm1;
	cu_X_kskm1_tmp = cu_X_kskm1;




	/*// X_kskm1_tmp = H_inf *  (Y_k - Y_kskm1)
	kp_gemv ('N',1,H_inf,*innovation,1,*X_kskm1_tmp); */
	kp_cu_gemv (cublasHandle, 'N', 1, cu_H_inf, cu_innovation, 1, cu_X_kskm1_tmp);



	// X_kskm1_tmp = X_kskm1_tmp + H_inf * innovation (= X_kskm1 + H_inf * (Y_k - Y_kskm1))

	//kp_cu_gemv (cublasHandle,'N',1, *cu_A1, cu_X_kskm1_tmp,0,cu_X_kp1sk);
	kernel_memset_real(cu_X_kp1sk.d_cu, 0.0, nb_az);
	kernel_memcpy_real(cu_A1_00_Xkdebut.d_cu, cu_X_kskm1_tmp.d_cu, nb_az);
	kernel_memcpy_real(cu_A1_01_Xkfin.d_cu, cu_X_kskm1_tmp.d_cu+nb_az, nb_az);

	cu_A1_00_Xkdebut *= cu_atur;

	cu_A1_01_Xkfin *= cu_btur;
	kernel_memcpy_real(cu_X_kp1sk.d_cu, cu_A1_00_Xkdebut.d_cu, nb_az);
	kernel_add_real(cu_X_kp1sk.d_cu, cu_A1_01_Xkfin.d_cu, nb_az);
	kernel_memcpy_real(cu_X_kp1sk.d_cu + nb_az, cu_X_kskm1_tmp.d_cu, nb_az);


	//init_from_kp_cu_vector2kp_vector(*X_kp1sk_debut, cu_X_kp1sk,0 , nb_az);
	kernel_memcpy_real(cu_X_kp1sk_debut.d_cu, cu_X_kp1sk.d_cu, nb_az);
	if (isZonal)
	{
		mean_Xkp1skdebut = kp_cu_reduce(cu_X_kp1sk_debut)/nb_az;
		// X_kp1sk_tmp = X_kp1sk(1:nb_az)-mean(X_kp1sk(1:nb_az))*ones(nb_az,1)
		cu_X_kp1sk_tmp = cu_X_kp1sk_debut; 
		cu_X_kp1sk_tmp -= mean_Xkp1skdebut; 
	}
//temps_op2.pause();




	
//temps_op3.start();
	//TENSION de CORRECTION


	//kp_gemv (1,PROJ, *X_kp1sk_tmp, 0, *U_k);
	//kp_cu_vector* cu_U_k = new kp_cu_vector(nb_act);

	if (isZonal)
		kp_cu_gemv(cublasHandle, 'N', 1, cu_PROJ, cu_X_kp1sk_tmp, 0, cu_U_k);
	else
		kp_cu_gemv(cublasHandle, 'N', 1, cu_PROJ, cu_X_kp1sk_debut, 0, cu_U_k);

	//kp_cu_gemv(cublasHandle, 'N', 1, *cu_PROJ_full, cu_X_kp1sk_tmp, 0, cu_U_k);
//temps_op3.pause();



	
	//MISE A JOUR
	cu_U_km2 = cu_U_km1;
	cu_U_km1 = cu_U_k;
	cu_X_kskm1 = cu_X_kp1sk;
	//*U_km2 = *U_km1; 
	//*U_km1 = *U_k; 
	//*X_kskm1 = *X_kp1sk;
		
	// On doit donner cu_U_k, mais comme on vient de faire la mise a jour juste 
	// avant au lieu de la faire apres, on donne cu_U_km1
	kp_cu2kp_vector(U_k, cu_U_km1);	








}
