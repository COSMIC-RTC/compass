//kp_kalman_core_full_CPU.cpp

#include "kp_kalman_core_full_CPU.h"

#include <fstream>
#include <iomanip>
#define __SP setprecision(20)<<

kp_kalman_core_full_CPU::kp_kalman_core_full_CPU(const kp_matrix<KFPP>& D_Mo_,
		const kp_matrix<KFPP>& N_Act_,
		const kp_matrix<KFPP>& PROJ_,
		bool isZonal_) : kp_kalman_core_full(D_Mo_, N_Act_, PROJ_, isZonal_)
{
	H_inf.resize(0,0);
	D_Mo = D_Mo_;

	N_Act = N_Act_;
	PROJ = PROJ_;
	U_km2.resize(nb_act);
	U_km2.zeros();
	U_km1.resize(nb_act);
	U_km1.zeros();
	X_kskm1.resize(nb_n);
	X_kskm1.zeros();
	
	// variables de next_step
	Nact_Ukm2.resize(nb_az);
	tmp_vec1.resize(nb_az);
	Y_kskm1.resize(nb_p);
	innovation.resize(nb_p);
	X_kskm1_tmp.resize(nb_n);
	A1_00_Xkdebut.resize(nb_az);
	A1_01_Xkfin.resize(nb_az);
	X_kp1sk_debut.resize(nb_az);
	X_kp1sk_tmp.resize(nb_az);
	X_kp1sk.resize(nb_n);
	
}


void kp_kalman_core_full_CPU::calculate_gain(float bruit_pix,
		float k_W, 
		const  kp_matrix<double>&  SigmaV,
		const  kp_vector<double>& atur_,
		const  kp_vector<double>& btur_)
{

	//atur = atur_;
	//btur = btur_;
	//

	if (atur_.size() != nb_az)
	{
		cerr<<"Error | kp_kalman_core_full_CPU::calculate_gain | size problem atur"<<endl;
		exit(EXIT_FAILURE);
	}
	if (btur_.size() != nb_az)
	{
		cerr<<"Error | kp_kalman_core_full_CPU::calculate_gain | size problem btur"<<endl;
		exit(EXIT_FAILURE);
	}
	if ((SigmaV.getDim1() != SigmaV.getDim2()) || (SigmaV.getDim1() != nb_az))
	{
		cerr<<"Error | kp_kalman_core_full_CPU::calculate_gain | size problem SigmaV"<<endl;
		exit(EXIT_FAILURE);
	}

	bool AR1 = true;
	for(int i = 0 ; i < btur_.size() ; i++) AR1 &= (btur_[i]==0);
	if (AR1) ordreAR = 1 ; else ordreAR = 2;
	int expression = 2;
	KFPP seuil = 1/pow(10,14);
	int boucle = 0;
	KFPP ecart = 1.0;
	const int boucle_max = 50;
	const KFPP SigmaW = k_W*bruit_pix;
	//KFPP Trac_T[boucle_max];
	KFPP trac1_tmp, trac2_tmp;

        kp_matrix<double> H_inf2(0,0);
        if (!gainComputed) 
		H_inf.resize(0,0);

//kp_timer temps_alphak, temps_betak, temps_Tk, temps_inversion ;



	vector<kp_matrix<double>*> ms ;
	
	kp_matrix<double> alpha_kp1;
	if (ordreAR == 1) alpha_kp1.resize(nb_az,nb_az);
	else alpha_kp1.resize(nb_n,nb_n);

	alpha_kp1.zeros();
	
	// pour ordreAR 1 et 2
	for (int i = 0 ; i < nb_az ; i++)
	{
		alpha_kp1(i,i) = atur_[i];
	}
	//pour ordreAR 2 uniquement
	if (ordreAR == 2)
	{
		for (int i = 0 ; i < nb_az ; i++)
		{
			alpha_kp1(i,i+nb_az) = 1.0;
			alpha_kp1(i+nb_az,i) = btur_[i];
		}

	}

	kp_matrix<double> beta_kp1;
	if (ordreAR == 1) beta_kp1.resize(nb_az, nb_az);
	else beta_kp1.resize(nb_n, nb_n);
	beta_kp1.zeros();

	// zeros_Dmo = [0 D_Mo]
	kp_matrix<double>* zeros_Dmo = new kp_matrix<double>();
	kp_matrix<double>* zeros_nbp_nbaz = new kp_matrix<double>();
	zeros_nbp_nbaz->resize(nb_p, nb_az);
	zeros_nbp_nbaz->zeros();
	ms.push_back(zeros_nbp_nbaz);
	kp_matrix<double> D_Mo_double(D_Mo);
	ms.push_back(&D_Mo_double);
	kp_horizcat(*zeros_Dmo,ms);
	ms.clear();
	delete zeros_nbp_nbaz ; zeros_nbp_nbaz = NULL;

	kp_matrix<double>  C1; 
	// si ordreAR 1 C1 = D_Mo
	if (ordreAR == 1) C1 = D_Mo_double;
	//si ordreAR 2 C1 = [zeros(nb_p,nb_az) D_Mo]
	else C1 = *zeros_Dmo;

	D_Mo_double.resize(0,0);

	// VRAI SEULEMENT SI SigmaW = k*Id (k reel)	
	// beta_kp1 = (C1)T * C1
	kp_gemm ('T','N', 1, C1, C1, 0, beta_kp1);
	// beta_kp1 = (C1)T * C1 / SigmaW (avec SigmaW reel)


	beta_kp1 *= 1/SigmaW ;

		



	kp_matrix<double> T_kp1sk;
	if (ordreAR == 1) T_kp1sk = SigmaV;
	else
	{
		// beta_kp1_tmp = [SigmaV ; zeros(nb_az, nb_az)];
		kp_matrix<double>* zeros_nbaz_nbaz = new kp_matrix<double>(nb_az, nb_az);
		zeros_nbaz_nbaz->zeros();
		kp_matrix<double>* T_kp1sk_tmp = new kp_matrix<double>(SigmaV);
		ms.push_back(T_kp1sk_tmp); ms.push_back(zeros_nbaz_nbaz);
		kp_matrix<double>* beta_kp1_tmp = new kp_matrix<double>(nb_n, nb_az);
		kp_vertcat(*beta_kp1_tmp, ms);
		delete zeros_nbaz_nbaz ; zeros_nbaz_nbaz = NULL;
		delete T_kp1sk_tmp; T_kp1sk_tmp = NULL;
		ms.clear();

		// T_kp1sk = [beta_kp1_tmp zeros_nbn_nbaz];
		kp_matrix<double>* zeros_nbn_nbaz = new kp_matrix<double>(nb_n, nb_az);
		zeros_nbn_nbaz->zeros();
		ms.push_back(beta_kp1_tmp); ms.push_back(zeros_nbn_nbaz);
		kp_horizcat(T_kp1sk, ms);
		delete beta_kp1_tmp ; beta_kp1_tmp = NULL;
		delete zeros_nbn_nbaz ; zeros_nbn_nbaz = NULL;
		ms.clear();
	}



	kp_matrix<double> alpha_k(alpha_kp1.getDim1(), alpha_kp1.getDim2());
	kp_matrix<double> T_k(T_kp1sk.getDim1(), T_kp1sk.getDim2());
	kp_matrix<double> beta_k(beta_kp1.getDim1(), beta_kp1.getDim2());
	
	kp_matrix<double> IBG_1(alpha_k.getDim2(), beta_k.getDim1());



	kp_matrix<double> Tk_IBG1(T_k.getDim1(), IBG_1.getDim2());
	kp_matrix<double> betak_Tk(beta_k.getDim1(), T_k.getDim2());
	kp_matrix<double> alphak_IBG1(alpha_k.getDim1(), IBG_1.getDim2());
	kp_matrix<double> alphak_IBG1_betak(alpha_k.getDim1(), beta_k.getDim2());
	kp_matrix<double> Tk_IBG1_alphak(T_k.getDim1(), alpha_k.getDim2());




	while( (ecart>seuil) && (boucle < boucle_max) )
	{
		// mise a jour
		alpha_k = alpha_kp1;
		beta_k = beta_kp1;
		T_k = T_kp1sk;


		// Calcul de IBG_1
		// betak_Tk = beta_k * T_k
		
		kp_gemm ('N', 'N', 1 , beta_k, T_k, 0, betak_Tk);

		
		if (ordreAR ==1)
			for (int i=0 ; i<nb_az ; i++) (betak_Tk(i,i)) += 1;
		else
		{
			for (int i=0 ; i<nb_n ; i++) (betak_Tk(i,i)) += 1;
		}

		IBG_1 = betak_Tk;

//temps_inversion.start();
		IBG_1.inverse();
//temps_inversion.pause();
		



	
//temps_alphak.start();
		// Calcul de alpha_kp1
		// alphak_IBG1 = alpha_k * IBG_1
		kp_gemm ('N', 'N', 1 , alpha_k, IBG_1 , 0, alphak_IBG1);


        	// alpha_kp1 = alphak_IBG1 * alphak (= alpha_k * IBG1 * alpha_k) 
		kp_gemm ('N', 'N', 1 , alphak_IBG1, alpha_k , 0, alpha_kp1);
	
//temps_alphak.pause();	








//temps_betak.start();
		// Calcul de beta_kp1

		// alphak_IBG1_betak = alphak_IBG1 * beta_k (= alpha_k * IBG1 * beta_k)
		kp_gemm ('N', 'N', 1 , alphak_IBG1, beta_k , 0, alphak_IBG1_betak);


		//beta_kp1 = alphak_IBG1_betak * (alpha_k)T (= alpha_k * IBG1 * beta_k * (alpkha_k)T)
		kp_gemm ('N', 'T', 1 , alphak_IBG1_betak, alpha_k , 1, beta_kp1);

//temps_betak.pause();






//temps_Tk.start();
		//Calcul de T_kp1sk

		// Tk_IBG1 = T_k * IBG1
		kp_gemm ('N', 'N', 1, T_k ,IBG_1, 0, Tk_IBG1);
		// Tk_IBG1_alphak = Tk_IBG1 * alpha_k (= T_k * IBG1 * alpha_k)
		kp_gemm ('N', 'N', 1 , Tk_IBG1 , alpha_k , 0, Tk_IBG1_alphak);		

		// T_kp1sk = (alpha_k)T * Tk_IBG1_alphak (= (alpha_k)T * T_k * IBG1 * alpha_k)
		kp_gemm ('T', 'N', 1 , alpha_k, Tk_IBG1_alphak , 1, T_kp1sk);		

//temps_Tk.pause();





		trac1_tmp=0;trac2_tmp=0;
		for (int i=0 ; i< T_kp1sk.getDim1() ; i++)
		{
			trac1_tmp += T_kp1sk(i,i);
			trac2_tmp += T_k(i,i);
		}
		//Trac_T[boucle]=trac1_tmp;
		ecart =fabs(trac1_tmp/trac2_tmp-1.0);
		boucle ++;

	}

        if (boucle >= boucle_max)
        {
                cerr << "Le gain de Kalman ne converge pas.";
                exit(EXIT_FAILURE);
        }


/*cout<<"nb boucle = "<<boucle<<endl;
cout<<"temps inversion " << temps_inversion.rez()<<endl;
cout<< "temps alpha_k = "<<temps_alphak.rez()<<endl;
cout<< "temps beta_k = "<<temps_betak.rez()<<endl;
cout<< "temps T_k = "<<temps_Tk.rez()<<endl;*/

	
	/*fichier.open("Tkp1sk_full_CPU.dat",ios::out);
	for (int i=0 ; i<T_kp1sk.getDim1() ; i++)
	{
		for(int j=0 ; j<T_kp1sk.getDim2() ; j++)
		{
			fichier << __SP T_kp1sk(i,j) << " ";
		}
		fichier<<endl;
	}
	fichier.close();*/

	


	if ( (expression == 1) && (ordreAR == 1))
	{
		//calcul de Sinf_0_0 (matrice superieure gauche de S_inf)
        	kp_matrix<double>* Atur_Tkp1sk = new kp_matrix<double>(nb_az,nb_az);
		kp_matrix<double>* Sinf_0_0 = new kp_matrix<double>(nb_az,nb_az);
		// Atur_Tkp1sk = Atur * T_kp1sk
		for(int i = 0 ; i < atur_.size() ; i++)
		{
			for (int j = 0 ; j < Atur_Tkp1sk->getDim2() ; j++)
			{
				(*Atur_Tkp1sk)(i,j) = atur_[i] * T_kp1sk(i,j);
				(*Sinf_0_0)(i,j) = atur_[i] * T_kp1sk(i,j) * atur_[j];
			}
		}
		
		// Sinf_0_0 = Sinf_0_0 + SigmaV
		*Sinf_0_0 += SigmaV ;
		
		
		// calcul de Sinf_0 (matrice superieure de S_inf)
		kp_matrix<double>* Sinf_0 =  new kp_matrix<double>(nb_az,nb_n);
		ms.push_back(Sinf_0_0) ; ms.push_back(Atur_Tkp1sk) ;
		// Sinf = [ Sinf_0_0  A_tur*T_kp1sk ] 
		kp_horizcat(*Sinf_0,ms);
		ms.clear();
		delete Sinf_0_0 ; Sinf_0_0 = NULL;
		delete Atur_Tkp1sk ; Atur_Tkp1sk = NULL;


		// calcul de Sinf_1 (matrice inferieure de S_inf)
		kp_matrix<double>* Sinf_1_0 = new kp_matrix<double>(nb_az,nb_az);
		// Sinf_1_0 = T_kp1sk * (Atur)T (=T_kp1sk * (Atur)T car Atur diagonale)
		for(int i = 0 ; i < atur_.size() ; i++)
			for (int j = 0 ; j < Sinf_1_0->getDim2() ; j++)
				(*Sinf_1_0)(i,j) = T_kp1sk(i,j) * atur_[j];

		
		kp_matrix<double>* Sinf_1 = new kp_matrix<double>(nb_az,nb_n);
		ms.push_back(Sinf_1_0) ; ms.push_back(&T_kp1sk) ;
		// Sinf_1 = [ Sinf_1_0  T_kp1sk] 
		kp_horizcat(*Sinf_1,ms);
		ms.clear();
		delete Sinf_1_0 ; Sinf_1_0 = NULL;
	
		// S_inf = [ Sinf_0 ; Sinf_1 ]
		kp_matrix<double>* S_inf = new kp_matrix<double>(nb_n, nb_n);
		ms.push_back(Sinf_0) ; ms.push_back(Sinf_1) ;	 
		kp_vertcat(*S_inf, ms);	
		ms.clear();
		delete Sinf_0 ; Sinf_0 = NULL;
		delete Sinf_1 ; Sinf_1 = NULL;
	







		kp_matrix<double>* Sinf_zerosDmot = new kp_matrix<double>(nb_n, nb_p);
		// Sinf_zerosDmot = Sinf * [0 D_Mo]T
		kp_gemm ('N', 'T', 1 , *S_inf, *zeros_Dmo , 0, *Sinf_zerosDmot);



		// calcul de Sigma_tot = [0 D_Mo] * S_inf * [0 D_Mo]' + SigmaW
		kp_matrix<double>* Sigma_tot = new kp_matrix<double>(nb_p,nb_p);


		// Sigma_tot = zeros_Dmo * Sinf_zerosDmot  ( = [0 D_Mo] * Sinf * [0 ; (D_Mo)T] )
		kp_gemm ('N', 'N', 1 , *zeros_Dmo, *Sinf_zerosDmot , 0, *Sigma_tot);


		// Sigma_tot = Sigma_tot + SigmaW*Id (avec SigmaW reel)  (= [0 D_Mo] * Sinf * [0 ; (D_Mo)T] + SigmaW*Id)
		for(int i = 0 ; i < Sigma_tot->getDim1() ; i++) (*Sigma_tot)(i,i) += SigmaW;
		delete zeros_Dmo ; zeros_Dmo = NULL;
		
		//inversion de Sigma_tot
		kp_matrix<double>* inv_Sigmatot = Sigma_tot;
		Sigma_tot = NULL;
		// inv_Sigmatot = inv(Sigma_tot)
		inv_Sigmatot->inverse();
		// ATTENTION !!!! Le pointeur Sigma_tot n'existe plus
	
		H_inf2.resize(nb_n,nb_p);		
		// H_inf = Sinf_zerosDmot * inv_Sigmatot  ( = Sinf * [0 (D_Mo)T] * inv(Sigma_tot))
		kp_gemm ('N', 'N', 1 , *Sinf_zerosDmot, *inv_Sigmatot , 0, H_inf2);

		delete Sinf_zerosDmot ; Sinf_zerosDmot = NULL;
		delete inv_Sigmatot ; inv_Sigmatot = NULL;
		
		
	}
	else
	{
        	//Calcul de T_kp1sk * (C1)T
        	kp_matrix<double>* Tkp1sk_C1t = new kp_matrix<double>();
		if(ordreAR==1) Tkp1sk_C1t->resize(nb_az,nb_p);
		else Tkp1sk_C1t->resize(nb_n,nb_p);
		kp_gemm ('N', 'T', 1 , T_kp1sk, C1 , 1, *Tkp1sk_C1t);


		//Calcul de C1 * T_kp1sk * C1' + SigmaW
        	kp_matrix<double>* Sigma_tot = new kp_matrix<double>(C1.getDim1(), Tkp1sk_C1t->getDim2());
		// Sigma_tot = SigmaW * Id (avec SigmaW reel)
		Sigma_tot->zeros();
		for(int i = 0 ; i < Sigma_tot->getDim1() ; i++) (*Sigma_tot)(i,i) += SigmaW;


		// Sigma_tot = Sigma_tot + C1 * Tkp1sk_C1t  ( = (C1 * T_kp1sk * (C1)T) + SigmaW*Id )
		kp_gemm ('N', 'N', 1 , C1, *Tkp1sk_C1t , 1, *Sigma_tot);

			
		// inversion de (Sigmatot)
        	kp_matrix<double>* inv_Sigmatot = Sigma_tot;
		Sigma_tot = NULL;
		//inv_Sigmatot = inv(inv_Sigmatot)
        	inv_Sigmatot->inverse();
		// ATTENTION !!!! Le pointeur Sigma_tot n'existe plus
        
		
		if (ordreAR == 1)
		{
			kp_matrix<double>* H_opt = new kp_matrix<double>(Tkp1sk_C1t->getDim1(), inv_Sigmatot->getDim2());


			// H_opt = Tkp1sk_C1t * inv_Sigmatot
			kp_gemm	('N', 'N', 1 , *Tkp1sk_C1t, *inv_Sigmatot , 0, *H_opt);


			delete Tkp1sk_C1t ; Tkp1sk_C1t = NULL;
			delete inv_Sigmatot ; inv_Sigmatot = NULL;
			
			kp_matrix<double>* Atur_Hopt = new kp_matrix<double>(nb_az, nb_p);
			// Atur_Hopt = Atur * H_opt
			for(int i = 0 ; i < atur_.size() ; i++)
				for (int j = 0 ; j < Atur_Hopt->getDim2() ; j++)
					(*Atur_Hopt)(i,j) = atur_[i] * (*H_opt)(i,j);

			ms.push_back(Atur_Hopt) ; ms.push_back(H_opt);
			// H_inf = [ Atur*H_opt ; H_opt]
			H_inf2.resize(nb_n,nb_p);		
			kp_vertcat(H_inf2,ms);
			ms.clear();
			delete Atur_Hopt ; Atur_Hopt = NULL;
			delete H_opt ; H_opt = NULL;
		}
		
		else
		{
			H_inf2.resize(nb_n,nb_p);		
			kp_gemm	('N', 'N', 1 , *Tkp1sk_C1t, *inv_Sigmatot , 0, H_inf2);
		}
	
	
	}
	
	atur = atur_;
	btur = btur_;
	H_inf = H_inf2;

	/*fichier.open("H_inf_full_CPU.dat",ios::out);
	for (int i=0 ; i<H_inf.getDim1() ; i++)
	{
		for(int j=0 ; j<H_inf.getDim2() ; j++)
		{
			fichier << __SP H_inf(i,j) << " ";
		}
		fichier<<endl;
	}
	fichier.close();*/

	gainComputed = true;




}




void kp_kalman_core_full_CPU::next_step(const kp_vector<KFPP>& Y_k, kp_vector<KFPP>& U_k)
{
	if(!gainComputed)
	{
		cerr << "Error | kp_kalman_core_full_CPU::next_step | gain has not been initialized"<<endl;
		exit(EXIT_FAILURE);
	}

	KFPP mean_Xkp1skdebut;

	X_kskm1_tmp.zeros();
	A1_00_Xkdebut.zeros();
	A1_01_Xkfin.zeros();


//temps_op1.start();
	// VECTEUR d'ESTIMATION de MESURE ( A l' INSTANT K )
	// Nact_Ukm2 = N_Act * U_km2 
	kp_gemv ('N', 1, N_Act, U_km2, 0, Nact_Ukm2);


	// tmp_vec1 = X_kskm1 - Nact_Ukm2 (= X_kskm1 - N_Act * U_km2)
	for(int i = 0 ; i < nb_az ; i++)
	{
		tmp_vec1[i] = X_kskm1[nb_az+i] + Nact_Ukm2[i] ;
	}

	// Y_kskm1 = D_Mo * tmp_vec1 (= D_Mo * (X_kskm1 - N_Act * U_km2))
	kp_gemv ('N', 1, D_Mo, tmp_vec1, 0, Y_kskm1); 
//temps_op1.pause();



//temps_op2.start();			
	// VECTEUR D'ESTIMATION de PREDICTION ( A l' INSTANT K )

	// innovation = Y_k - Y_kskm1
	innovation = Y_k;

	innovation += Y_kskm1;
	X_kskm1_tmp = X_kskm1; 

	// X_kskm1_tmp = H_inf *  (Y_k - Y_kskm1)
	kp_gemv ('N', 1, H_inf, innovation, 1, X_kskm1_tmp); 

	// X_kskm1_tmp = X_kskm1_tmp + H_inf * innovation (= X_kskm1 + H_inf * (Y_k - Y_kskm1))
		
	X_kp1sk.zeros();
	A1_00_Xkdebut.init_from_vector(X_kskm1_tmp, 0, nb_az);
	A1_01_Xkfin.init_from_vector(X_kskm1_tmp, nb_az, nb_n);
	A1_00_Xkdebut *= atur;
	A1_01_Xkfin *= btur;
	
	for (int i=0;i<nb_az;i++)
	{
		X_kp1sk[i] = A1_00_Xkdebut[i] + A1_01_Xkfin[i];
		X_kp1sk[i+nb_az] =  X_kskm1_tmp[i];
	}
	

	//X_kp1sk_debut = X_kp1sk(1:nb_az)
	for(int i = 0 ; i < nb_az ; i++) X_kp1sk_debut[i] = X_kp1sk[i];

	// si on est en zonal, on retire le piston
	if (isZonal)
	{
		// mean_Xkp1skdebut = mean (X_kp1sk_debut) (= mean(X_kp1sk(1:nb_az)))
		mean_Xkp1skdebut = X_kp1sk_debut.mean();
	
		// X_kp1sk_tmp = X_kp1sk(1:nb_az)-mean(X_kp1sk(1:nb_az))*ones(nb_az,1)
		X_kp1sk_tmp = X_kp1sk_debut; 
		X_kp1sk_tmp -= mean_Xkp1skdebut; 
	}
//temps_op2.pause();
	


U_k.resize(nb_act);
//temps_op3.start();
	//TENSION de CORRECTION
	if (isZonal)
		kp_gemv ('N', -1, PROJ, X_kp1sk_tmp, 0, U_k); 
	else
		kp_gemv ('N', -1, PROJ, X_kp1sk_debut, 0, U_k);  


//temps_op3.pause();

	//MISE A JOUR
	
	U_km2 = U_km1; 
	U_km1 = U_k; 
	X_kskm1 = X_kp1sk;
}
