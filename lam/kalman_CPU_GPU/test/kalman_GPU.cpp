#include "kalman_GPU.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#define __SP setprecision(20)<<

//#define RAND (rand())/((real) RAND_MAX)
//#define RANDN (sqrt(-2.0*log(RAND))*cos(2*M_PI*RAND))


int main()
{
	cudaSetDevice(0);	
	//cudaDeviceReset();
	const bool afficher_temps = 1;
	const bool comparaison_matlab = 0;
	const bool sparse_matrix = 1  ;
	bool isZonal = true;

	cusparseHandle_t cusparseHandle;
  	cublasHandle_t cublasHandle;

	cublasStatus_t cublasStat;
	cusparseStatus_t cusparseStat;

	cublasStat = cublasCreate(&cublasHandle);
	if (cublasStat != CUBLAS_STATUS_SUCCESS) 
	{
		cerr<<"Error | kp_kalman_core_full_GPU::kp_kalman_core_full_GPU | CUBLAS initialization failed"<<endl;
		exit(EXIT_FAILURE);
	}

	cusparseStat = cusparseCreate(&cusparseHandle);
	if (cusparseStat != CUSPARSE_STATUS_SUCCESS) 
	{
		cerr<<"Error | kp_kalman_core_sparse_GPU::kp_kalman_core_sparse_GPU | CUSPARSE Library initialization failed"<<endl;
		exit(EXIT_FAILURE);
	}
	magma_init();

	

	gsl_rng* rng;
	int i,j;
	real E_coh=0;
	
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, time(NULL));
	
	kp_timer temps_total, temps_boucle;

	kp_cu_timer temps_riccati_cu;

	kp_cu_timer temps_effectif;

	if (afficher_temps) temps_total.start();

	cout <<("---DEBUT EXECUTION KALMAN---")<<endl;
	
	ofstream fichier,fichier2;
	mat_t *mat,*mat2;

	
	//CHARGEMENT DES MATRICES (load de Matlab)
	
	kp_smatrix* D_Mo_s  = new kp_smatrix(); 
	kp_smatrix* N_Act_s = new kp_smatrix();
	kp_smatrix* PROJ_s  = new kp_smatrix();
	kp_matrix* D_Mo_f   = new kp_matrix(); 
	kp_matrix* N_Act_f  = new kp_matrix();
	kp_matrix* PROJ_f   = new kp_matrix();


	kp_smatrix D_Sys;
	kp_smatrix* A1 = new kp_smatrix();
	kp_smatrix N_Ph;

	kp_matrix SigmaV;

	int nb_z  = -1;
	int nb_az = -1;

	if (isZonal) mat = Mat_Open("/home/tgautrais/Documents/v2_kalman_kp/data/KALM_08P16_AR1Za980.mat",MAT_ACC_RDONLY);
	else mat = Mat_Open("/home/tgautrais/Documents/v2_kalman_kp/data/KALM_08P16_AR1Mo30.mat",MAT_ACC_RDONLY);
	//LECTURE DE D_Sys
	kp_init4matio_smatrix(D_Sys,mat,"D_Sys");
	//LECTURE DE A1
	kp_init4matio_smatrix(*A1,mat,"A1");
	//LECTURE DE SigmaV
	kp_init4matio_matrix(SigmaV,mat,"SigmaV");
	//LECTURE DE D_Mo, N_Act et PROJ
	if (isZonal) 
	{
		kp_init4matio_smatrix(*D_Mo_s, mat, "D_Mo");
		kp_init4matio_smatrix(*N_Act_s, mat, "N_Act");
		kp_init4matio_smatrix(*PROJ_s, mat, "PROJ");
	}
	else 
	{
		kp_init4matio_matrix(*D_Mo_f, mat, "D_Mo");
		kp_init4matio_matrix(*N_Act_f, mat, "N_Act");
		kp_init4matio_matrix(*PROJ_f, mat, "PROJ");
	}
	//LECTURE DE N_Ph
	kp_init4matio_smatrix(N_Ph,mat,"N_Ph");
	//LECTURE DES DIMENSIONS
	const int NB_ACT = kp_init4matio_int(mat, "nb_act");
	const int NB_N = kp_init4matio_int(mat, "nb_n");
	const int NB_P = kp_init4matio_int(mat, "nb_p");
	const int NB_PX = kp_init4matio_int(mat, "nb_px");
	if (!isZonal) 
	{
		nb_z   = kp_init4matio_int(mat, "nb_z");
		nb_az = nb_z;
	}
	else nb_az = NB_ACT;
	Mat_Close(mat);
	const int NB_Z  = nb_z;
	const int NB_AZ = nb_az;
	//INITIALISATIONS et PARAMETRES
	const int L0 = 25;
	const int i0var = 51;

	const real BRUIT_RAD = 0.04;
	const real k_W = 5;
	const real BRUIT_PIX = BRUIT_RAD/(M_PI*M_PI);
	const real SQRT_BRUIT_PIX = sqrt(BRUIT_PIX);
	const int NB_BOUCLE = 5000;

	if (!sparse_matrix && isZonal)
	{
		D_Mo_f->init_from_smatrix(*D_Mo_s);
		N_Act_f->init_from_smatrix(*N_Act_s);
		PROJ_f->init_from_smatrix(*PROJ_s);
	}
	if (sparse_matrix && (!isZonal))
	{
		D_Mo_s->init_from_matrix( *D_Mo_f);
		N_Act_s->init_from_matrix( *N_Act_f);
		PROJ_s->init_from_matrix( *PROJ_f);
	}

	if (sparse_matrix)
	{
		D_Mo_s ->resize2rowMajor();
		N_Act_s->resize2rowMajor();
		PROJ_s ->resize2rowMajor();
	}

	kp_kalman_core_sparse_GPU kalman_GPU(*D_Mo_s, *N_Act_s, *PROJ_s, isZonal, cublasHandle, cusparseHandle);
	//kp_kalman_core_full_GPU kalman_GPU(*D_Mo_f, *N_Act_f, *PROJ_f, isZonal, cublasHandle);

	
	// Atur
	kp_vector A1_00(NB_AZ) ; A1_00.zeros();
	// Btur (different de zero si AR2)
	kp_vector A1_01(NB_AZ) ; A1_01.zeros();
	A1_2vec(A1_00, A1_01, *A1);
	delete A1 ; A1 = NULL;


	real* Trac_T;
	int boucle_riccati = 0;
	
	
	//Initialisation boucle OA
	kp_vector Var_PhRes(NB_BOUCLE) ; Var_PhRes.zeros();
	kp_vector U_km1(NB_ACT) ; U_km1.zeros();
	kp_vector U_k(NB_ACT) ; U_k.zeros();

	kp_vector CPC_km1(NB_PX) ; CPC_km1.zeros();
	kp_vector CPC_k(NB_PX) ; CPC_k.zeros();
	kp_vector CPC_kp1(NB_PX) ; CPC_kp1.zeros();
	kp_vector X_kskm1(NB_N) ; X_kskm1.zeros();
	kp_vector X_kp1sk(NB_N) ; X_kp1sk.zeros();
	kp_vector CPT_km1(NB_PX);
	kp_vector CPR_km1(NB_PX);
	kp_vector Y_k(NB_P);
	kp_vector randn(NB_P);
	kp_vector racSigmaW_randn(NB_P);
	kp_vector Y_kskm1(NB_P);

	
	real mean_CPCkp1; 





	


	if (afficher_temps) temps_riccati_cu.start();
		kalman_GPU.calculate_gain(BRUIT_PIX, k_W, SigmaV, A1_00, A1_01);

	if (afficher_temps)
	{
		temps_riccati_cu.pause();
		cout << "   Temps riccati = " << temps_riccati_cu.rez() << endl;
	}


	/*kp_matrix H_inf_tmp;
	mat2 = Mat_Open("/home/tgautrais/Documents/v2_kalman_kp/data/H_inf.mat",MAT_ACC_RDONLY);
	kp_init4matio_matrix(H_inf_tmp, mat2, "H_inf_8m");
	kalman_GPU.cu_H_inf = H_inf_tmp;
	Mat_Close(mat2);
	H_inf_tmp.resize(0,0);*/



	kp_multif_phasetur mf_phasetur;
	vector<string> fichiers_phase;

	fichiers_phase.push_back("/home/tgautrais/Documents/v2_kalman_kp/data/PHASTUR_08_160_VKv14Pix.mat");
	/*fichiers_phase.push_back("/home/tgautrais/Documents/v2_kalman_kp/data/PHASTUR_40_800_VKv14Pix_part1.mat");
	fichiers_phase.push_back("/home/tgautrais/Documents/v2_kalman_kp/data/PHASTUR_40_800_VKv14Pix_part2.mat");
	fichiers_phase.push_back("/home/tgautrais/Documents/v2_kalman_kp/data/PHASTUR_40_800_VKv14Pix_part3.mat");
	fichiers_phase.push_back("/home/tgautrais/Documents/v2_kalman_kp/data/PHASTUR_40_800_VKv14Pix_part4.mat");
	fichiers_phase.push_back("/home/tgautrais/Documents/v2_kalman_kp/data/PHASTUR_40_800_VKv14Pix_part5.mat");*/

	mf_phasetur.init(fichiers_phase);

	//utilise pour comparaison avec MATLAB
	stringstream nom_rand,nom_fichier;
	int boucle = 0;	


	delete PROJ_s ; PROJ_s = NULL;
	delete N_Act_s ; N_Act_s = NULL;
	delete D_Mo_s ; D_Mo_s = NULL;
	
	delete PROJ_f ; PROJ_f = NULL;
	delete N_Act_f ; N_Act_f = NULL;
	delete D_Mo_f ; D_Mo_f = NULL;
	if (afficher_temps) temps_boucle.start();
	
	if (comparaison_matlab) mat = Mat_Open("/home/tgautrais/Documents/v2_kalman_kp/data/rand0.04.mat",MAT_ACC_RDONLY);



	//for (boucle=0 ; boucle<NB_BOUCLE ; boucle++)
	while (mf_phasetur.set_next_phase(CPT_km1) && ( boucle < NB_BOUCLE || NB_BOUCLE <= 0) )
	{
		//PHASE RESIDUELLE ( ENTRE  les INSTANTS  K-1 & K )
		//CPR_km1 = CPT_km1 - CPC_km1

		CPR_km1 = CPT_km1;
		CPR_km1 -= CPC_km1;


		// VARIANCE PHASE RESIDUELLE ( INSTANT K-1 )
		// Var_PhRes = var(CPR_km1)
		Var_PhRes.d[boucle] = CPR_km1.var();


		// randn = randn(NB_P);
		kp_randn(randn, rng);

		

		// CALCUL MESURE  ( A l' INSTANT K )
		// Y_k = D_Sys * CPR_km1
		kp_gemv (1, D_Sys, CPR_km1, 0, Y_k);
		

		//pour comparer le resultat avec matlab, on lit la matrice rand au lieu de la generer
		if (comparaison_matlab)
		{
			nom_rand.str("");
			nom_rand << "rand" << boucle;
			kp_init4matio_vector(racSigmaW_randn,mat,nom_rand.str());
		}
		

		// Y_k = Y_k + randn * SQRT_BRUIT_PIX
		if (!comparaison_matlab)
		{
			randn *=  SQRT_BRUIT_PIX;
			Y_k += randn;
		}
		else Y_k += racSigmaW_randn;




		//kalman_GPU.temps_init_1boucle.reset();kalman_GPU.temps_1boucle.reset();
		temps_effectif.start();

		kalman_GPU.next_step(Y_k, U_k);
		
		temps_effectif.pause();


		
		// PHASE de CORRECTION
		// CPC_kp1 = N_Ph * U_k
		kp_gemv (1, N_Ph, U_k, 0, CPC_kp1);



		// CPC_kp1 = CPC_kp1 - mean(CPC_kp1)
		mean_CPCkp1 = CPC_kp1.mean();
		CPC_kp1 -= mean_CPCkp1; 
	
		//MISE A JOUR
		
		CPC_km1 = CPC_k; 
		CPC_k = CPC_kp1; 
		
		//if (boucle%200 == 0) cout << "boucle " << boucle <<endl;
	
		// ecriture d'un fichier pour comparer avec matlab
		/*if (comparaison_matlab)
		{		
			kp_cu2kp_vector(X_kskm1, kalman_GPU.cu_X_kskm1);
			kp_cu2kp_vector(Y_kskm1, kalman_GPU.cu_Y_kskm1);
			kp_cu2kp_vector(U_k, kalman_GPU.cu_U_k);
			fichier2.open("kalman.dat",ios::app);
			if(fichier2)
			{
				fichier2 << boucle << " " <<  __SP U_k.d[212] << " " << __SP CPR_km1.d[boucle] \
					<< " " << __SP CPC_k.d[boucle] << " " << __SP Y_kskm1.d[212] << " " <<\
				       	__SP X_kskm1.d[212] << " "<< __SP Var_PhRes.d[boucle]<< endl;
				fichier2.close();
			}
		}*/

		/*if (comparaison_matlab)
		{
		
			kp_cu2kp_vector(X_kp1sk, kalman_GPU.cu_X_kp1sk);
			kp_cu2kp_vector(Y_kskm1, kalman_GPU.cu_Y_kskm1);
			//kp_cu2kp_vector(U_k, kalman_GPU.cu_U_k);
			fichier2.open("U_kD.dat",ios::app);
			if(fichier2)
			{
				for (int i=0 ; i<U_k.size() ; i++)
					fichier2 << __SP U_k.d[i] << " ";
			}
			fichier2<<endl;
			fichier2.close();

			fichier2.open("Y_kskm1D.dat",ios::app);
			if(fichier2)
			{
				for (int i=0 ; i<Y_kskm1.size() ; i++)
					fichier2 << __SP Y_kskm1.d[i] << " ";
			}
			fichier2<<endl;
			fichier2.close();

			fichier2.open("X_kp1skD.dat",ios::app);
			if(fichier2)
			{
				for (int i=0 ; i<X_kp1sk.size() ; i++)
					fichier2 << __SP X_kp1sk.d[i] << " ";
			}
			fichier2<<endl;
			fichier2.close();


		}*/

		boucle++;
	}
	
		//cout<<"ratio="<<kalman_GPU.temps_init_1boucle.rez()<<" "<<kalman_GPU.temps_1boucle.rez()<<" "<<kalman_GPU.temps_init_1boucle.rez()/kalman_GPU.temps_1boucle.rez()<<endl;
	if (comparaison_matlab) Mat_Close(mat);


	
	
	cout << "---FIN EXECUTION KALMAN---" << endl;


	if (afficher_temps) 
	{
		temps_boucle.pause();
		temps_total.pause();
		cout << "Temps boucle = " << temps_boucle.rez() << " s" << endl;
	
		/*cout << "Temps op1 = " << kalman_GPU.temps_op1.rez() << " s" << endl;
		cout << "Temps op2 = " << kalman_GPU.temps_op2.rez() << " s" << endl;
		cout << "Temps op3 = " << kalman_GPU.temps_op3.rez() << " s" << endl;*/
		
		cout << "Temps effectif = " << temps_effectif.rez() << " s" << endl;
		cout << "Temps total = " << temps_total.rez() << endl;

		/*fichier.open("temps2.dat",ios_base::app);
		fichier << temps_ecoule(debut1,fin1).tv_sec<<"."<<temps_ecoule(debut1,fin1).tv_nsec << " " <<\
		temps_ecoule(debut2,fin2).tv_sec<<"."<< temps_ecoule(debut2,fin2).tv_nsec <<" "<<\
		temps_ecoule(debut3,fin3).tv_sec<<"."<< temps_ecoule(debut3,fin3).tv_nsec <<" "<< \
		temps5.tv_sec<<"."<< temps5.tv_nsec <<" "<<temps6.tv_sec<<"."<< temps6.tv_nsec <<" "<<\
		temps7.tv_sec<<"."<< temps7.tv_nsec <<endl;
		fichier.close();*/
	}
	kp_vector* VarPhRes_fin = new kp_vector(NB_BOUCLE-i0var+1) ;
	//fichier.open("Var_Ph_Res.dat",ios_base::app);
	for (i=0 ; i<NB_BOUCLE-i0var+1 ; i++)
	{
		VarPhRes_fin->d[i] = Var_PhRes.d[i+i0var-1];
	}
	cout <<"length_th="<<NB_BOUCLE-i0var+1<<endl;
	cout <<"length_reel=" << VarPhRes_fin->size()<<endl;
	//fichier.close();*/
	E_coh = exp(- VarPhRes_fin->mean())*100;
	
	cout << "E_coh = " << E_coh<<endl;
	cout << "mean_var = "<<	VarPhRes_fin->mean()<<endl;
	
	delete VarPhRes_fin ; VarPhRes_fin = NULL;


	
	cusparseDestroy(cusparseHandle);
	cublasDestroy(cublasHandle);
	magma_finalize();


	return 0;
}


	
//retourne 0 si OK
//retourne 1 si pb de lecture de fichier
//retourne 3 si probleme de multiplication matrice-vecteur
//retourne 4 si probleme d'addidtion de vecteurs
//retourne 5 si probleme de copie de vecteurs
//retourne 6 si probleme de soustraction de vecteurs


