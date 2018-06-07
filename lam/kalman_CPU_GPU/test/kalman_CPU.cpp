#include "kalman_CPU.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>

#define __SP setprecision(20) <<

int main(int argc, char* argv[]) {
  const bool afficher_temps = 1;
  const bool comparaison_matlab = 0;
  const bool sparse_matrix = 1;
  bool isZonal = true;

  gsl_rng* rng;
  int i;
  KFPP E_coh = 0;

  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng, time(NULL));

  kp_timer temps_total, temps_riccati, temps_boucle, temps_Ykskm1_1,
      temps_Ykskm1_2, temps_Xkp1sk_1, temps_Xkp1sk_2, temps_Uk, temps_effectif;
  if (afficher_temps) temps_total.start();

  cout << ("---DEBUT EXECUTION KALMAN---") << endl;

  ofstream fichier, fichier2;
  mat_t* mat;

  // CHARGEMENT DES MATRICES (load de Matlab)
  kp_smatrix<KFPP>* D_Mo_s = new kp_smatrix<KFPP>();
  kp_smatrix<KFPP>* N_Act_s = new kp_smatrix<KFPP>();
  kp_smatrix<KFPP>* PROJ_s = new kp_smatrix<KFPP>();
  kp_matrix<KFPP>* D_Mo_f = new kp_matrix<KFPP>();
  kp_matrix<KFPP>* N_Act_f = new kp_matrix<KFPP>();
  kp_matrix<KFPP>* PROJ_f = new kp_matrix<KFPP>();

  kp_smatrix<KFPP> D_Sys;
  kp_smatrix<double>* A1 = new kp_smatrix<double>();
  kp_smatrix<KFPP> N_Ph;

  kp_matrix<double> SigmaV;

  int nb_z = -1;
  int nb_az = -1;

  if (isZonal)
    mat = Mat_Open(
        "/home/tgautrais/Documents/v2_kalman_kp/data/KALM_08P16_AR1Za980.mat",
        MAT_ACC_RDONLY);
  else
    mat = Mat_Open(
        "/home/tgautrais/Documents/v2_kalman_kp/data/KALM_08P16_AR1Mo30.mat",
        MAT_ACC_RDONLY);
  // LECTURE DE D_Sys
  kp_init4matio_smatrix(D_Sys, mat, "D_Sys");
  // LECTURE DE A1
  kp_init4matio_smatrix(*A1, mat, "A1");
  // LECTURE DE SigmaV
  kp_init4matio_matrix(SigmaV, mat, "SigmaV");
  // LECTURE DE D_Mo, N_Act et PROJ
  if (isZonal) {
    kp_init4matio_smatrix(*D_Mo_s, mat, "D_Mo");
    kp_init4matio_smatrix(*N_Act_s, mat, "N_Act");
    kp_init4matio_smatrix(*PROJ_s, mat, "PROJ");
  } else {
    kp_init4matio_matrix(*D_Mo_f, mat, "D_Mo");
    kp_init4matio_matrix(*N_Act_f, mat, "N_Act");
    kp_init4matio_matrix(*PROJ_f, mat, "PROJ");
  }
  // LECTURE DE N_Ph
  kp_init4matio_smatrix(N_Ph, mat, "N_Ph");
  // LECTURE DES DIMENSIONS
  const int NB_ACT = kp_init4matio_int(mat, "nb_act");
  const int NB_P = kp_init4matio_int(mat, "nb_p");
  const int NB_PX = kp_init4matio_int(mat, "nb_px");
  if (!isZonal) {
    nb_z = kp_init4matio_int(mat, "nb_z");
    nb_az = nb_z;
  } else
    nb_az = NB_ACT;
  Mat_Close(mat);
  const int NB_AZ = nb_az;
  // INITIALISATIONS et PARAMETRES
  const int i0var = 51;
  const float BRUIT_RAD = 0.04;

  const float k_W = 5;
  const float BRUIT_PIX = BRUIT_RAD / (M_PI * M_PI);
  const float SQRT_BRUIT_PIX = sqrt(BRUIT_PIX);
  const int NB_BOUCLE = 5000;

  if (!sparse_matrix && isZonal) {
    D_Mo_f->init_from_smatrix(*D_Mo_s);
    N_Act_f->init_from_smatrix(*N_Act_s);
    PROJ_f->init_from_smatrix(*PROJ_s);
  }
  if (sparse_matrix && !isZonal) {
    D_Mo_s->init_from_matrix(*D_Mo_f);
    N_Act_s->init_from_matrix(*N_Act_f);
    PROJ_s->init_from_matrix(*PROJ_f);
  }

  if (sparse_matrix) {
    D_Mo_s->resize2rowMajor();
    N_Act_s->resize2rowMajor();
    PROJ_s->resize2rowMajor();
  }

  kp_kalman_core_sparse_CPU kalman_CPU(*D_Mo_s, *N_Act_s, *PROJ_s, isZonal);
  // kp_kalman_core_full_CPU kalman_CPU(*D_Mo_f, *N_Act_f, *PROJ_f, isZonal);

  // Atur
  kp_vector<double> A1_00(NB_AZ);
  A1_00.zeros();
  // Btur (different de zero si AR2)
  kp_vector<double> A1_01(NB_AZ);
  A1_01.zeros();
  A1_2vec(A1_00, A1_01, *A1);

  /*fichier.open("A1_00.dat",ios::out);
  fichier2.open("A1_01.dat",ios::out);
  for (i=0;i<NB_AZ;i++)
  {
          fichier<< A1_00->d[i]<<" ";
          fichier2<< A1_01->d[i]<<" ";
  }
  fichier.close();
  fichier2.close();*/

  int s = 0;

  // Initialisation boucle OA
  kp_vector<KFPP> Var_PhRes(NB_BOUCLE);
  Var_PhRes.zeros();
  kp_vector<KFPP> CPC_km1(NB_PX);
  CPC_km1.zeros();
  kp_vector<KFPP> CPC_k(NB_PX);
  CPC_k.zeros();
  kp_vector<KFPP> CPC_kp1(NB_PX);
  CPC_kp1.zeros();
  kp_vector<KFPP> CPT_km1(NB_PX);
  kp_vector<KFPP> CPR_km1(NB_PX);
  kp_vector<KFPP> Y_k(NB_P);
  kp_vector<KFPP> randn(NB_P);
  kp_vector<KFPP> racSigmaW_randn(NB_P);
  kp_vector<KFPP> U_k(NB_ACT);
  U_k.zeros();

  // kp_vector<KFPP>* X_kskm1_fin = new kp_vector<KFPP>(NB_AZ);

  KFPP mean_CPCkp1;

  if (afficher_temps) temps_riccati.start();

  kalman_CPU.calculate_gain(BRUIT_PIX, k_W, SigmaV, A1_00, A1_01);

  if (afficher_temps) {
    temps_riccati.pause();
    cout << "      Temps riccati = " << temps_riccati.rez() << endl;
  }

  kp_multif_phasetur<KFPP> mf_phasetur;
  vector<string> fichiers_phase;

  fichiers_phase.push_back(
      "/home/tgautrais/Documents/v2_kalman_kp/data/"
      "PHASTUR_08_160_VKv14Pix.mat");
  /*fichiers_phase.push_back("/home/tgautrais/Documents/v2_kalman_kp/data/PHASTUR_40_800_VKv14Pix_part1.mat");
  fichiers_phase.push_back("/home/tgautrais/Documents/v2_kalman_kp/data/PHASTUR_40_800_VKv14Pix_part2.mat");
  fichiers_phase.push_back("/home/tgautrais/Documents/v2_kalman_kp/data/PHASTUR_40_800_VKv14Pix_part3.mat");
  fichiers_phase.push_back("/home/tgautrais/Documents/v2_kalman_kp/data/PHASTUR_40_800_VKv14Pix_part4.mat");
  fichiers_phase.push_back("/home/tgautrais/Documents/v2_kalman_kp/data/PHASTUR_40_800_VKv14Pix_part5.mat");*/

  mf_phasetur.init(fichiers_phase);

  if (s != 0) return 2;

  // utilise pour comparaison avec MATLAB
  stringstream nom_rand, nom_fichier;
  int boucle = 0;

  delete PROJ_s;
  PROJ_s = NULL;
  delete N_Act_s;
  N_Act_s = NULL;
  delete D_Mo_s;
  D_Mo_s = NULL;

  delete PROJ_f;
  PROJ_f = NULL;
  delete N_Act_f;
  N_Act_f = NULL;
  delete D_Mo_f;
  D_Mo_f = NULL;

  if (afficher_temps) temps_boucle.start();

  if (comparaison_matlab)
    mat = Mat_Open("/home/tgautrais/Documents/v2_kalman_kp/data/rand0.04.mat",
                   MAT_ACC_RDONLY);

  // for (boucle=0 ; boucle<NB_BOUCLE ; boucle++)
  while (mf_phasetur.set_next_phase(CPT_km1) &&
         (boucle < NB_BOUCLE || NB_BOUCLE <= 0)) {
    // PHASE RESIDUELLE ( ENTRE  les INSTANTS  K-1 & K )
    // CPR_km1 = CPT_km1 - CPC_km1

    CPR_km1 = CPT_km1;
    CPR_km1 += CPC_km1;

    // VARIANCE PHASE RESIDUELLE ( INSTANT K-1 )
    // Var_PhRes = var(CPR_km1)
    *(Var_PhRes.getData(boucle)) = CPR_km1.var();

    // randn = randn(NB_P);
    kp_randn(randn, rng);

    // CALCUL MESURE  ( A l' INSTANT K )
    // Y_k = D_Sys * CPR_km1
    kp_gemv(1, D_Sys, CPR_km1, 0, Y_k);

    // pour comparer le resultat avec matlab, on lit la matrice rand au lieu de
    // la generer
    if (comparaison_matlab) {
      nom_rand.str("");
      nom_rand << "rand" << boucle;
      kp_init4matio_vector(racSigmaW_randn, mat, nom_rand.str());
    }

    // Y_k = Y_k + randn * SQRT_BRUIT_PIX
    if (!comparaison_matlab) {
      randn *= SQRT_BRUIT_PIX;
      Y_k += randn;
    } else
      Y_k += racSigmaW_randn;

    if (afficher_temps) temps_effectif.start();

    kalman_CPU.next_step(Y_k, U_k);

    if (afficher_temps) temps_effectif.pause();

    // PHASE de CORRECTION
    // CPC_kp1 = N_Ph * U_k
    kp_gemv(1, N_Ph, U_k, 0, CPC_kp1);

    // CPC_kp1 = CPC_kp1 - mean(CPC_kp1)
    mean_CPCkp1 = CPC_kp1.mean();
    CPC_kp1 -= mean_CPCkp1;

    // MISE A JOUR
    CPC_km1 = CPC_k;
    CPC_k = CPC_kp1;

    // if (boucle%200 == 0) cout << "boucle " << boucle <<endl;

    // ecriture d'un fichier pour comparer avec matlab
    /*if (comparaison_matlab)
    {
            fichier2.open("kalman.dat",ios::app);
            if(fichier2)
            {
                    fichier2 <<
                            boucle << " " <<  __SP
                            kalman_CPU.U_km1.d[212] << " " << __SP
                            CPR_km1.d[boucle] << " " << __SP
                            CPC_k.d[boucle] << " " << __SP
                            kalman_CPU.Y_kskm1.d[212] << " " << __SP
                            kalman_CPU.X_kskm1.d[212] << " "<< __SP
                            Var_PhRes.d[boucle]<< endl;
                    fichier2.close();
            }
    }*/

    /*if (comparaison_matlab)
    {
            fichier2.open("Uk_40mDS.dat",ios_base::app);
            for (int i=0 ; i<NB_ACT; i++) fichier2 << __SP U_k.d[i] << " ";
            fichier2 << endl;
            fichier2.close();
    }*/

    // if (boucle == 2500)
    //      kalman_CPU.calculate_gain(BRUIT_PIX, k_W, SigmaV, A1_00, A1_01);

    boucle++;
  }

  if (comparaison_matlab) Mat_Close(mat);

  cout << "---FIN EXECUTION KALMAN---" << endl;

  if (afficher_temps) {
    temps_boucle.pause();
    temps_total.pause();
    cout << "Temps boucle = " << temps_boucle.rez() << " s" << endl;

    /*cout << "Temps op1 = " << kalman_CPU.temps_op1.rez() << " s" << endl;
    cout << "Temps op2 = " << kalman_CPU.temps_op2.rez() << " s" << endl;
    cout << "Temps op3 = " << kalman_CPU.temps_op3.rez() << " s" << endl;*/

    cout << "Temps total = " << temps_total.rez() << " s" << endl;
    cout << "Temps effectif = " << temps_effectif.rez() << " s" << endl;

    /*fichier.open("temps2.dat",ios_base::app);
    fichier <<
    temps_ecoule(debut1,fin1).tv_sec<<"."<<temps_ecoule(debut1,fin1).tv_nsec <<
    " " <<\
    temps_ecoule(debut2,fin2).tv_sec<<"."<< temps_ecoule(debut2,fin2).tv_nsec
    <<" "<<\
    temps_ecoule(debut3,fin3).tv_sec<<"."<< temps_ecoule(debut3,fin3).tv_nsec
    <<" "<< \
    temps5.tv_sec<<"."<< temps5.tv_nsec <<" "<<temps6.tv_sec<<"."<<
    temps6.tv_nsec <<" "<<\ temps7.tv_sec<<"."<< temps7.tv_nsec <<endl;
    fichier.close();*/
  }

  kp_vector<KFPP>* VarPhRes_fin = new kp_vector<KFPP>(NB_BOUCLE - i0var + 1);
  // fichier.open("Var_Ph_Res.dat",ios_base::app);
  for (i = 0; i < NB_BOUCLE - i0var + 1; i++) {
    *(VarPhRes_fin->getData(i)) = Var_PhRes[i + i0var - 1];
    // fichier << VarPhRes_fin->d[i] <<endl;
  }
  cout << "length_th=" << NB_BOUCLE - i0var + 1 << endl;
  cout << "length_reel=" << VarPhRes_fin->size() << endl;
  // fichier.close();*/
  E_coh = exp(-VarPhRes_fin->mean()) * 100;

  cout << "E_coh = " << E_coh << endl;
  cout << "mean_var = " << VarPhRes_fin->mean() << endl;

  delete (VarPhRes_fin);

  return 0;
}

// retourne 0 si OK
// retourne 1 si pb de lecture de fichier
// retourne 3 si probleme de multiplication matrice-vecteur
// retourne 4 si probleme d'addidtion de vecteurs
// retourne 5 si probleme de copie de vecteurs
// retourne 6 si probleme de soustraction de vecteurs
