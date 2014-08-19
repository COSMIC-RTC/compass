//sutra_controller_kalman.cpp

#include "sutra_controller_kalman.h"

#ifdef COMPILATION_LAM

#include "kp_kalman_core_sparse_GPU.h"
#include "kp_kalman_core_full_GPU.h"
#include "kp_kalman_core_sparse_CPU.h"
#include "kp_kalman_core_full_CPU.h"
#include "kp_smatrix.h"
#include "kp_carma_tools.h"
#include <fstream>

#include <cuda.h>
carma_obj<float>* calculate_D_Mo(carma_context* context, int n_slopes, int n_actu_zern, bool is_zonal) {
  long dims[] = { 2, n_slopes, n_actu_zern }; //877 pour cas particulier en 16m
  ifstream* vec_stream;
  if (is_zonal)
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/D_Mo_08P16_AR1Za980.dat");
  else
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/D_Mo_08P16_AR1Mo30.dat");
  istream_iterator<float> start(*vec_stream), eos;
  vector<float> vec(start, eos); 
  carma_obj<float>* D_Mo = new carma_obj<float>(context, dims, &vec[0]);
  //delete vec_stream;
  return D_Mo;
}
carma_obj<float>* calculate_N_Act(carma_context* context, int n_actu_zern, int n_actus, bool is_zonal) { 

  long dims[] = { 2, n_actu_zern, n_actus };
  ifstream* vec_stream;
  if (is_zonal)
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/N_Act_08P16_AR1Za980.dat");
  else
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/N_Act_08P16_AR1Mo30.dat");
  istream_iterator<float> start(*vec_stream), eos;
  vector<float> vec(start, eos); 
  carma_obj<float>* N_Act = new carma_obj<float>(context, dims, &vec[0]);  
  delete vec_stream;
  return N_Act;
}
carma_obj<float>* calculate_PROJ(carma_context* context, int n_actus, int n_actu_zern, bool is_zonal) {

  long dims[] = { 2, n_actus, n_actu_zern };
  ifstream* vec_stream;
  if (is_zonal)
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/PROJ_08P16_AR1Za980.dat");
  else
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/PROJ_08P16_AR1Mo30.dat");
  istream_iterator<float> start(*vec_stream), eos;
  vector<float> vec(start, eos); 
  carma_obj<float>* PROJ = new carma_obj<float>(context, dims, &vec[0]);
  delete vec_stream;
  return PROJ;
}
carma_obj<float>* calculate_btur(carma_context* context, int n_actu_zern, bool is_zonal) {
  /*long dims[] = { 1, n_actu_zern };
  float* zeros_tmp = new float[n_actu_zern];
  for (int i=0 ; i<n_actu_zern ; i++) zeros_tmp[i] = 0;
  carma_obj<float>* btur = new carma_obj<float>(context, dims, zeros_tmp);
  delete [] zeros_tmp ; zeros_tmp=NULL;*/
  return NULL;//btur;
}
carma_obj<float>* calculate_SigmaV(carma_context* context, int n_actu_zern, bool is_zonal) {

  long dims[] = { 2, n_actu_zern, n_actu_zern };
  ifstream* vec_stream;
  if (is_zonal)
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/SigmaV_08P16_AR1Za980.dat");
  else
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/SigmaV_08P16_AR1Mo30.dat");
  istream_iterator<float> start(*vec_stream), eos;
  vector<float> vec(start, eos);
 //for (int i=0 ; i<vec.size() ; i++) vec[i] = vec[i]/100; 
  carma_obj<float>* SigmaV = new carma_obj<float>(context, dims, &vec[0]);
  delete vec_stream;
  return SigmaV;
}
carma_obj<float>* calculate_atur(carma_context* context, int n_actu_zern, bool is_zonal) {
  /*long dims[] = { 1, n_actu_zern };
  ifstream* vec_stream;
  if (is_zonal)
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/atur_08P16_AR1Za980.dat");
  else
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/atur_08P16_AR1Mo30.dat");
  istream_iterator<float> start(*vec_stream), eos;
  vector<float> vec(start, eos); */
  carma_obj<float>* atur = NULL;//new carma_obj<float>(context, dims, &vec[0]);
  //delete vec_stream;
  return atur;
}


sutra_controller_kalman::sutra_controller_kalman(carma_context* context_, int nslope_, int nactu_) : sutra_controller(context_, nslope_, nactu_) {
   core_sparse = NULL;
   core_full = NULL;
   cusparseHandle = NULL;
   isGPU = true;
   isZonal = true;
   isSparse = true;
   isInit = false;
   isGainSet = false;

   matrices_matlab = false;
   pentes_matlab = false;
   sigmaVmatlab = false;
  
   /*int nvalid = nslope_/2;
   delay = 1;
   if (delay > 0) {
   long dims_data2[3];
    dims_data2[1] = 2 * nvalid;
    dims_data2[2] = delay + 1;
    d_cenbuff = new carma_obj<float>(this->current_context, dims_data2);
  }*/

   if(pentes_matlab)
   {
      Yk = new real* [5000];
      for (int i=0 ; i<5000 ; i++) Yk[i] = new real [416];

       ind_Yk=0;
       string line;
       int row,col;
       ifstream pFile ("Yk_8m.dat");
       if (pFile.is_open())
       {
           row=0;
           while(!pFile.eof())
           {
               getline(pFile, line);
               stringstream ss(line);
               col=0;
               while(ss >> Yk[row][col])
               {
                   col++;
               }
               row++;
           } 
           pFile.close();
       }
       else
       { 
           cerr << "Unable to open file"; 
	   exit(EXIT_FAILURE);
       }
   }

}



void sutra_controller_kalman::init_kalman(carma_host_obj<float>& chD_Mo, carma_host_obj<float>& chN_Act, carma_host_obj<float>& chPROJ, bool is_zonal, bool is_sparse, bool is_GPU) {
   
   if (isInit)
   {
           cerr << "Error |sutra_controller_kalman::init_kalman | Kalman controller has already been initialized"<<endl;
           exit(EXIT_FAILURE);
   }

   core_sparse = NULL;
   core_full = NULL;
   isGPU = is_GPU;
   isZonal = is_zonal;
   isSparse = is_sparse;


   //convert from carma_host_obj to kp_matrix
   kp_matrix kD_Mo, kN_Act, kPROJ;

   if (!matrices_matlab)
   {
      kp_carma_host_obj_to_kp_matrix(chD_Mo,  kD_Mo);
      kp_carma_host_obj_to_kp_matrix(chN_Act, kN_Act);
      kp_carma_host_obj_to_kp_matrix(chPROJ,  kPROJ);
     
      /*ofstream fichier;
      fichier.open("Matrice_D_Mo_plateforme.dat",ios::out);
      for (int i=0 ; i < kD_Mo.dim1 ; i++)
      {
         for (int j=0 ; j<kD_Mo.dim2 ; j++)
         {
            fichier << kD_Mo(i,j) << " ";
         }
         fichier << endl;
      }
      fichier.close();*/

   }
   else
   {
      //Utilisation des matrices MATLAB
      carma_obj<float> *cN_Act, *cPROJ, *cD_Mo;
      cN_Act =  calculate_N_Act(current_context, nactu(), nactu(), isZonal);
      cPROJ = calculate_PROJ(current_context, nactu(), nactu(), isZonal);
      cD_Mo = calculate_D_Mo(current_context, nslope(), nactu(), isZonal);
      kp_carma_obj_to_kp_matrix(*cN_Act, kN_Act);
      kp_carma_obj_to_kp_matrix(*cPROJ,  kPROJ);
      kp_carma_obj_to_kp_matrix(*cD_Mo,  kD_Mo);
   }



   cudaSetDevice(2);

   if (is_sparse)
   {
      //convert from kp_matrix to kp_smatrix (full -> sparse)
      kp_smatrix sD_Mo, sN_Act, sPROJ;
      sD_Mo.init_from_matrix(kD_Mo);sD_Mo.resize2rowMajor();
      sN_Act.init_from_matrix(kN_Act);sN_Act.resize2rowMajor();
      sPROJ.init_from_matrix(kPROJ);sPROJ.resize2rowMajor();
      
      if (is_GPU)
      {
         cusparseStatus_t cusparseStat = cusparseCreate(&cusparseHandle);
         if (cusparseStat != CUSPARSE_STATUS_SUCCESS)
         { 
            cerr<<"Error | sutra_controller_kalman::sutra_controller_kalman  | cusparseCreate failed "<<endl;
            exit(EXIT_FAILURE);
         }
         core_sparse = new kp_kalman_core_sparse_GPU(sD_Mo, sN_Act, sPROJ,
					       is_zonal, 
					       current_context->get_cublasHandle(), 
					       cusparseHandle);
      }
      else
         core_sparse = new kp_kalman_core_sparse_CPU(sD_Mo, sN_Act, sPROJ,
					       is_zonal);

   }  
   else if (is_GPU)
      core_full = new kp_kalman_core_full_GPU(kD_Mo, kN_Act, kPROJ, is_zonal, current_context->get_cublasHandle());
   else
      core_full = new kp_kalman_core_full_CPU(kD_Mo, kN_Act, kPROJ, is_zonal);

   cudaSetDevice(0);
   
   isInit = true;
}

sutra_controller_kalman::~sutra_controller_kalman() {
   if(core_full)
   {
      delete core_full;
      core_full=NULL;
   }
   if(core_sparse)
   {
      delete core_sparse;
      core_sparse=NULL;
   }

   if (cusparseHandle)
   {
      cusparseDestroy(cusparseHandle);
      cusparseHandle = NULL;
   }
   if (pentes_matlab)
   {
      for (int i=0;i<5000 ; i++) delete [] Yk[i];
      delete[] Yk;
   }
   
   //if (delay > 0) delete d_cenbuff;

}

void sutra_controller_kalman::calculate_gain(double bruit,
    carma_host_obj<float>& chSigmaV, carma_host_obj<float>& chatur,
    carma_host_obj<float>& chbtur) {
	
   if (!isInit)
   {
      cerr << "Error | sutra_controller_kalman::calculate_gain | Kalman controller has not been initialiez"<<endl;
      exit(EXIT_FAILURE);
   }
   //convert carma_obj to kp_matrix
   kp_matrix kSigmaV;
   if (!sigmaVmatlab)
      kp_carma_host_obj_to_kp_matrix(chSigmaV, kSigmaV);
   else 
   {
      //Utilisation tmp de SigmaV Matlab
      carma_obj<float>* chtestSigmaV;
      if (isZonal)
          chtestSigmaV = calculate_SigmaV(current_context, nactu(), isZonal);
      else
          chtestSigmaV = calculate_SigmaV(current_context, 495, isZonal);
      kp_carma_obj_to_kp_matrix(*chtestSigmaV, kSigmaV); 
   }

   //cout<<"SigmaV : "<<kSigmaV.dim1<<"x"<<kSigmaV.dim2<<endl;
   /*ofstream fichier;
   fichier.open("Matrice_test_SigmaV.dat",ios::out);
   for (int i=0 ; i < kSigmaV.dim1 ; i++)
   {
      for (int j=0 ; j<kSigmaV.dim2 ; j++)
      {
         fichier << kSigmaV(i,j) << " ";
      }
      fichier << endl;
   }
   fichier.close();*/

   
   //convert carma_obj to kp_vector
   kp_vector katur, kbtur;
   kp_carma_host_obj_to_kp_vector(chatur, katur);
   kp_carma_host_obj_to_kp_vector(chbtur, kbtur);

   cudaSetDevice(2);
   //gain (attribut de la classe) correspond a k_W
   if (!isGainSet)
   {
      cerr << "k_W (gain parameter) has not been set" << endl;
      exit(EXIT_FAILURE);
   }
   if (core_sparse)
      core_sparse->calculate_gain(bruit, gain, kSigmaV, katur, kbtur);
   else if (core_full)
      core_full->calculate_gain(bruit, gain, kSigmaV, katur, kbtur);
   cudaSetDevice(0);
}
//                                                              
/*int sutra_controller_kalman::frame_delay() {
  // here we place the content of d_centroids into cenbuf and get
  // the actual centroid frame for error computation depending on delay value

  if (delay > 0) {
    for (int cc = 0; cc < delay; cc++)
      shift_buf(&((this->d_cenbuff->getData())[cc * this->nslope()]), 1,
          this->nslope(), this->device);

    cutilSafeCall(
        cudaMemcpy(&(this->d_cenbuff->getData()[delay * this->nslope()]),
            this->d_centroids->getData(), sizeof(float) * this->nslope(),
            cudaMemcpyDeviceToDevice));
    
    cutilSafeCall(
        cudaMemcpy(this->d_centroids->getData(), this->d_cenbuff->getData(),
            sizeof(float) * this->nslope(), cudaMemcpyDeviceToDevice));
  }

  return EXIT_SUCCESS;
}*/                       
int sutra_controller_kalman::comp_com() {
   //frame_delay();
   //d_com2->copy(this->d_com1, 1, 1);
   //d_com1->copy(this->d_com, 1, 1);
   
   kp_vector Y_k,Y_k_tmp, U_k;

   //d_centroids->axpy(1/0.3, d_centroids, 1, 1);

   //conversion des pentes d'arcsec en px : Y_k[px] = Y_k[arcsec]/pixsize[arcsec/px]
   //carma_obj<float> cd_centroids_tmp(current_context, d_centroids->getDims());
   //cd_centroids_tmp.axpy(1/0.135, d_centroids, 1, 1);  // 0.135 = pixsize :  angle en arcsec par rapport a l'axe optique de la microlentille permettant d'obtenir un deplacement du spot sur detecteur de 1 pixel

   kp_carma_obj_to_kp_vector(*d_centroids, Y_k); 

   ofstream fichier;
   fichier.open("Yk_auto.dat",ofstream::app);
   for (int i=0 ; i<Y_k.size() ; i++)
      fichier << __SP Y_k.d[i]<<" ";
   fichier<<endl;
   fichier.close();

   if (pentes_matlab) Y_k.d = Yk[ind_Yk];//SUPPR

   cudaSetDevice(2);
   if (core_sparse)
   {
      core_sparse->next_step(Y_k, U_k);
   }
   else if (core_full)
   {
      core_full->next_step(Y_k, U_k);
   }
   cudaSetDevice(0);

   //conversion des tensions de rad en V : U_k[V] = U_k[rad]*lamda[um]/(2*pi)/unitpervolt[um/V] 
   //U_k *= 1.654/(2*M_PI)/0.01;
   //U_k *= 1.654/(2*M_PI)/100000;
   //U_k *= -13.3802;
  U_k *= -1; 

   fichier.open("Uk_auto.dat",ofstream::app);
   for (int i=0 ; i<U_k.size() ; i++)
      fichier << __SP U_k.d[i]<<" ";
   fichier<<endl;
   fichier.close();

 
   //U_k.zeros();
   kp_kp_vector_to_carma_obj(U_k, *d_com);


   if (pentes_matlab) ind_Yk++;//SUPPR
  return -378;
}

#else

carma_obj<float>* calculate_D_Mo(carma_context* context, int n_slopes, int n_actu_zern, bool is_zonal) {
  return NULL;
}
carma_obj<float>* calculate_N_Act(carma_context* context, int n_actu_zern, int n_actus, bool is_zonal) {
  return NULL;
}
carma_obj<float>* calculate_PROJ(carma_context* context, int n_actus, int n_actu_zern, bool is_zonal) {
  return NULL;
}
carma_obj<float>* calculate_btur(carma_context* context, int n_actu_zern, bool is_zonal) {
  return NULL;
}
carma_obj<float>* calculate_SigmaV(carma_context* context, int n_actu_zern, bool is_zonal) {
  return NULL;
}
carma_obj<float>* calculate_atur(carma_context* context, int n_actu_zern, bool is_zonal) {
  return NULL;
}

sutra_controller_kalman::sutra_controller_kalman(carma_context* context_, int nslope_, int nactu_) : sutra_controller(context_, nslope_, nactu_) {
   core_sparse = NULL;
   core_full = NULL;
   cusparseHandle = NULL;
   isGPU = true;
   isZonal = true;
   isSparse = true;
   isInit = false;
}

void sutra_controller_kalman::init_kalman(carma_host_obj<float>& chD_Mo, 
		carma_host_obj<float>& chN_Act, carma_host_obj<float>& chPROJ, 
		bool is_zonal, bool is_sparse, bool is_GPU) {
}
sutra_controller_kalman::~sutra_controller_kalman() {
}

void sutra_controller_kalman::calculate_gain(double bruit,
    carma_host_obj<float>& chSigmaV, carma_host_obj<float>& chatur,
    carma_host_obj<float>& chbtur) {
}
//
int sutra_controller_kalman::comp_com() {
  return -378;
}

//int sutra_controller_kalman::frame_delay() {}
#endif

int sutra_controller_kalman::set_gain(float k_W) {
  this->gain = k_W;
   isGainSet = true;
  return EXIT_SUCCESS;
}


//                                                                                     
