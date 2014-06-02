#include  <stdio.h>
#include  <stdlib.h>
#include  <sys/types.h>
#include  <unistd.h>

#include <mpi.h>

#include<cuda_runtime.h>
//#include"CudaTimer.h"

#include <carma_context.h>
#include <sutra_wfs_mpi.h>

#include <cfitsio/fitsio2.h>

#include "svipc.h"

template <typename T>
int read_shm(int shm_key, int *dims, slot_type YSVIPC_TYPE, char *id, T* data){
  long nb_elements=1;
  slot_array arr;
  arr.countdims = dims[0];
  arr.number = (int *) malloc(arr.countdims * sizeof(*arr.number));
  for(int i=0; i<dims[0]-1; i++){
    arr.number[i] = dims[i+1];
    nb_elements*=dims[i+1];
  }
  arr.typeID = YSVIPC_TYPE;

  arr.data = data;

  int status = svipc_shm_read(shm_key, id, &arr, 0);

  free(arr.number);
  status = svipc_shm_free(shm_key, id);
  return status;
}

int main(int argc, char *argv[]) {
  //int provided;
  //MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  //printf ("provided %d\n",provided);
  MPI_Init(&argc, &argv);

  int my_rank;
  int size;
  /* DETERMINE RANK OF THIS PROCESSOR*/
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /*DETERMINE TOTAL NUMBER OF PROCESSORS*/
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  long *d_data;
  cudaMalloc((void**)&d_data, 10000*sizeof(long));

  int shm_key = 1234;
  int sem_key = 2345;
  int status;
  if (my_rank == 0) {
    int i = 0;

    status = svipc_shm_init(shm_key, 32);
    status = svipc_sem_init(sem_key, 10);

    char* arg_list[] = { "yorick", /* argv[0], le nom du programme. */
    "-batch", "/home/sevin/compass/wrapper_mpi/ao_rtc-main.i",
    "/home/sevin/compass/yoga_ao/data/par/1wfs8x8_1layer_rtc_dm.par", NULL /* La liste d'arguments doit se terminer par NULL.  */
    };

    // intercommunicator
    MPI_Comm everyone;
    MPI_Comm_spawn("mpi_yorick", arg_list, 1, MPI_INFO_NULL, 0, MPI_COMM_SELF,
        &everyone, MPI_ERRCODES_IGNORE);
  }

  long nwfs;
  if (my_rank == 0) {
    printf("on attend Yorick [1]\n");
    status = svipc_semtake(sem_key, 0, 1, -1);

    char id_nwfs[] = "nwfs";
    int dims_wfs[]={1,1};
    status=read_shm<long>(shm_key, dims_wfs,SVIPC_LONG, id_nwfs, &nwfs);
  }

  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&nwfs, 1, MPI_LONG, 0, MPI_COMM_WORLD);

//  char type[3*nwfs];
  long nbsub[nwfs];
  long nvalid[nwfs];
  long npix[nwfs];
  long nphase[nwfs];
  long nrebin[nwfs];
  long ntot[nwfs];
  long npup[nwfs];
  long nfft[nwfs];
  float pdiam[nwfs];
  float nphot[nwfs];
  int lgs[nwfs];
  if (my_rank == 0) {
//    char id_type[] = "type";
//    int dims_type[]={1,128};
//    status=read_shm<char>(shm_key, dims_type,SVIPC_CHAR, id_type, type);

    char id_nbsub[] = "nbsub";
    int dims_nbsub[]={1,nwfs};
    status=read_shm<long>(shm_key, dims_nbsub,SVIPC_LONG, id_nbsub, nbsub);

    char id_nvalid[] = "nvalid";
    int dims_nvalid[]={1,nwfs};
    status=read_shm<long>(shm_key, dims_nvalid,SVIPC_LONG, id_nvalid, nvalid);

    char id_npix[] = "npix";
    int dims_npix[]={1,nwfs};
    status=read_shm<long>(shm_key, dims_npix,SVIPC_LONG, id_npix, npix);

    char id_nphase[] = "nphase";
    int dims_nphase[]={1,nwfs};
    status=read_shm<long>(shm_key, dims_nphase,SVIPC_LONG, id_nphase, nphase);

    char id_nrebin[] = "nrebin";
    int dims_nrebin[]={1,nwfs};
    status=read_shm<long>(shm_key, dims_nrebin,SVIPC_LONG, id_nrebin, nrebin);

    char id_nfft[] = "nfft";
    int dims_nfft[]={1,nwfs};
    status=read_shm<long>(shm_key, dims_nfft,SVIPC_LONG, id_nfft, nfft);

    char id_ntot[] = "ntot";
    int dims_ntot[]={1,nwfs};
    status=read_shm<long>(shm_key, dims_ntot,SVIPC_LONG, id_ntot, ntot);

    char id_npup[] = "npup";
    int dims_npup[]={1,nwfs};
    status=read_shm<long>(shm_key, dims_npup,SVIPC_LONG, id_npup, npup);

    char id_pdiam[] = "pdiam";
    int dims_pdiam[]={1,nwfs};
    status=read_shm<float>(shm_key, dims_pdiam,SVIPC_FLOAT, id_pdiam, pdiam);

    char id_nphot[] = "nphot";
    int dims_nphot[]={1,nwfs};
    status=read_shm<float>(shm_key, dims_nphot,SVIPC_FLOAT, id_nphot, nphot);

    char id_lgs[] = "lgs";
    int dims_lgs[]={1,nwfs};
    status=read_shm<int>(shm_key, dims_lgs,SVIPC_INT, id_lgs, lgs);

//    fprintf(stderr,"read SHM -> [%d] %d %d %d %d %d %d %d %d %d %f %f %d\n",
//        my_rank, nwfs, nbsub[0], nvalid[0], npix[0],
//        nphase[0], nrebin[0], ntot[0],
//        npup[0], nfft[0], pdiam[0],
//        nphot[0], lgs[0]);

  }

  //MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(nbsub, nwfs, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(nvalid, nwfs, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(npix, nwfs, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(nphase, nwfs, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(nrebin, nwfs, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(ntot, nwfs, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(npup, nwfs, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(nfft, nwfs, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(pdiam, nwfs, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(nphot, nwfs, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(lgs, nwfs, MPI_INT, 0, MPI_COMM_WORLD);

//  for(int i=0; i<size; i++){
//    if(i==my_rank){
//      fprintf(stderr,"[%d] %d %d %d %d %d %d %d %d %d %f %f %d\n",
//          my_rank, nwfs, nbsub[0], nvalid[0], npix[0],
//          nphase[0], nrebin[0], ntot[0],
//          npup[0], nfft[0], pdiam[0],
//          nphot[0], lgs[0]);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//  }

  carma_context current_context;
  int device = my_rank&1;
  current_context.set_activeDevice(device, 1);
  const char type[] = "sh";
  sutra_sensors_mpi sensors(&current_context, type, nwfs, nbsub, nvalid, npix,
      nphase, nrebin, nfft, ntot, npup[0], pdiam, nphot, lgs, device);

  float xpos[nwfs];
  float ypos[nwfs];
  float lambda[nwfs];
  float gsmag[nwfs];
  long gssize[nwfs];
  float noise[nwfs];
  if (my_rank == 0) {
    printf("on attend Yorick [2]\n");
    status = svipc_semtake(sem_key, 0, 1, -1);
     char id_xpos[] = "xpos";
     int dims_xpos[]={1,nwfs};
     status=read_shm<float>(shm_key, dims_xpos,SVIPC_FLOAT, id_xpos, xpos);

     char id_ypos[] = "ypos";
     int dims_ypos[]={1,nwfs};
     status=read_shm<float>(shm_key, dims_ypos,SVIPC_FLOAT, id_ypos, ypos);

     char id_lambda[] = "lambda";
     int dims_lambda[]={1,nwfs};
     status=read_shm<float>(shm_key, dims_lambda,SVIPC_FLOAT, id_lambda, lambda);

     char id_gsmag[] = "gsmag";
     int dims_gsmag[]={1,nwfs};
     status=read_shm<float>(shm_key, dims_gsmag,SVIPC_FLOAT, id_gsmag, gsmag);

     char id_gssize[] = "size";
     int dims_gssize[]={1,nwfs};
     status=read_shm<long>(shm_key, dims_gssize,SVIPC_LONG, id_gssize, gssize);

     char id_noise[] = "noise";
     int dims_noise[]={1,nwfs};
     status=read_shm<float>(shm_key, dims_noise,SVIPC_FLOAT, id_noise, noise);

  }

  MPI_Bcast(xpos, nwfs, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(ypos, nwfs, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(lambda, nwfs, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(gsmag, nwfs, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(gssize, nwfs, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(noise, nwfs, MPI_FLOAT, 0, MPI_COMM_WORLD);

  sensors.sensors_initgs(xpos, ypos, lambda, gsmag, gssize); //, noise);

  int phasemap[sensors.d_wfs[0]->d_phasemap->getNbElem()];
  int hrmap_size=1;
  if (sensors.d_wfs[0]->ntot != sensors.d_wfs[0]->nfft)
    hrmap_size=sensors.d_wfs[0]->d_hrmap->getNbElem();
  int hrmap[hrmap_size];
  int binmap[sensors.d_wfs[0]->d_binmap->getNbElem()];
  float offsets[sensors.d_wfs[0]->d_offsets->getNbElem()];
  float pupil[sensors.d_wfs[0]->d_pupil->getNbElem()];
  float fluxPerSub[sensors.d_wfs[0]->d_fluxPerSub->getNbElem()];
  int isvalid[sensors.d_wfs[0]->d_isvalid->getNbElem()];
  int validsubsx[sensors.d_wfs[0]->d_validsubsx->getNbElem()];
  int validsubsy[sensors.d_wfs[0]->d_validsubsy->getNbElem()];
  int istart[sensors.d_wfs[0]->d_istart->getNbElem()];
  int jstart[sensors.d_wfs[0]->d_jstart->getNbElem()];
  float ftkernel[2*sensors.d_wfs[0]->d_ftkernel->getNbElem()];

  if (my_rank == 0) {
    printf("on attend Yorick [3]\n");
    status = svipc_semtake(sem_key, 0, 1, -1);
    char id_phasemap[] = "phasemap";
    int dims_phasemap[]={1,sensors.d_wfs[0]->d_phasemap->getNbElem()};
    status=read_shm<int>(shm_key, dims_phasemap,SVIPC_INT, id_phasemap, phasemap);

    char id_hrmap[] = "hrmap";
    int dims_hrmap[]={1,hrmap_size};
    status=read_shm<int>(shm_key, dims_hrmap,SVIPC_INT, id_hrmap, hrmap);

    char id_binmap[] = "binmap";
    int dims_binmap[]={1,sensors.d_wfs[0]->d_binmap->getNbElem()};
    status=read_shm<int>(shm_key, dims_binmap,SVIPC_INT, id_binmap, binmap);

    char id_offsets[] = "offsets";
    int dims_offsets[]={1,sensors.d_wfs[0]->d_offsets->getNbElem()};
    status=read_shm<float>(shm_key, dims_offsets,SVIPC_FLOAT, id_offsets, offsets);

    char id_pupil[] = "pupil";
    int dims_pupil[]={1,sensors.d_wfs[0]->d_pupil->getNbElem()};
    status=read_shm<float>(shm_key, dims_pupil,SVIPC_FLOAT, id_pupil, pupil);

    char id_fluxPerSub[] = "fluxPerSub";
    int dims_fluxPerSub[]={1,sensors.d_wfs[0]->d_fluxPerSub->getNbElem()};
    status=read_shm<float>(shm_key, dims_fluxPerSub,SVIPC_FLOAT, id_fluxPerSub, fluxPerSub);

    char id_isvalid[] = "isvalid";
    int dims_isvalid[]={1,sensors.d_wfs[0]->d_isvalid->getNbElem()};
    status=read_shm<int>(shm_key, dims_isvalid,SVIPC_INT, id_isvalid, isvalid);

    char id_validsubsx[] = "validsubsx";
    int dims_validsubsx[]={1,sensors.d_wfs[0]->d_validsubsx->getNbElem()};
    status=read_shm<int>(shm_key, dims_validsubsx,SVIPC_INT, id_validsubsx, validsubsx);

    char id_validsubsy[] = "validsubsy";
    int dims_validsubsy[]={1,sensors.d_wfs[0]->d_validsubsy->getNbElem()};
    status=read_shm<int>(shm_key, dims_validsubsy,SVIPC_INT, id_validsubsy, validsubsy);

    char id_istart[] = "istart";
    int dims_istart[]={1,sensors.d_wfs[0]->d_istart->getNbElem()};
    status=read_shm<int>(shm_key, dims_istart,SVIPC_INT, id_istart, istart);

    char id_jstart[] = "jstart";
    int dims_jstart[]={1,sensors.d_wfs[0]->d_jstart->getNbElem()};
    status=read_shm<int>(shm_key, dims_jstart,SVIPC_INT, id_jstart, jstart);

    char id_ftkernel[] = "ftkernel";
    int dims_ftkernel[]={1,sensors.d_wfs[0]->d_ftkernel->getNbElem()};
    status=read_shm<float>(shm_key, dims_ftkernel,SVIPC_FLOAT, id_ftkernel, ftkernel);
  }

  MPI_Bcast(phasemap, sensors.d_wfs[0]->d_phasemap->getNbElem(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(hrmap, hrmap_size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(binmap, sensors.d_wfs[0]->d_binmap->getNbElem(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(offsets, sensors.d_wfs[0]->d_offsets->getNbElem(), MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(pupil, sensors.d_wfs[0]->d_pupil->getNbElem(), MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(fluxPerSub, sensors.d_wfs[0]->d_fluxPerSub->getNbElem(), MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(isvalid, sensors.d_wfs[0]->d_isvalid->getNbElem(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(validsubsx, sensors.d_wfs[0]->d_validsubsx->getNbElem(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(validsubsy, sensors.d_wfs[0]->d_validsubsy->getNbElem(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(istart, sensors.d_wfs[0]->d_istart->getNbElem(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(jstart, sensors.d_wfs[0]->d_jstart->getNbElem(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(ftkernel, 2*sensors.d_wfs[0]->d_ftkernel->getNbElem(), MPI_FLOAT, 0, MPI_COMM_WORLD);

  sensors.d_wfs[0]->wfs_initarrays(phasemap, hrmap, binmap, offsets, pupil, fluxPerSub, isvalid,
      validsubsx, validsubsy, istart, jstart, (cuFloatComplex *)ftkernel);

  cudaMemset(sensors.d_wfs[0]->d_gs->d_phase->d_screen->getData(), 0,
      sensors.d_wfs[0]->d_gs->d_phase->d_screen->getNbElem());
  sensors.d_wfs[0]->comp_image();
  const long *dims=sensors.d_wfs[0]->d_bincube->getDims();
  long nPixSppX=dims[1], nPixSppY=dims[2], nSppValid=dims[3];

  if (my_rank == 0) {
    status = svipc_shm_cleanup(shm_key);
    status = svipc_sem_cleanup(sem_key);
  }

  float h_bincube[sensors.d_wfs[0]->d_bincube->getNbElem()];
  sensors.d_wfs[0]->d_bincube->device2host(h_bincube);

  fitsfile *fptr; /* FITS file pointer, defined in fitsio.h */
  fits_create_file(&fptr, "image.fits", &status);
  fits_report_error(stderr, status);

  long naxis = 3;
  long fpixel[] = { 1, 1 , 1 };
  long naxes_ipct[] = { nPixSppX, nPixSppY, nSppValid };
  long nelements = nPixSppX * nPixSppY * nSppValid;
  fits_create_img(fptr, FLOAT_IMG, naxis, naxes_ipct, &status);
  fits_report_error(stderr, status);
  unsigned char pixtype = TFLOAT;
  fits_write_pix(fptr, pixtype, fpixel, nelements, h_bincube, &status);
  fits_report_error(stderr, status);

  fits_close_file(fptr, &status);
  fits_report_error(stderr, status);

  MPI_Finalize();

  return 0;
}
