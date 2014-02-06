#include <sutra_rtc.h>


sutra_rtc::sutra_rtc(carma_context *context)
{
  this->current_context=context;
  this->device = context->get_activeDevice();

}

sutra_rtc::~sutra_rtc()
{
  for (size_t idx = 0; idx < (this->d_centro).size(); idx++) {
   delete this->d_centro[(this->d_centro).size()-1];
   d_centro.pop_back();;
  } 
  for (size_t idx = 0; idx < (this->d_control).size(); idx++) {
   delete this->d_control[(this->d_control).size()-1];
   d_control.pop_back();
  } 

  //delete this->current_context;
}

int sutra_rtc::add_centroider(long nwfs, long nvalid, float offset, float scale, long device, char *typec)
{
  d_centro.push_back(new sutra_centroider(current_context, nwfs, nvalid,offset,scale,device,typec));

  return EXIT_SUCCESS;
}

int sutra_rtc::add_controler(long nactu, long delay, long device,const char *typec)
{
  int ncentroids = 0;

  for (size_t idx = 0; idx < (this->d_centro).size(); idx++) 
   ncentroids += this->d_centro[idx]->nvalid;

  d_control.push_back(new sutra_controler(current_context, ncentroids,nactu,delay,device,typec));

  return EXIT_SUCCESS;
}

int sutra_rtc::rm_controler()
{
  delete this->d_control[(this->d_control).size()-1];
  d_control.pop_back();

  return EXIT_SUCCESS;
}

int sutra_rtc::do_imat(int ncntrl, sutra_sensors *sensors, sutra_dms *ydm)
{
  map<type_screen,sutra_dm *>::iterator p;
  p = ydm->d_dms.begin();
  int inds1;
  inds1 = 0;
  cout << "starting imat measurement" << endl;
  while (p != ydm->d_dms.end()) {
    for (int j = 0; j < p->second->ninflu; j++) {
      p->second->comp_oneactu(j,p->second->push4imat); 
      for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
	int nwfs = this->d_centro[idx_cntr]->nwfs;
	float tmp_noise = sensors->d_wfs[nwfs]->noise;
	sensors->d_wfs[nwfs]->noise = -1;
	sensors->d_wfs[nwfs]->kernconv = true;
	sensors->d_wfs[nwfs]->sensor_trace(ydm, 1);
	sensors->d_wfs[nwfs]->comp_image();
	sensors->d_wfs[nwfs]->noise = tmp_noise;
	sensors->d_wfs[nwfs]->kernconv = false;
      }

      //cout << "actu # " << j << endl;
      do_centroids(ncntrl,sensors,true);
 
      convert_centro(this->d_control[ncntrl]->d_centroids->getData(),
		     this->d_control[ncntrl]->d_centroids->getData(),
		     0,1.0f/p->second->push4imat,
		     this->d_control[ncntrl]->d_centroids->getNbElem(),
		     this->d_control[ncntrl]->device);
      
      cutilSafeCall( cudaMemcpy(&((this->d_control[ncntrl]->d_imat->getData())[inds1]),
				this->d_control[ncntrl]->d_centroids->getData(),
				sizeof(float)* 2 * this->d_control[ncntrl]->nvalid, 
				cudaMemcpyDeviceToDevice));

      p->second->reset_shape();
      inds1 +=  2 * this->d_control[ncntrl]->nvalid;
    }
    p++;
  }

  return EXIT_SUCCESS;
}

int sutra_rtc::do_imat_geom(int ncntrl, sutra_sensors *sensors, sutra_dms *ydm, int type)
{
  map<type_screen,sutra_dm *>::iterator p;
  p = ydm->d_dms.begin();
  int inds1,inds2;
  inds1 = 0;
  while (p != ydm->d_dms.end()) {
    for (int j = 0; j < p->second->ninflu; j++) {
      p->second->comp_oneactu(j, p->second->push4imat); //
      inds2 = 0;
      for (size_t idx_cntr = 0; idx_cntr < (this->d_control).size(); idx_cntr++) {
	int nwfs = this->d_centro[idx_cntr]->nwfs;

	sensors->d_wfs[nwfs]->sensor_trace(ydm, 1);
	//sensors->d_wfs[nwfs]->comp_image();

	float *data=this->d_control[ncntrl]->d_imat->getData();
	sensors->d_wfs[nwfs]->slopes_geom(type,&(data[inds1+inds2]));
	inds2 += 2*sensors->d_wfs[nwfs]->nvalid;
      }
    p->second->reset_shape();
    inds1 +=  2 * this->d_control[ncntrl]->nvalid;
    }
    p++;
  }

  return EXIT_SUCCESS;
}

int sutra_rtc::do_centroids(sutra_sensors *sensors)
{
  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    int nwfs = this->d_centro[idx_cntr]->nwfs;

    if (this->d_centro[idx_cntr]->typec.compare("cog")==0)
      this->d_centro[idx_cntr]->get_cog(sensors->d_wfs[nwfs]);

    if (this->d_centro[idx_cntr]->typec.compare("tcog")==0)
      this->d_centro[idx_cntr]->get_tcog(this->d_centro[idx_cntr]->threshold,sensors->d_wfs[nwfs]);

    if (this->d_centro[idx_cntr]->typec.compare("bpcog")==0)
      this->d_centro[idx_cntr]->get_bpcog(this->d_centro[idx_cntr]->nmax,sensors->d_wfs[nwfs]);

    if (this->d_centro[idx_cntr]->typec.compare("wcog")==0)
      this->d_centro[idx_cntr]->get_wcog(sensors->d_wfs[nwfs]);

    if (this->d_centro[idx_cntr]->typec.compare("corr")==0)
      this->d_centro[idx_cntr]->get_corr(sensors->d_wfs[nwfs]);

    if (this->d_centro[idx_cntr]->typec.compare("pyr")==0)
      this->d_centro[idx_cntr]->get_pyr(sensors->d_wfs[nwfs]);
  }

  return EXIT_SUCCESS;
}

int sutra_rtc::do_centroids(int ncntrl,sutra_sensors *sensors)
{
  return do_centroids(ncntrl,sensors,false);
}

int sutra_rtc::do_centroids(int ncntrl,sutra_sensors *sensors, bool imat)
{
  int inds2 = 0;

  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    int nwfs = this->d_centro[idx_cntr]->nwfs;

    float *odata=this->d_control[ncntrl]->d_centroids->getData();

    if (sensors->d_wfs[nwfs]->nstreams > 1) {
    //if (0) {
      /*
      cutilSafeCall(cudaMemcpy(this->d_centro[idx_cntr]->d_bincube->getData(), sensors->d_wfs[nwfs]->d_bincube->getData(),
			       sizeof(float)*sensors->d_wfs[nwfs]->nvalid*sensors->d_wfs[nwfs]->npix * sensors->d_wfs[nwfs]->npix,
			       cudaMemcpyDeviceToDevice));
      */
      int nstreams = sensors->d_wfs[nwfs]->nstreams;
      for (int i=0;i<nstreams;i++) {
	int istart =  i * sensors->d_wfs[nwfs]->npix * sensors->d_wfs[nwfs]->npix * sensors->d_wfs[nwfs]->nvalid/nstreams;

	
	cutilSafeCall(cudaMemcpyAsync(&((this->d_centro[idx_cntr]->d_bincube->getData())[istart]),
				      &((sensors->d_wfs[nwfs]->d_bincube->getData())[istart]),
				      sizeof(float)*sensors->d_wfs[nwfs]->npix * sensors->d_wfs[nwfs]->npix * sensors->d_wfs[nwfs]->nvalid/nstreams,
				      cudaMemcpyDeviceToDevice, sensors->d_wfs[nwfs]->streams->get_stream(i)));
	
	if (this->d_centro[idx_cntr]->typec.compare("cog")==0)
	  this->d_centro[idx_cntr]->get_cog_async(sensors->d_wfs[nwfs]->streams, 
						  this->d_centro[idx_cntr]->d_bincube->getData(),
						  this->d_centro[idx_cntr]->d_subsum->getData(), 
						  this->d_centro[idx_cntr]->d_centroids->getData(), 
						  sensors->d_wfs[nwfs]->nvalid,sensors->d_wfs[nwfs]->npix);
	
	/*
	  if (this->d_centro[idx_cntr]->typec.compare("tcog")==0)
	  this->d_centro[idx_cntr]->get_tcog(this->d_centro[idx_cntr]->threshold,sensors->d_wfs[nwfs]);
	  
	  if (this->d_centro[idx_cntr]->typec.compare("bpcog")==0)
	  this->d_centro[idx_cntr]->get_bpcog(this->d_centro[idx_cntr]->nmax,sensors->d_wfs[nwfs]);
	  
	  if (this->d_centro[idx_cntr]->typec.compare("wcog")==0)
	  this->d_centro[idx_cntr]->get_wcog(sensors->d_wfs[nwfs]);
	  
	  if (this->d_centro[idx_cntr]->typec.compare("corr")==0)
	  this->d_centro[idx_cntr]->get_corr(sensors->d_wfs[nwfs]);
	*/
	
	istart =  i * sensors->d_wfs[nwfs]->nvalid/nstreams;

	cutilSafeCall(cudaMemcpyAsync(&(odata[inds2+istart]),
				      &((this->d_centro[idx_cntr]->d_centroids->getData())[istart]),
				      sizeof(float)* sensors->d_wfs[nwfs]->nvalid/nstreams,
				      cudaMemcpyDeviceToDevice,sensors->d_wfs[nwfs]->streams->get_stream(i)));

	cutilSafeCall(cudaMemcpyAsync(&(odata[inds2+istart+sensors->d_wfs[nwfs]->nvalid]),
				      &((this->d_centro[idx_cntr]->d_centroids->getData())[istart+sensors->d_wfs[nwfs]->nvalid]),
				      sizeof(float)* sensors->d_wfs[nwfs]->nvalid/nstreams,
				      cudaMemcpyDeviceToDevice,sensors->d_wfs[nwfs]->streams->get_stream(i)));

	
      }
      sensors->d_wfs[nwfs]->streams->wait_all_streams();
      /*
      cutilSafeCall(cudaMemcpy(&(odata[inds2]), this->d_centro[idx_cntr]->d_centroids->getData(),
			       sizeof(float)*2*sensors->d_wfs[nwfs]->nvalid,
			       cudaMemcpyDeviceToDevice));
      */
      inds2 += 2*sensors->d_wfs[nwfs]->nvalid;
    } else {
      if (this->d_centro[idx_cntr]->typec.compare("cog")==0)
	this->d_centro[idx_cntr]->get_cog(sensors->d_wfs[nwfs]);
      
      if (this->d_centro[idx_cntr]->typec.compare("tcog")==0) {
	if (imat) this->d_centro[idx_cntr]->get_cog(sensors->d_wfs[nwfs]);
	else this->d_centro[idx_cntr]->get_tcog(this->d_centro[idx_cntr]->threshold,sensors->d_wfs[nwfs]);
      }
      
      if (this->d_centro[idx_cntr]->typec.compare("bpcog")==0) {
	if (imat) this->d_centro[idx_cntr]->get_cog(sensors->d_wfs[nwfs]);
	else this->d_centro[idx_cntr]->get_bpcog(this->d_centro[idx_cntr]->nmax,sensors->d_wfs[nwfs]);
      }
      
      if (this->d_centro[idx_cntr]->typec.compare("wcog")==0)
	this->d_centro[idx_cntr]->get_wcog(sensors->d_wfs[nwfs]);
      
      if (this->d_centro[idx_cntr]->typec.compare("corr")==0)
	this->d_centro[idx_cntr]->get_corr(sensors->d_wfs[nwfs]);
      
      //cout << this->d_centro[idx_cntr]->typec << endl;

      if (this->d_centro[idx_cntr]->typec.compare("pyr")==0) {
	this->d_centro[idx_cntr]->get_pyr(sensors->d_wfs[nwfs]);
      }

      cutilSafeCall(cudaMemcpy(&(odata[inds2]), sensors->d_wfs[nwfs]->d_slopes->getData(),
			       sizeof(float)*2*sensors->d_wfs[nwfs]->nvalid,
			       cudaMemcpyDeviceToDevice));
      inds2 += 2*sensors->d_wfs[nwfs]->nvalid;
    }
  }

  return EXIT_SUCCESS;
}

int sutra_rtc::do_control(int ncntrl, sutra_dms *ydm)
{

  this->d_control[ncntrl]->frame_delay();

  this->d_control[ncntrl]->comp_com();

  map<type_screen,sutra_dm *>::iterator p;
  p = ydm->d_dms.begin();
  int idx = 0;
  while (p != ydm->d_dms.end()) {
    int nstreams = this->d_control[ncntrl]->nstreams;
    if (nstreams > p->second->ninflu) {
      for (int i=0;i<nstreams;i++) {
	int istart =  i * p->second->ninflu / nstreams;
	cutilSafeCall(cudaMemcpyAsync(&((p->second->d_comm->getData())[istart]),
				      &((this->d_control[ncntrl]->d_com->getData())[idx + istart]),
				      sizeof(float)*p->second->ninflu / nstreams,
				      cudaMemcpyDeviceToDevice,this->d_control[ncntrl]->streams->get_stream(i)));
	p->second->comp_shape();
      }
    } else {
      float *data;
      data=this->d_control[ncntrl]->d_com->getData();
      p->second->comp_shape(&(data[idx]));
    }
    idx += p->second->ninflu;
    p++;
  }


  return EXIT_SUCCESS;
}

