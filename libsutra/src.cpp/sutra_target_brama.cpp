#ifdef USE_BRAMA

#include<sutra_target_brama.h>

sutra_target_brama::sutra_target_brama(carma_context *context, ACE_TCHAR* name,
                                       int subsample_, int ntargets,
                                       float *xpos, float *ypos, float *lambda,
                                       float *mag, float zerop, long *sizes, float *pupil,
                                       int Npts, int device) :
    sutra_target(context, ntargets, xpos, ypos, lambda, mag, zerop, sizes, pupil, Npts,
                 device) {
  brama = new BRAMA_supervisor(name);
  frame_handle = 0;
  framecounter = 0;
  subsample = subsample_;

  buff_pixels = NULL;
  dims_pixels = NULL;

  string topics[] = BRAMA_TOPICS;

  try {
    // Create a subscriber for the command topic
    brama->create_subscriber();
    // Create a publisher for the frame topic
    brama->create_publisher();

    // Register the BRAMA types
    brama->register_all_data_types();

    // Create an BRAMA Command listener
    cmd_listener = (new sutra_target_bramaListenerImpl);
    cmd_listener_servant =
        dynamic_cast<sutra_target_bramaListenerImpl*>(cmd_listener.in());

    if (CORBA::is_nil(cmd_listener.in())) {
      throw "BRAMA Command listener is nil.";
    }
    cmd_listener_servant->attach_target(this);

    cmd_dr = brama->create_datareader(topics[3], cmd_listener);

    frame_base_dw = brama->create_datawriter(topics[0]);
    if (CORBA::is_nil(frame_base_dw.in())) {
      cerr << "create_datawriter for " << topics[0] << " failed." << endl;
      ACE_OS::exit(1);
    }
    frame_dw = BRAMA::FrameDataWriter::_narrow(frame_base_dw.in());
    if (CORBA::is_nil(frame_dw.in())) {
      throw "FrameDataWriter could not be narrowed";
    }

    BRAMA::Frame zFrame;
    frame_handle = frame_dw->register_instance(zFrame);

  } catch (CORBA::Exception& e) {
    cerr << "Exception caught in main.cpp:" << endl << e << endl;
    ACE_OS::exit(1);
  }
}

sutra_target_brama::~sutra_target_brama() {

  if (buff_pixels)
    BRAMA::Values::freebuf(buff_pixels);
  if (dims_pixels)
    BRAMA::Dims::freebuf(dims_pixels);

  delete brama;
}

void sutra_target_brama::allocateBuffers() {
  try {
    //TODO : handle targets with different supports...
    dims_pixels = BRAMA::Dims::allocbuf(4);
    dims_pixels[0] = 4;
    dims_pixels[1] = d_targets.size();
    dims_pixels[2] = d_targets[0]->d_leimage->getDims(1);
    dims_pixels[3] = d_targets[0]->d_leimage->getDims(2);

    buff_pixels = BRAMA::Values::allocbuf(
        d_targets.size() * d_targets[0]->d_leimage->getNbElem()
        * sizeof(float));
  } catch (CORBA::Exception& e) {
    cerr << "Exception caught in main.cpp:" << endl << e << endl;
    ACE_OS::exit(1);
  }
}

void sutra_target_brama::set_subsample(int ntarget, int subsample_){
  ACE_Guard<ACE_Mutex> guard(this->lock_);
  this->d_targets[ntarget]->reset_strehlmeter();
  this->framecounter=1;
  this->subsample=subsample_;
}

void sutra_target_brama::publish() {
  ACE_Guard<ACE_Mutex> guard(this->lock_);
  if (framecounter%subsample!=0 || subsample<=0) {
    framecounter++;
    return;
  }

  if (buff_pixels == NULL)
    allocateBuffers();

  CORBA::Float* buff_pixels_servant = (CORBA::Float*) buff_pixels;

  long idx = 0;
  carma_obj<float> tmp_img(d_targets[0]->current_context, d_targets[0]->d_leimage->getDims());
  for (size_t target = 0; target < d_targets.size(); target++) {
    float flux = d_targets[target]->zp
        * powf(10, -0.4 * d_targets[target]->mag);
    roll_mult<float>(
        tmp_img.getData(),
        d_targets[target]->d_leimage->getData(),
        d_targets[target]->d_leimage->getDims(1),
        d_targets[target]->d_leimage->getDims(2),
        flux,
        d_targets[target]->current_context->get_device(
            d_targets[target]->device));
    tmp_img.device2host(buff_pixels_servant + idx);

    idx += d_targets[target]->d_leimage->getNbElem();
  }

  BRAMA::Frame zFrame;
  zFrame.framecounter = framecounter;
  zFrame.timestamp = get_timestamp();
  zFrame.from = CORBA::string_dup("BRAMA TARGET");
  zFrame.dimensions = BRAMA::Dims(4, 4, dims_pixels, 0);
  zFrame.data = BRAMA::Values(idx * sizeof(float), idx * sizeof(float),
                              buff_pixels, 0);
  zFrame.sizeofelements = sizeof(float);

//cout << "Publishing zFrame: " << zFrame.framecounter << endl;
  DDS::ReturnCode_t ret = frame_dw->write(zFrame, frame_handle);

  if (ret != DDS::RETCODE_OK) {
    ACE_ERROR(
        (LM_ERROR, ACE_TEXT("(%P|%t)ERROR: frame write returned %d.\n"), ret));
    return;
  }

  for (size_t target = 0; target < d_targets.size(); target++) {
    d_targets[target]->reset_strehlmeter();
    framecounter=0;
  }

  framecounter++;
//ACE_Time_Value ace_wait(0, 25);
//ACE_OS::sleep(ace_wait);
}

#endif /* USE_BRAMA */