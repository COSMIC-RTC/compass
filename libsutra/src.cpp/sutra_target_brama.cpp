#ifdef USE_BRAMA

#include <sutra_target_brama.h>
#include <sutra_telescope.h>

sutra_target_brama::sutra_target_brama(carma_context *context, ACE_TCHAR *name,
                                       sutra_telescope *d_tel, int subsample_,
                                       int ntargets, float *xpos, float *ypos,
                                       float *lambda, float *mag, float zerop,
                                       long *sizes, int Npts, int device)
  : sutra_target(context, d_tel, ntargets, xpos, ypos, lambda, mag, zerop,
                 sizes, Npts, device) {
  DEBUG_TRACE("init %s", name);
  BRAMA::BRAMA_context brama = BRAMA::BRAMA_context::get_instance(name);
  frame_handle = 0;
  framecounter = 0;
  samplecounter = 0;
  subsample = subsample_;
  is_initialised = 0;

  buff_pixels = NULL;
  dims_pixels = NULL;

  string topics[] = BRAMA_TOPICS;

  if (!brama.is_initialised()) {
    cerr << "brama initialisation failed!" << endl;
    //    throw "brama initialisation failed!";
    return;
  }

  try {
    // Create a subscriber for the command topic
    sub = brama.create_subscriber();
    // Create a publisher for the frame topic
    pub = brama.create_publisher();

    // Create an BRAMA Command listener
    brama.register_command_type(topics[BRAMA::CommandType]);
    cmd_listener = (new sutra_target_bramaListenerImpl);
    cmd_listener_servant =
      dynamic_cast<sutra_target_bramaListenerImpl *>(cmd_listener.in());

    if (CORBA::is_nil(cmd_listener.in())) {
      throw "BRAMA Command listener is nil.";
    }
    cmd_listener_servant->attach_target(this);

    cmd_dr =
      brama.create_datareader(sub, topics[BRAMA::CommandType], cmd_listener);

    // Create an BRAMA Frame writer
    brama.register_frame_type(topics[BRAMA::FrameType]);
    frame_base_dw = brama.create_datawriter(pub, topics[BRAMA::FrameType]);
    if (CORBA::is_nil(frame_base_dw.in())) {
      cerr << "create_datawriter for " << topics[BRAMA::FrameType] << " failed."
           << endl;
      return;
    }
    frame_dw = BRAMA::FrameDataWriter::_narrow(frame_base_dw.in());
    if (CORBA::is_nil(frame_dw.in())) {
      throw "FrameDataWriter could not be narrowed";
    }

    BRAMA::Frame zFrame;
    frame_handle = frame_dw->register_instance(zFrame);

    is_initialised = 1;

  } catch (CORBA::Exception &e) {
    cerr << "Exception caught in main.cpp:" << endl << e << endl;
    ACE_OS::exit(1);
  }
}

sutra_target_brama::~sutra_target_brama() {
  if (!is_initialised) {
    return;
  }

  if (buff_pixels)
    BRAMA::Values::freebuf(buff_pixels);
  if (dims_pixels)
    BRAMA::Dims::freebuf(dims_pixels);
}

void sutra_target_brama::allocateBuffers() {
  if (!is_initialised) {
    return;
  }

  try {
    // TODO : handle targets with different supports...
    dims_pixels = BRAMA::Dims::allocbuf(3);
    dims_pixels[0] = d_targets.size();
    dims_pixels[1] = d_targets[0]->d_leimage->getDims(1);
    dims_pixels[2] = d_targets[0]->d_leimage->getDims(2);

    buff_pixels = BRAMA::Values::allocbuf(d_targets.size() *
                                          d_targets[0]->d_leimage->getNbElem() *
                                          sizeof(float));
  } catch (CORBA::Exception &e) {
    cerr << "Exception caught in main.cpp:" << endl << e << endl;
    ACE_OS::exit(1);
  }
}

void sutra_target_brama::set_subsample(int ntarget, int subsample_) {
  if (!is_initialised) {
    return;
  }

  ACE_Guard<ACE_Mutex> guard(this->lock_);
  this->d_targets[ntarget]->reset_strehlmeter();
  this->samplecounter = 1;
  this->subsample = subsample_;
}

void sutra_target_brama::publish() {
  if (!is_initialised) {
    cerr << "brama not initialised!" << endl;
    return;
  }

  ACE_Guard<ACE_Mutex> guard(this->lock_);
  if (samplecounter % subsample != 0 || subsample <= 0) {
    samplecounter++;
    return;
  }

  if (buff_pixels == NULL)
    allocateBuffers();

  CORBA::Float *buff_pixels_servant = (CORBA::Float *)buff_pixels;

  long idx = 0;
  carma_obj<float> tmp_img(d_targets[0]->current_context,
                           d_targets[0]->d_leimage->getDims());
  for (size_t target = 0; target < d_targets.size(); target++) {
    float flux =
      d_targets[target]->zp * powf(10, -0.4 * d_targets[target]->mag);
    roll_mult<float>(tmp_img.getData(), d_targets[target]->d_leimage->getData(),
                     d_targets[target]->d_leimage->getDims(1),
                     d_targets[target]->d_leimage->getDims(2), flux,
                     d_targets[target]->current_context->get_device(
                       d_targets[target]->device));
    tmp_img.device2host(buff_pixels_servant + idx);

    idx += d_targets[target]->d_leimage->getNbElem();
  }

  BRAMA::Frame zFrame;
  zFrame.framecounter = framecounter;
  zFrame.timestamp = BRAMA::get_timestamp();
  zFrame.source = CORBA::string_dup("BRAMA TARGET");
  zFrame.dimensions = BRAMA::Dims(3, 3, dims_pixels, 0);
  zFrame.data =
    BRAMA::Values(idx * sizeof(float), idx * sizeof(float), buff_pixels, 0);
  zFrame.datatype = BRAMA::BRAMA_float32_t;
  zFrame.sizeofelements = sizeof(float);

  // cout << "Publishing zFrame: " << zFrame.framecounter << endl;
  DDS::ReturnCode_t ret = frame_dw->write(zFrame, frame_handle);

  if (ret != DDS::RETCODE_OK) {
    ACE_ERROR(
      (LM_ERROR, ACE_TEXT("(%P|%t)ERROR: frame write returned %d.\n"), ret));
    return;
  }

  for (size_t target = 0; target < d_targets.size(); target++) {
    d_targets[target]->reset_strehlmeter();
    samplecounter = 0;
  }

  framecounter++;
  // ACE_Time_Value ace_wait(0, 25);
  // ACE_OS::sleep(ace_wait);
}

#endif /* USE_BRAMA */
