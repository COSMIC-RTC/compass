#ifdef USE_BRAMA

#include<sutra_rtc_brama.h>
#include<sutra_rtc_bramaListenerImpl.h>

sutra_rtc_brama::sutra_rtc_brama(carma_context *context, sutra_sensors *wfs_, sutra_target *target_, ACE_TCHAR* name) :
sutra_rtc(context), wfs(wfs_), target(target_) {
  brama = new BRAMA_supervisor(name);
  cmd_listener_servant = NULL;
  superframe_handle = 0;
  megaframe_handle = 0;
  framecounter = 0;

  buff_wfs = NULL;
  buff_intensities = NULL;
  buff_slopes = NULL;
  buff_commands = NULL;
  buff_target = NULL;

  dims_wfs = NULL;
  dims_intensities = NULL;
  dims_slopes = NULL;
  dims_commands = NULL;
  dims_target = NULL;

  string topics[] = BRAMA_TOPICS;

  try {
    // Create a subscriber for the command topic
    brama->create_subscriber();
    // Create a publisher for the megaframe topic
    brama->create_publisher();

    // Register the BRAMA types
    brama->register_all_data_types();

    // Create an BRAMA Command listener
    cmd_listener = (new sutra_rtc_bramaListenerImpl);
    cmd_listener_servant =
    dynamic_cast<sutra_rtc_bramaListenerImpl*>(cmd_listener.in());

    if (CORBA::is_nil(cmd_listener.in())) {
      throw "BRAMA Command listener is nil.";
    }
    cmd_listener_servant->attach_rtc(this);

    cmd_dr = brama->create_datareader(topics[CommandType], cmd_listener);

    superframe_base_dw = brama->create_datawriter(topics[SuperFrameType]);
    if (CORBA::is_nil(superframe_base_dw.in())) {
      cerr << "create_datawriter for " << topics[SuperFrameType] << " failed." << endl;
      ACE_OS::exit(1);
    }
    superframe_dw = BRAMA::SuperFrameDataWriter::_narrow(
        superframe_base_dw.in());
    if (CORBA::is_nil(superframe_dw.in())) {
      throw "SuperFrameDataWriter could not be narrowed";
    }

    BRAMA::SuperFrame xFrame;
    superframe_handle = superframe_dw->register_instance(xFrame);

    megaframe_base_dw = brama->create_datawriter(topics[MegaFrameType]);
    if (CORBA::is_nil(megaframe_base_dw.in())) {
      cerr << "create_datawriter for " << topics[MegaFrameType] << " failed." << endl;
      ACE_OS::exit(1);
    }
    megaframe_dw = BRAMA::MegaFrameDataWriter::_narrow(
        megaframe_base_dw.in());
    if (CORBA::is_nil(megaframe_dw.in())) {
      throw "MegaFrameDataWriter could not be narrowed";
    }

    BRAMA::MegaFrame zFrame;
    megaframe_handle = megaframe_dw->register_instance(zFrame);

  } catch (CORBA::Exception& e) {
    cerr << "Exception caught in main.cpp:" << endl << e << endl;
    ACE_OS::exit(1);
  }
}

sutra_rtc_brama::~sutra_rtc_brama() {

  if (buff_wfs)
  BRAMA::Values::freebuf(buff_wfs);
  if (buff_intensities)
  BRAMA::Values::freebuf(buff_intensities);
  if (buff_slopes)
  BRAMA::Values::freebuf(buff_slopes);
  if (buff_commands)
  BRAMA::Values::freebuf(buff_commands);
  if (buff_target)
  BRAMA::Values::freebuf(buff_target);

  if (dims_wfs)
  BRAMA::Dims::freebuf(dims_wfs);
  if (dims_intensities)
  BRAMA::Dims::freebuf(dims_intensities);
  if (dims_slopes)
  BRAMA::Dims::freebuf(dims_slopes);
  if (dims_commands)
  BRAMA::Dims::freebuf(dims_commands);
  if (dims_target)
  BRAMA::Dims::freebuf(dims_target);

  delete brama;
}

void sutra_rtc_brama::allocateBuffers() {
  try {
    wfs_size = 0;
    for (unsigned int i = 0; i < wfs->d_wfs.size(); i++) {
      wfs_size += wfs->d_wfs[i]->d_binimg->getNbElem();
    }

    target_size = 0;
    for (unsigned int i = 0; i < target->d_targets.size(); i++) {
      target_size += target->d_targets[i]->d_image->getNbElem();
    }

    nslp = 0;
    ncmd = 0;
    for (unsigned int i = 0; i < d_control.size(); i++) {
      nslp += d_control[i]->nslope();
      ncmd += d_control[i]->nactu();
    }
    nvalid = nslp/2;

    buff_wfs = BRAMA::Values::allocbuf(wfs_size * sizeof(float));
    buff_intensities = BRAMA::Values::allocbuf(nvalid * sizeof(float));
    buff_slopes = BRAMA::Values::allocbuf(nslp * sizeof(float));
    buff_commands = BRAMA::Values::allocbuf(ncmd * sizeof(float));
    buff_target = BRAMA::Values::allocbuf(target_size * sizeof(float));

    dims_wfs = BRAMA::Dims::allocbuf(2);
    dims_wfs[0] = 1;
    dims_wfs[1] = wfs_size;

    dims_intensities = BRAMA::Dims::allocbuf(2);
    dims_intensities[0] = 1;
    dims_intensities[1] = nvalid;

    dims_slopes = BRAMA::Dims::allocbuf(2);
    dims_slopes[0] = 1;
    dims_slopes[1] = nslp;

    dims_commands = BRAMA::Dims::allocbuf(2);
    dims_commands[0] = 1;
    dims_commands[1] = ncmd;

    dims_target = BRAMA::Dims::allocbuf(2);
    dims_target[0] = 1;
    dims_target[1] = target_size;

  } catch (CORBA::Exception& e) {
    cerr << "Exception caught in main.cpp:" << endl << e << endl;
    ACE_OS::exit(1);
  }
}

void sutra_rtc_brama::publish() {
  if(buff_intensities == NULL)
  allocateBuffers();
  current_context->set_activeDevice(device,1);

  CORBA::Float* buff_wfs_servant = (CORBA::Float*) buff_wfs;
  CORBA::Float* buff_intensities_servant = (CORBA::Float*) buff_intensities;
  CORBA::Float* buff_slopes_servant = (CORBA::Float*) buff_slopes;
  CORBA::Float* buff_commands_servant = (CORBA::Float*) buff_commands;
  CORBA::Float* buff_target_servant = (CORBA::Float*) buff_target;

  long idx = 0;
  for (size_t wfs_i = 0; wfs_i < wfs->d_wfs.size(); wfs_i++) {
    wfs->d_wfs[wfs_i]->fill_binimage(0);
    wfs->d_wfs[wfs_i]->d_binimg->device2host(buff_wfs_servant + idx);
    idx += wfs->d_wfs[wfs_i]->d_binimg->getNbElem();
  }

  int nslp_current = 0;
  int ncmd_current = 0;
  int nvalid_current = 0;

  for (unsigned int i = 0; i < d_control.size(); i++) {
    d_control[i]->d_subsum->device2host(buff_intensities_servant + nvalid_current);
    d_control[i]->d_centroids->device2host(buff_slopes_servant + nslp_current);
    d_control[i]->d_voltage->device2host(buff_commands_servant + ncmd_current);
    nvalid_current += d_control[i]->nslope()/2;
    nslp_current += d_control[i]->nslope();
    ncmd_current += d_control[i]->nactu();
  }

  idx = 0;
  for (size_t i = 0; i < target->d_targets.size(); i++) {
    target->d_targets[i]->comp_image(0);
    carma_obj<float> tmp_img(target->d_targets[i]->current_context, target->d_targets[i]->d_image->getDims());
    float flux = target->d_targets[i]->zp
        * powf(10, -0.4 * target->d_targets[i]->mag);
    roll_mult<float>(
        tmp_img.getData(),
        target->d_targets[i]->d_image->getData(),
        target->d_targets[i]->d_image->getDims(1),
        target->d_targets[i]->d_image->getDims(2),
        flux,
        target->d_targets[i]->current_context->get_device(
            target->d_targets[i]->device));
    tmp_img.device2host(buff_target_servant + idx);

    idx += target->d_targets[i]->d_image->getNbElem();
  }

  BRAMA::MegaFrame zFrame;
  zFrame.framecounter = framecounter;
  zFrame.timestamp = get_timestamp();
  zFrame.from = CORBA::string_dup("COMPASS RTC");

  zFrame.wfs.framecounter = framecounter;
  zFrame.wfs.timestamp = get_timestamp();
  zFrame.wfs.from = CORBA::string_dup("COMPASS WFSs");
  zFrame.wfs.dimensions = BRAMA::Dims(2, 2, dims_wfs, 0);
  zFrame.wfs.data = BRAMA::Values(wfs_size * sizeof(float), wfs_size * sizeof(float),
      buff_wfs, 0);
  zFrame.wfs.sizeofelements = sizeof(float);

  zFrame.loopData.ints.from = CORBA::string_dup("COMPASS intensities");
  zFrame.loopData.ints.typeofelements = CORBA::string_dup("intensities");
  zFrame.loopData.ints.sizeofelements = sizeof(float);
  zFrame.loopData.ints.dimensions = BRAMA::Dims(2, 2, dims_intensities, 0);
  zFrame.loopData.ints.framecounter = framecounter;
  zFrame.loopData.ints.data = BRAMA::Values(nvalid * sizeof(float),
      nvalid * sizeof(float), buff_intensities, 0);
  zFrame.loopData.ints.timestamp = get_timestamp();

  zFrame.loopData.slps.from = CORBA::string_dup("COMPASS slopes");
  zFrame.loopData.slps.typeofelements = CORBA::string_dup("slopes");
  zFrame.loopData.slps.sizeofelements = sizeof(float);
  zFrame.loopData.slps.dimensions = BRAMA::Dims(2, 2, dims_slopes, 0);
  zFrame.loopData.slps.framecounter = framecounter;
  zFrame.loopData.slps.data = BRAMA::Values(nslp * sizeof(float), nslp * sizeof(float),
      buff_slopes, 0);
  zFrame.loopData.slps.timestamp = get_timestamp();

  zFrame.loopData.cmds.from = CORBA::string_dup("COMPASS commands");
  zFrame.loopData.cmds.typeofelements = CORBA::string_dup("commands");
  zFrame.loopData.cmds.sizeofelements = sizeof(float);
  zFrame.loopData.cmds.dimensions = BRAMA::Dims(2, 2, dims_commands, 0);
  zFrame.loopData.cmds.framecounter = framecounter;
  zFrame.loopData.cmds.data = BRAMA::Values(ncmd * sizeof(float), ncmd * sizeof(float),
      buff_commands, 0);
  zFrame.loopData.cmds.timestamp = get_timestamp();

  zFrame.target.framecounter = framecounter;
  zFrame.target.timestamp = get_timestamp();
  zFrame.target.from = CORBA::string_dup("COMPASS Targets");
  zFrame.target.dimensions = BRAMA::Dims(2, 2, dims_target, 0);
  zFrame.target.data = BRAMA::Values(target_size * sizeof(float), target_size * sizeof(float),
      buff_target, 0);
  zFrame.target.sizeofelements = sizeof(float);

//cout << "Publishing zFrame: " << zFrame.framecounter << endl;
  DDS::ReturnCode_t ret = megaframe_dw->write(zFrame, megaframe_handle);

  if (ret != DDS::RETCODE_OK) {
    ACE_ERROR(
        (LM_ERROR, ACE_TEXT("(%P|%t)ERROR: megaframe write returned %d.\n"), ret));
    return;
  }

  ret = superframe_dw->write(zFrame.loopData, superframe_handle);

  if (ret != DDS::RETCODE_OK) {
    ACE_ERROR(
        (LM_ERROR, ACE_TEXT("(%P|%t)ERROR: superframe write returned %d.\n"), ret));
    return;
  }

  framecounter++;
//ACE_Time_Value ace_wait(0, 25);
//ACE_OS::sleep(ace_wait);
}

#endif /* USE_BRAMA */
