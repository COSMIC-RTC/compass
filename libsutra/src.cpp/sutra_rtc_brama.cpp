#ifdef USE_BRAHMA

#include <sutra_rtc_brahma.h>
#include <sutra_rtc_brahmaListenerImpl.h>

sutra_rtc_brahma::sutra_rtc_brahma(carma_context *context, sutra_sensors *wfs_, sutra_target *target_, ACE_TCHAR *name) : sutra_rtc(context), wfs(wfs_), target(target_) {
  DEBUG_TRACE("init %s", name);
  BRAHMA::BRAHMA_context brahma = BRAHMA::BRAHMA_context::get_instance(name);
  cmd_listener_servant = NULL;
  superframe_handle = 0;
  megaframe_handle = 0;
  framecounter = 0;
  is_initialised = 0;

  buff_wfs = NULL;
  buff_wfs_phase = NULL;
  buff_intensities = NULL;
  buff_slopes = NULL;
  buff_commands = NULL;
  buff_target = NULL;
  buff_target_phase = NULL;

  dims_wfs = NULL;
  dims_wfs_phase = NULL;
  dims_intensities = NULL;
  dims_slopes = NULL;
  dims_commands = NULL;
  dims_target = NULL;
  dims_target_phase = NULL;

  string topics[] = BRAHMA_TOPICS;

  if (!brahma.is_initialised()) {
    cerr << "brahma initialisation failed!" << endl;
    //    throw "brahma initialisation failed!";
    return;
  }

  try {
    // Create a subscriber for the command topic
    sub = brahma.create_subscriber();
    // Create a publisher for the megaframe topic
    pub = brahma.create_publisher();

    // Create an BRAHMA Command listener
    brahma.register_command_type(topics[BRAHMA::CommandType]);
    cmd_listener = (new sutra_rtc_brahmaListenerImpl);
    cmd_listener_servant =
      dynamic_cast<sutra_rtc_brahmaListenerImpl *>(cmd_listener.in());

    if (CORBA::is_nil(cmd_listener.in())) {
      throw "BRAHMA Command listener is nil.";
    }
    cmd_listener_servant->attach_rtc(this);

    cmd_dr = brahma.create_datareader(sub, topics[BRAHMA::CommandType], cmd_listener);

    // Create an BRAHMA SuperFrame writer
    brahma.register_superframe_type(topics[BRAHMA::SuperFrameType]);
    superframe_base_dw = brahma.create_datawriter(pub, topics[BRAHMA::SuperFrameType]);
    if (CORBA::is_nil(superframe_base_dw.in())) {
      cerr << "create_datawriter for " << topics[BRAHMA::SuperFrameType] << " failed." << endl;
      return;
    }
    superframe_dw = BRAHMA::SuperFrameDataWriter::_narrow(
                      superframe_base_dw.in());
    if (CORBA::is_nil(superframe_dw.in())) {
      throw "SuperFrameDataWriter could not be narrowed";
    }

    BRAHMA::SuperFrame xFrame;
    superframe_handle = superframe_dw->register_instance(xFrame);

    if (target != NULL) {
      // Create an BRAHMA MegaFrame writer
      brahma.register_megaframe_type(topics[BRAHMA::MegaFrameType]);
      megaframe_base_dw = brahma.create_datawriter(pub, topics[BRAHMA::MegaFrameType]);
      if (CORBA::is_nil(megaframe_base_dw.in())) {
        cerr << "create_datawriter for " << topics[BRAHMA::MegaFrameType] << " failed." << endl;
        return;
      }
      megaframe_dw = BRAHMA::MegaFrameDataWriter::_narrow(
                       megaframe_base_dw.in());
      if (CORBA::is_nil(megaframe_dw.in())) {
        throw "MegaFrameDataWriter could not be narrowed";
      }

      BRAHMA::MegaFrame zFrame;
      megaframe_handle = megaframe_dw->register_instance(zFrame);
    }

    is_initialised = 1;
  } catch (CORBA::Exception &e) {
    cerr << "Exception caught in main.cpp:" << endl
         << e << endl;
    ACE_OS::exit(1);
  }
}

sutra_rtc_brahma::~sutra_rtc_brahma() {
  if (!is_initialised) {
    return;
  }

  if (buff_wfs)
    BRAHMA::Values::freebuf(buff_wfs);
  if (buff_wfs_phase)
    BRAHMA::Values::freebuf(buff_wfs_phase);
  if (buff_intensities)
    BRAHMA::Values::freebuf(buff_intensities);
  if (buff_slopes)
    BRAHMA::Values::freebuf(buff_slopes);
  if (buff_commands)
    BRAHMA::Values::freebuf(buff_commands);
  if (buff_target)
    BRAHMA::Values::freebuf(buff_target);
  if (buff_target_phase)
    BRAHMA::Values::freebuf(buff_target_phase);

  if (dims_wfs)
    BRAHMA::Dims::freebuf(dims_wfs);
  if (dims_wfs_phase)
    BRAHMA::Dims::freebuf(dims_wfs_phase);
  if (dims_intensities)
    BRAHMA::Dims::freebuf(dims_intensities);
  if (dims_slopes)
    BRAHMA::Dims::freebuf(dims_slopes);
  if (dims_commands)
    BRAHMA::Dims::freebuf(dims_commands);
  if (dims_target)
    BRAHMA::Dims::freebuf(dims_target);
  if (dims_target_phase)
    BRAHMA::Dims::freebuf(dims_target_phase);
}

void sutra_rtc_brahma::allocateBuffers() {
  if (!is_initialised) {
    return;
  }

  try {
    wfs_size = 0;
    wfs_phase_size = 0;
    target_size = 0;
    target_phase_size = 0;
    if (target != 0L) {
      for (unsigned int i = 0; i < wfs->d_wfs.size(); i++) {
        wfs_size += wfs->d_wfs[i]->d_binimg->getNbElem();
        wfs_phase_size += wfs->d_wfs[i]->d_gs->d_phase->d_screen->getNbElem();
      }
      for (unsigned int i = 0; i < target->d_targets.size(); i++) {
        target_size += target->d_targets[i]->d_image->getNbElem();
        target_phase_size += target->d_targets[i]->d_phase->d_screen->getNbElem();
      }
    }

    nslp = 0;
    ncmd = 0;
    for (unsigned int i = 0; i < d_control.size(); i++) {
      nslp += d_control[i]->nslope();
      ncmd += d_control[i]->nactu();
    }
    nvalid = nslp / 2;

    buff_intensities = BRAHMA::Values::allocbuf(nvalid * sizeof(float));
    buff_slopes = BRAHMA::Values::allocbuf(nslp * sizeof(float));
    buff_commands = BRAHMA::Values::allocbuf(ncmd * sizeof(float));
    if (target != NULL) {
      buff_wfs = BRAHMA::Values::allocbuf(wfs_size * sizeof(float));
      buff_wfs_phase = BRAHMA::Values::allocbuf(wfs_phase_size * sizeof(float));
      buff_target = BRAHMA::Values::allocbuf(target_size * sizeof(float));
      buff_target_phase = BRAHMA::Values::allocbuf(target_phase_size * sizeof(float));
    } else {
      buff_wfs = NULL;
      buff_target = NULL;
      buff_target_phase = NULL;
    }
    dims_wfs = BRAHMA::Dims::allocbuf(1);
    dims_wfs[0] = wfs_size;

    dims_wfs_phase = BRAHMA::Dims::allocbuf(1);
    dims_wfs_phase[0] = wfs_phase_size;

    dims_intensities = BRAHMA::Dims::allocbuf(1);
    dims_intensities[0] = nvalid;

    dims_slopes = BRAHMA::Dims::allocbuf(1);
    dims_slopes[0] = nslp;

    dims_commands = BRAHMA::Dims::allocbuf(1);
    dims_commands[0] = ncmd;

    dims_target = BRAHMA::Dims::allocbuf(1);
    dims_target[0] = target_size;

    dims_target_phase = BRAHMA::Dims::allocbuf(1);
    dims_target_phase[0] = target_phase_size;

  } catch (CORBA::Exception &e) {
    cerr << "Exception caught in main.cpp:" << endl
         << e << endl;
    ACE_OS::exit(1);
  }
}

void sutra_rtc_brahma::publish() {
  if (!is_initialised) {
    cerr << "brahma not initialised!" << endl;
    return;
  }

  if (buff_intensities == NULL)
    allocateBuffers();
  current_context->set_activeDevice(device, 1);

  CORBA::Float *buff_wfs_servant = (CORBA::Float *)buff_wfs;
  CORBA::Float *buff_wfs_phase_servant = (CORBA::Float *)buff_wfs_phase;
  CORBA::Float *buff_intensities_servant = (CORBA::Float *)buff_intensities;
  CORBA::Float *buff_slopes_servant = (CORBA::Float *)buff_slopes;
  CORBA::Float *buff_commands_servant = (CORBA::Float *)buff_commands;
  CORBA::Float *buff_target_servant = (CORBA::Float *)buff_target;
  CORBA::Float *buff_target_phase_servant = (CORBA::Float *)buff_target_phase;

  int nslp_current = 0;
  int ncmd_current = 0;
  int nvalid_current = 0;

  for (unsigned int i = 0; i < d_control.size(); i++) {
    d_control[i]->d_subsum->device2host(buff_intensities_servant + nvalid_current);
    d_control[i]->d_centroids->device2host(buff_slopes_servant + nslp_current);
    d_control[i]->d_voltage->device2host(buff_commands_servant + ncmd_current);
    nvalid_current += d_control[i]->nslope() / 2;
    nslp_current += d_control[i]->nslope();
    ncmd_current += d_control[i]->nactu();
  }

  if (target != NULL) {
    long idx = 0;
    long idx_phase = 0;
    for (size_t wfs_i = 0; wfs_i < wfs->d_wfs.size(); wfs_i++) {
      if (wfs->d_wfs[wfs_i]->type == "sh")
        wfs->d_wfs[wfs_i]->fill_binimage(0);
      wfs->d_wfs[wfs_i]->d_binimg->device2host(buff_wfs_servant + idx);
      idx += wfs->d_wfs[wfs_i]->d_binimg->getNbElem();
      wfs->d_wfs[wfs_i]->d_gs->d_phase->d_screen->device2host(buff_wfs_phase_servant + idx_phase);
      idx_phase += wfs->d_wfs[wfs_i]->d_gs->d_phase->d_screen->getNbElem();
    }

    idx = 0;
    idx_phase = 0;
    carma_obj<float> tmp_img(target->d_targets[0]->current_context, target->d_targets[0]->d_image->getDims());
    for (size_t i = 0; i < target->d_targets.size(); i++) {
      target->d_targets[i]->comp_image(0, true);
      float flux = 1.0f;
      // target->d_targets[i]->zp * powf(10, -0.4 * target->d_targets[i]->mag);
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

      target->d_targets[i]->d_phase->d_screen->device2host(buff_target_phase_servant + idx_phase);

      idx_phase += target->d_targets[i]->d_phase->d_screen->getNbElem();
    }
  }

  BRAHMA::MegaFrame zFrame;
  zFrame.source = CORBA::string_dup("COMPASS RTC");
  zFrame.framecounter = framecounter;
  zFrame.timestamp = BRAHMA::get_timestamp();

  zFrame.loopData.source = CORBA::string_dup("COMPASS WFSs");

  zFrame.loopData.slps.source = CORBA::string_dup("COMPASS slopes");
  zFrame.loopData.slps.typeofelements = CORBA::string_dup("slopes");
  zFrame.loopData.slps.datatype = BRAHMA::BRAHMA_float32_t;
  zFrame.loopData.slps.sizeofelements = sizeof(float);
  zFrame.loopData.slps.dimensions = BRAHMA::Dims(1, 1, dims_slopes, 0);
  zFrame.loopData.slps.framecounter = framecounter;
  zFrame.loopData.slps.data = BRAHMA::Values(nslp * sizeof(float), nslp * sizeof(float),
                              buff_slopes, 0);
  zFrame.loopData.slps.timestamp = BRAHMA::get_timestamp();

  zFrame.loopData.ints.source = CORBA::string_dup("COMPASS intensities");
  zFrame.loopData.ints.typeofelements = CORBA::string_dup("intensities");
  zFrame.loopData.ints.datatype = BRAHMA::BRAHMA_float32_t;
  zFrame.loopData.ints.sizeofelements = sizeof(float);
  zFrame.loopData.ints.dimensions = BRAHMA::Dims(1, 1, dims_intensities, 0);
  zFrame.loopData.ints.framecounter = framecounter;
  zFrame.loopData.ints.data = BRAHMA::Values(nvalid * sizeof(float),
                              nvalid * sizeof(float), buff_intensities, 0);
  zFrame.loopData.ints.timestamp = BRAHMA::get_timestamp();

  zFrame.loopData.cmds.source = CORBA::string_dup("COMPASS commands");
  zFrame.loopData.cmds.typeofelements = CORBA::string_dup("commands");
  zFrame.loopData.cmds.datatype = BRAHMA::BRAHMA_float32_t;
  zFrame.loopData.cmds.sizeofelements = sizeof(float);
  zFrame.loopData.cmds.dimensions = BRAHMA::Dims(1, 1, dims_commands, 0);
  zFrame.loopData.cmds.framecounter = framecounter;
  zFrame.loopData.cmds.data = BRAHMA::Values(ncmd * sizeof(float), ncmd * sizeof(float),
                              buff_commands, 0);
  zFrame.loopData.cmds.timestamp = BRAHMA::get_timestamp();

  zFrame.loopData.framecounter = framecounter;
  zFrame.loopData.timestamp = BRAHMA::get_timestamp();

  //cout << "Publishing zFrame: " << zFrame.framecounter << endl;
  if (target != NULL) {
    // DEBUG_TRACE("wfs_size %d, target_size %d, target_phase_size %d\n", wfs_size, target_size , target_phase_size);
    zFrame.wfs.source = CORBA::string_dup("COMPASS WFSs");
    zFrame.wfs.typeofelements = CORBA::string_dup("wfs image");
    zFrame.wfs.datatype = BRAHMA::BRAHMA_float32_t;
    zFrame.wfs.sizeofelements = sizeof(float);
    zFrame.wfs.dimensions = BRAHMA::Dims(1, 1, dims_wfs, 0);
    zFrame.wfs.framecounter = framecounter;
    zFrame.wfs.data = BRAHMA::Values(wfs_size * sizeof(float), wfs_size * sizeof(float),
                                     buff_wfs, 0);
    zFrame.wfs.timestamp = BRAHMA::get_timestamp();

    zFrame.wfs_phase.source = CORBA::string_dup("COMPASS WFSs");
    zFrame.wfs_phase.typeofelements = CORBA::string_dup("wfs phase");
    zFrame.wfs_phase.datatype = BRAHMA::BRAHMA_float32_t;
    zFrame.wfs_phase.sizeofelements = sizeof(float);
    zFrame.wfs_phase.dimensions = BRAHMA::Dims(1, 1, dims_wfs_phase, 0);
    zFrame.wfs_phase.framecounter = framecounter;
    zFrame.wfs_phase.data = BRAHMA::Values(wfs_phase_size * sizeof(float), wfs_phase_size * sizeof(float),
                                           buff_wfs_phase, 0);
    zFrame.wfs_phase.timestamp = BRAHMA::get_timestamp();

    zFrame.target.source = CORBA::string_dup("COMPASS Targets");
    zFrame.target.typeofelements = CORBA::string_dup("target image");
    zFrame.target.datatype = BRAHMA::BRAHMA_float32_t;
    zFrame.target.sizeofelements = sizeof(float);
    zFrame.target.dimensions = BRAHMA::Dims(1, 1, dims_target, 0);
    zFrame.target.framecounter = framecounter;
    zFrame.target.data = BRAHMA::Values(target_size * sizeof(float), target_size * sizeof(float),
                                        buff_target, 0);
    zFrame.target.timestamp = BRAHMA::get_timestamp();

    zFrame.target_phase.source = CORBA::string_dup("COMPASS Targets Phase");
    zFrame.target_phase.typeofelements = CORBA::string_dup("target phase");
    zFrame.target_phase.datatype = BRAHMA::BRAHMA_float32_t;
    zFrame.target_phase.sizeofelements = sizeof(float);
    zFrame.target_phase.dimensions = BRAHMA::Dims(1, 1, dims_target_phase, 0);
    zFrame.target_phase.framecounter = framecounter;
    zFrame.target_phase.data = BRAHMA::Values(target_phase_size * sizeof(float), target_phase_size * sizeof(float),
                               buff_target_phase, 0);
    zFrame.target_phase.timestamp = BRAHMA::get_timestamp();

    DDS::ReturnCode_t ret = megaframe_dw->write(zFrame, megaframe_handle);
    if (ret != DDS::RETCODE_OK) {
      ACE_ERROR(
        (LM_ERROR, ACE_TEXT("(%P|%t)ERROR: megaframe_dw write returned %d.\n"), ret));
      return;
    }

    DDS::Duration_t dds_wait = {10, 0};
    ret = megaframe_dw->wait_for_acknowledgments(dds_wait);
    if (ret != DDS::RETCODE_OK) {
      ACE_ERROR(
        (LM_ERROR, ACE_TEXT("(%P|%t)ERROR: megaframe_dw wait_for_acknowledgments returned %d.\n"), ret));
    }
  }

  DDS::ReturnCode_t ret = superframe_dw->write(zFrame.loopData, superframe_handle);

  if (ret != DDS::RETCODE_OK) {
    ACE_ERROR(
      (LM_ERROR, ACE_TEXT("(%P|%t)ERROR: superframe_dw write returned %d.\n"), ret));
    return;
  }

  DDS::Duration_t dds_wait = {10, 0};
  ret = superframe_dw->wait_for_acknowledgments(dds_wait);
  if (ret != DDS::RETCODE_OK) {
    ACE_ERROR(
      (LM_ERROR, ACE_TEXT("(%P|%t)ERROR: superframe_dw wait_for_acknowledgments returned %d.\n"), ret));
  }

  framecounter++;
  //ACE_Time_Value ace_wait(0, 25);
  //ACE_OS::sleep(ace_wait);
}

#endif /* USE_BRAHMA */
