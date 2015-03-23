#ifdef USE_BRAMA

#include<sutra_rtc_brama.h>

sutra_rtc_brama::sutra_rtc_brama(carma_context *context, ACE_TCHAR* name) :
    sutra_rtc(context) {
  brama = new BRAMA_supervisor(name);
  cmd_listener_servant = NULL;
  superframe_handle = 0;
  framecounter = 0;

  buff_intensities = NULL;
  buff_slopes = NULL;
  buff_commands = NULL;

  dims_intensities = NULL;
  dims_slopes = NULL;
  dims_commands = NULL;

}

sutra_rtc_brama::~sutra_rtc_brama() {

  if (buff_intensities)
    BRAMA::Values::freebuf(buff_intensities);
  if (buff_slopes)
    BRAMA::Values::freebuf(buff_slopes);
  if (buff_commands)
    BRAMA::Values::freebuf(buff_commands);

  if (dims_intensities)
    BRAMA::Dims::freebuf(dims_intensities);
  if (dims_slopes)
    BRAMA::Dims::freebuf(dims_slopes);
  if (dims_commands)
    BRAMA::Dims::freebuf(dims_commands);

  delete brama;
}

void sutra_rtc_brama::update_param() {
  try {
    BRAMA::Command cmd = cmd_listener_servant->getlastCmd();
    const char *type = cmd.name.in();
    DEBUG_TRACE("get command type %s", type);
    if (strncmp(type, "CM_", 3) == 0) {
      int ncontrol = atoi(&type[3]);

      DEBUG_TRACE("Updating control matrix num %d", ncontrol);

      unsigned int nslope = d_control[ncontrol]->nslope();
      unsigned int ncmd = d_control[ncontrol]->nactu();

      if (cmd.dimensions[0] != 2 || cmd.dimensions[1] != ncmd
          || cmd.dimensions[2] != nslope) {
        BRAMA_DEBUG_TRACE("wrong dimensions : %d %d %d", cmd.dimensions[0],
                          cmd.dimensions[1], cmd.dimensions[2]);
        BRAMA_DEBUG_TRACE("it should be : 2 %d %d", ncmd, nslope);
        throw CORBA::BAD_PARAM();
      }

      CORBA::Float *data = (CORBA::Float*)cmd.data.get_buffer();
      if (d_control[ncontrol]->get_type() == "ls") {
        sutra_controller_ls *control =
            dynamic_cast<sutra_controller_ls *>(d_control[ncontrol]);
        control->d_cmat->host2device(data);
        return;
      } else if (d_control[ncontrol]->get_type() == "ls") {
        sutra_controller_mv *control =
            dynamic_cast<sutra_controller_mv *>(d_control[ncontrol]);
        control->d_cmat->host2device(data);
        return;
      } else {
        BRAMA_DEBUG_TRACE("controller %d must be a ls or mv controller",
                          ncontrol);
        throw CORBA::BAD_PARAM();
      }
    } else if (strncmp(type, "PertVolt_",9) == 0) {
      int ncontrol = atoi(&type[9]);
      DEBUG_TRACE("Updating perturbation voltages on controller %d", ncontrol);

      unsigned int ncmd = d_control[ncontrol]->nactu();

      if (cmd.dimensions[0] != 2 || cmd.dimensions[1] != ncmd) {
        BRAMA_DEBUG_TRACE("wrong dimensions : %d %d %d", cmd.dimensions[0],
                          cmd.dimensions[1], cmd.dimensions[2]);
        BRAMA_DEBUG_TRACE("it should be : 2 %d nb_elem", ncmd);
        throw CORBA::BAD_PARAM();
      }

      CORBA::Float *data = (CORBA::Float*)cmd.data.get_buffer();
      for(int idx=0; idx<ncmd; idx++){
        BRAMA_DEBUG_TRACE("cmd %d: %f", idx, data[idx]);
      }
      BRAMA_DEBUG_TRACE("nb_elem %d", cmd.dimensions[2]);
      d_control[ncontrol]->set_perturbcom(data, cmd.dimensions[2]);
      return;
    } else {
      throw "Unknown parameter";
    }
  } catch (char const*msg) {
    //BRAMA_DEBUG_TRACE("%s",msg);
  } catch (CORBA::BAD_PARAM &p) {
    std::ostringstream stm;
    stm << p;
    BRAMA_DEBUG_TRACE("%s", stm.str().c_str());
  }
}

void sutra_rtc_brama::initDDS() {
  string topics[] = BRAMA_TOPICS;
  current_context->set_activeDevice(device,1);

  try {
    // Create a publisher and subscriber for the two topics
    brama->create_publisher();
    brama->create_subscriber();

    // Register the BRAMA types
    brama->register_all_data_types();

    // Create an BRAMA Event listener
    cmd_listener = (new CommandDataReaderListenerImpl);
    cmd_listener_servant =
        dynamic_cast<CommandDataReaderListenerImpl*>(cmd_listener.in());

    if (CORBA::is_nil(cmd_listener.in())) {
      throw "BRAMA Command listener is nil.";
    }

    cmd_dr = brama->create_datareader(topics[3], cmd_listener);

    superframe_base_dw = brama->create_datawriter(topics[1]);
    if (CORBA::is_nil(superframe_base_dw.in())) {
      cerr << "create_datawriter for " << topics[1] << " failed." << endl;
      ACE_OS::exit(1);
    }
    superframe_dw = BRAMA::SuperFrameDataWriter::_narrow(
        superframe_base_dw.in());
    if (CORBA::is_nil(superframe_dw.in())) {
      throw "SuperFrameDataWriter could not be narrowed";
    }

    BRAMA::SuperFrame superframe;
    superframe_handle = superframe_dw->register_instance(superframe);

    int nvalid = 0;
    for (unsigned int i = 0; i < d_centro.size(); i++) {
      nvalid += d_centro[i]->nvalid;
    }

    int nslp = 0;
    int ncmd = 0;
    for (unsigned int i = 0; i < d_control.size(); i++) {
      nslp += d_control[i]->nslope();
      ncmd += d_control[i]->nactu();
    }

    //cerr << nvalid << " " << ncmd << endl;

    buff_intensities = BRAMA::Values::allocbuf(nvalid * sizeof(float));
    buff_slopes = BRAMA::Values::allocbuf(nslp * sizeof(float));
    buff_commands = BRAMA::Values::allocbuf(ncmd * sizeof(float));

    dims_intensities = BRAMA::Dims::allocbuf(2);
    dims_intensities[0] = 1;
    dims_intensities[1] = nvalid;

    dims_slopes = BRAMA::Dims::allocbuf(2);
    dims_slopes[0] = 1;
    dims_slopes[1] = nslp;

    dims_commands = BRAMA::Dims::allocbuf(2);
    dims_commands[0] = 1;
    dims_commands[1] = ncmd;

  } catch (CORBA::Exception& e) {
    cerr << "Exception caught in main.cpp:" << endl << e << endl;
    ACE_OS::exit(1);
  }
}

void sutra_rtc_brama::publish() {

  current_context->set_activeDevice(device,1);

  update_param();

  BRAMA::SuperFrame zFrame;
  zFrame.from = CORBA::string_dup("BRAMA RTC");
  zFrame.framecounter = framecounter;

  CORBA::Float* buff_intensities_servant = (CORBA::Float*) buff_intensities;
  CORBA::Float* buff_slopes_servant = (CORBA::Float*) buff_slopes;
  CORBA::Float* buff_commands_servant = (CORBA::Float*) buff_commands;

  int nvalid = 0;
  for (unsigned int i = 0; i < d_centro.size(); i++) {
    d_centro[i]->wfs->d_subsum->device2host(buff_intensities_servant + nvalid);
    nvalid += d_centro[i]->nvalid;
  }

  int nslp = 0;
  int ncmd = 0;
  for (unsigned int i = 0; i < d_control.size(); i++) {
    d_control[i]->d_centroids->device2host(buff_slopes_servant + nslp);
    d_control[i]->d_com->device2host(buff_commands_servant + ncmd);
    nslp += d_control[i]->nslope();
    ncmd += d_control[i]->nactu();
  }

  zFrame.ints.from = CORBA::string_dup("BRAMA intensities");
  zFrame.ints.typeofelements = CORBA::string_dup("intensities");
  zFrame.ints.sizeofelements = sizeof(float);
  zFrame.ints.dimensions = BRAMA::Dims(2, 2, dims_intensities, 0);
  zFrame.ints.framecounter = framecounter;
  zFrame.ints.data = BRAMA::Values(nvalid, nvalid, buff_intensities, 0);
  zFrame.ints.timestamp = get_timestamp();

  zFrame.slps.from = CORBA::string_dup("BRAMA slopes");
  zFrame.slps.typeofelements = CORBA::string_dup("slopes");
  zFrame.slps.sizeofelements = sizeof(float);
  zFrame.slps.dimensions = BRAMA::Dims(2, 2, dims_slopes, 0);
  zFrame.slps.framecounter = framecounter;
  zFrame.slps.data = BRAMA::Values(nslp, nslp, buff_slopes, 0);
  zFrame.slps.timestamp = get_timestamp();

  zFrame.cmds.from = CORBA::string_dup("BRAMA commands");
  zFrame.cmds.typeofelements = CORBA::string_dup("commands");
  zFrame.cmds.sizeofelements = sizeof(float);
  zFrame.cmds.dimensions = BRAMA::Dims(2, 2, dims_commands, 0);
  zFrame.cmds.framecounter = framecounter;
  zFrame.cmds.data = BRAMA::Values(ncmd, ncmd, buff_commands, 0);
  zFrame.cmds.timestamp = get_timestamp();

 //cout << "Publishing zFrame: " << zFrame.framecounter << endl;
  DDS::ReturnCode_t ret = superframe_dw->write(zFrame, superframe_handle);

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
