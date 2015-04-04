#ifdef USE_BRAMA

#include<sutra_rtc_brama.h>
#include<sutra_rtc_bramaListenerImpl.h>

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

  string topics[] = BRAMA_TOPICS;

  try {
    // Create a subscriber for the command topic
    brama->create_subscriber();
    // Create a publisher for the superframe topic
    brama->create_publisher();

    // Register the BRAMA types
    brama->register_all_data_types();

    // Create an BRAMA Event listener
    cmd_listener = (new sutra_rtc_bramaListenerImpl);
    cmd_listener_servant =
        dynamic_cast<sutra_rtc_bramaListenerImpl*>(cmd_listener.in());

    if (CORBA::is_nil(cmd_listener.in())) {
      throw "BRAMA Command listener is nil.";
    }
    cmd_listener_servant->attach_rtc(this);

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

    BRAMA::SuperFrame zFrame;
    superframe_handle = superframe_dw->register_instance(zFrame);

} catch (CORBA::Exception& e) {
    cerr << "Exception caught in main.cpp:" << endl << e << endl;
    ACE_OS::exit(1);
  }
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

void sutra_rtc_brama::allocateBuffers() {
  try {
    int nslp = 0;
    int ncmd = 0;
    for (unsigned int i = 0; i < d_control.size(); i++) {
      nslp   += d_control[i]->nslope();
      ncmd   += d_control[i]->nactu();
    }
    int nvalid = nslp/2;

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
  if(buff_intensities == NULL)
    allocateBuffers();
  current_context->set_activeDevice(device,1);


  CORBA::Float* buff_intensities_servant = (CORBA::Float*) buff_intensities;
  CORBA::Float* buff_slopes_servant = (CORBA::Float*) buff_slopes;
  CORBA::Float* buff_commands_servant = (CORBA::Float*) buff_commands;

  int nslp = 0;
  int ncmd = 0;
  int nvalid = 0;

  for (unsigned int i = 0; i < d_control.size(); i++) {
    d_control[i]->d_subsum->device2host(buff_intensities_servant + nvalid);
    d_control[i]->d_centroids->device2host(buff_slopes_servant + nslp);
    d_control[i]->d_voltage->device2host(buff_commands_servant + ncmd);
    nvalid += d_control[i]->nslope()/2;
    nslp += d_control[i]->nslope();
    ncmd += d_control[i]->nactu();
  }

  BRAMA::SuperFrame zFrame;
  zFrame.framecounter = framecounter;
  zFrame.timestamp = get_timestamp();
  zFrame.from = CORBA::string_dup("BRAMA RTC");

  zFrame.ints.from = CORBA::string_dup("BRAMA intensities");
  zFrame.ints.typeofelements = CORBA::string_dup("intensities");
  zFrame.ints.sizeofelements = sizeof(float);
  zFrame.ints.dimensions = BRAMA::Dims(2, 2, dims_intensities, 0);
  zFrame.ints.framecounter = framecounter;
  zFrame.ints.data = BRAMA::Values(nvalid * sizeof(float),
                                   nvalid * sizeof(float), buff_intensities, 0);
  zFrame.ints.timestamp = get_timestamp();

  zFrame.slps.from = CORBA::string_dup("BRAMA slopes");
  zFrame.slps.typeofelements = CORBA::string_dup("slopes");
  zFrame.slps.sizeofelements = sizeof(float);
  zFrame.slps.dimensions = BRAMA::Dims(2, 2, dims_slopes, 0);
  zFrame.slps.framecounter = framecounter;
  zFrame.slps.data = BRAMA::Values(nslp * sizeof(float), nslp * sizeof(float),
                                   buff_slopes, 0);
  zFrame.slps.timestamp = get_timestamp();

  zFrame.cmds.from = CORBA::string_dup("BRAMA commands");
  zFrame.cmds.typeofelements = CORBA::string_dup("commands");
  zFrame.cmds.sizeofelements = sizeof(float);
  zFrame.cmds.dimensions = BRAMA::Dims(2, 2, dims_commands, 0);
  zFrame.cmds.framecounter = framecounter;
  zFrame.cmds.data = BRAMA::Values(ncmd * sizeof(float), ncmd * sizeof(float),
                                   buff_commands, 0);
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
