#ifdef USE_BRAMA

#include<sutra_rtc_brama.h>

sutra_rtc_brama::sutra_rtc_brama(carma_context *context, ACE_TCHAR* name) :
    sutra_rtc(context) {
  brama = new BRAMA_supervisor(name);
  framecounter = 0;
  initDDS();
}

sutra_rtc_brama::~sutra_rtc_brama() {
}

void sutra_rtc_brama::initDDS() {
  string topics[] = BRAMA_TOPICS;

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
      cerr << "BRAMA Command listener is nil." << endl;
      ACE_OS::exit(1);
    }

    cmd_dr = brama->create_datareader(topics[2], cmd_listener);

    superframe_base_dw = brama->create_datawriter(topics[1]);
    if (CORBA::is_nil(superframe_base_dw.in())) {
      cerr << "create_datawriter for " << topics[1] << " failed." << endl;
      ACE_OS::exit(1);
    }
    superframe_dw = BRAMA::SuperFrameDataWriter::_narrow(
        superframe_base_dw.in());
    if (CORBA::is_nil(superframe_dw.in())) {
      cerr << "SuperFrameDataWriter could not be narrowed" << endl;
      ACE_OS::exit(1);
    }

    superframe_handle = superframe_dw->register_instance(superframe);

  } catch (CORBA::Exception& e) {
    cerr << "Exception caught in main.cpp:" << endl << e << endl;
    ACE_OS::exit(1);
  }
}

void sutra_rtc_brama::publish() {
  ACE_Time_Value ace_wait(0, 25);


  int nvalid = 0;
  for (unsigned int i = 0; i < d_centro.size(); i++) {
    nvalid += d_centro[i]->nvalid;
  }

  int ncmd = 0;
  for (unsigned int i = 0; i < d_control.size(); i++) {
    ncmd += d_control[i]->nactu();
  }

  buff_intensities = new CORBA::Octet[nvalid * sizeof(float)];
  buff_slopes = new CORBA::Octet[nvalid * 2 * sizeof(float)];
  buff_commands = new CORBA::Octet[ncmd * sizeof(float)];

  dims_intensities = new CORBA::ULong[2];
  dims_intensities[0] = 1;
  dims_intensities[1] = nvalid;

  dims_slopes = new CORBA::ULong[2];
  dims_slopes[0] = 1;
  dims_slopes[1] = nvalid * 2;

  dims_commands = new CORBA::ULong[2];
  dims_commands[0] = 1;
  dims_commands[1] = ncmd;

  //cerr << nvalid << " " << ncmd << endl;

  BRAMA::SuperFrame zFrame;
  zFrame.from = CORBA::string_dup("BRAMA RTC");
  zFrame.framecounter = framecounter;

  d_centro[0]->wfs->d_subsum->device2host((float*) buff_intensities);
  zFrame.ints.from = CORBA::string_dup("BRAMA intensities");
  zFrame.ints.typeofelements = CORBA::string_dup("intensities");
  zFrame.ints.sizeofelements = sizeof(float);
  zFrame.ints.dimensions = BRAMA::Dims(2, 2, dims_intensities, 0);
  zFrame.ints.framecounter = framecounter;
  zFrame.ints.data = BRAMA::Values(dims_intensities[1], dims_intensities[1],
      buff_intensities, 0);
  zFrame.ints.timestamp = get_timestamp();

  d_control[0]->d_centroids->device2host((float*) buff_slopes);
  zFrame.slps.from = CORBA::string_dup("BRAMA slopes");
  zFrame.slps.typeofelements = CORBA::string_dup("slopes");
  zFrame.slps.sizeofelements = sizeof(float);
  zFrame.slps.dimensions = BRAMA::Dims(2, 2, dims_slopes, 0);
  zFrame.slps.framecounter = framecounter;
  zFrame.slps.data = BRAMA::Values(dims_slopes[1], dims_slopes[1], buff_slopes,
      0);
  zFrame.slps.timestamp = get_timestamp();

  d_control[0]->d_com->device2host((float*) buff_commands);
  BRAMA::Frame zCommands;
  zFrame.cmds.from = CORBA::string_dup("BRAMA commands");
  zFrame.cmds.typeofelements = CORBA::string_dup("commands");
  zFrame.cmds.sizeofelements = sizeof(float);
  zFrame.cmds.dimensions = BRAMA::Dims(2, 2, dims_commands, 0);
  zFrame.cmds.framecounter = framecounter;
  zFrame.cmds.data = BRAMA::Values(dims_commands[1], dims_commands[1],
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
  //ACE_OS::sleep(ace_wait);

  delete[] buff_intensities;
  delete[] buff_slopes;
  delete[] buff_commands;

  delete[] dims_intensities;
  delete[] dims_slopes;
  delete[] dims_commands;
}

#endif /* USE_BRAMA */
