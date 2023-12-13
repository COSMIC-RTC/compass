import shesha.constants as scons
import shesha.config as conf

# FOR BOTH SIMU AND RTC STANDALONE
simul_name = "PYRCADO_39m_20190122"

p_loop = conf.ParamLoop()
p_loop.set_niter(1000)
p_loop.set_ittime(1 / 500.)  # =1/500
p_loop.set_devices([4, 5])


# For the RTC standalone only
# fakecam to pass to bench
class Param_camera:
    pass


p_cams = [Param_camera()]
p_cams[0].type = 'FakeCam'
p_cams[0].camAddr = ""
p_cams[0].width = 800
p_cams[0].height = 800
p_cams[0].offset_w = 0
p_cams[0].offset_h = 0
p_cams[0].expo_usec = 1000
p_cams[0].framerate = 100

# centroider
p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)
p_centroider0.set_type(scons.CentroiderType.MASKEDPIX)

# controller
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

p_controller0.set_type(scons.ControllerType.GENERIC)
p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(150.)
p_controller0.set_delay(2.0)
p_controller0.set_gain(0.4)
