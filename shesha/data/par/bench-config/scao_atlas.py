import shesha.config as conf


class Param_camera:

    def __init__(self) -> None:
        self.type = "FakeCam"  # "FakeCam", "Vimbacam", "VimbacamBrahma"  or "EVTcam"
        self.camAddr = ""
        self.width = 0
        self.height = 0
        self.offset_w = 0
        self.offset_h = 0
        self.expo_usec = 0
        self.framerate = 0


# Camera params
# TODO set CameraType
p_cams = [Param_camera()]
p_cams[0].type = "VimbacamBrahma"
p_cams[0].camAddr = "169.254.1.14"
p_cams[0].width = 1936
p_cams[0].height = 1216
p_cams[0].offset_w = 0
p_cams[0].offset_h = 0
p_cams[0].expo_usec = 1000
p_cams[0].framerate = 100

# WFS params
p_wfss = [conf.ParamWfs()]
p_wfss[0].set_type("sh")
p_wfss[0].set_nxsub(16)
p_wfss[0].set_npix(8)

# p_wfss[0].set_validsubsx( validsubsx )
# for SH p_nvalid.shape(Nsubap) lowerleft corner for each ssp
# 	[X0, X1, ..., XN]
# for PYR p_nvalid.shape(4*nPixPerPupil)
#   [P0X0, P0X1, ..., P0XN, P1X0, ..., P4XN]
# for PYRGEN p_nvalid.shape(nPixTot)
#   [X0, X1, ..., XN]

# p_wfss[0].set_validsubsy( validsubsx )
# for SH p_nvalid.shape(Nsubap) lowerleft corner for each ssp
# 	[Y0, Y1, ..., YN]
# for PYR p_nvalid.shape(4*nPixPerPupil)
#   [P0Y0, P0Y1, ..., P0YN, P1Y0, ..., P4YN]
# for PYRGEN p_nvalid.shape(nPixTot)
#   [Y0, Y1, ..., YN]

#DM params
p_dms = [conf.ParamDm()]
p_dms[0].set_type("pzt")
p_dms[0].set_alt(0.)
p_dms[0].set_nact(4096)

# Centroiders params
p_centroiders = [conf.ParamCentroider()]
p_centroiders[0].set_nwfs(0)
p_centroiders[0].set_type("cog")
# p_centroiders[0].set_type("corr")
# p_centroiders[0].set_type_fct("model")

# Controller params
p_controllers = [conf.ParamController()]
p_controllers[0].set_type("ls")
p_controllers[0].set_nwfs([0])
p_controllers[0].set_ndm([0, 1])
p_controllers[0].set_maxcond(100.)
p_controllers[0].set_delay(1)
p_controllers[0].set_gain(0.4)
