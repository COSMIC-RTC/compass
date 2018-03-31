import shesha_config as conf

class Param_camera:
    def __init__(self) -> None:
        self.camAddr = ""
        self.width = 0
        self.height = 0
        self.offset_w = 0
        self.offset_h = 0
        self.expo_usec = 0
        self.framerate = 0
        
p_cams = [Param_camera()]
p_cams[0].camAddr = ""
p_cams[0].width = 512
p_cams[0].height = 512
p_cams[0].offset_w = 0
p_cams[0].offset_h = 0
p_cams[0].expo_usec = 1000
p_cams[0].framerate = 100

p_wfss= [conf.Param_wfs()]
p_wfss[0].set_type("sh")
p_wfss[0].set_nxsub(16)
p_wfss[0].set_npix(16)
