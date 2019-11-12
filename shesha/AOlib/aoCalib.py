import compassConfigToFile as cf


class adoptCalib_class():

    def __init__(self, config, wfs, tel, atm, dms, tar, rtc, ao):
        self.wfs = wfs
        self.tel = tel
        self.atm = atm
        self.dms = dms
        self.tar = tar
        self.rtc = rtc
        self.config = config
        self.ao = ao

    def getConfig(self, wao, path):
        return cf.returnConfigfromWao(wao, filepath=path)

    def returnkl2V(self):
        KL2V = self.ao.compute_KL2V(self.config.p_controllers[0], self.dms,
                                    self.config.p_dms, self.config.p_geom,
                                    self.config.p_atmos, self.config.p_tel)
        return KL2V

    def doRefslopes(self):
        self.rtc.do_centroids_ref(0)
        print("refslopes done")
