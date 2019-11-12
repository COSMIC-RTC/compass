## @package   shesha.sim.benchBrahma
## @brief     Benchmark class for COMPASS with BRAHMA simulation timing (Not used, incomplete)
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   4.3.0
## @date      2011/01/28
## @copyright GNU Lesser General Public License
#
#  This file is part of COMPASS <https://anr-compass.github.io/compass/>
#
#  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
#  All rights reserved.
#  Distributed under GNU - LGPL
#
#  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
#  General Public License as published by the Free Software Foundation, either version 3 of the License,
#  or any later version.
#
#  COMPASS: End-to-end AO simulation tool using GPU acceleration
#  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
#
#  The final product includes a software package for simulating all the critical subcomponents of AO,
#  particularly in the context of the ELT and a real-time core based on several control approaches,
#  with performances consistent with its integration into an instrument. Taking advantage of the specific
#  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
#  conduct large simulation campaigns called to the ELT.
#
#  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
#  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
#  various systems configurations such as multi-conjugate AO.
#
#  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
#  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
#  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.

from time import time, sleep

from .simulatorBrahma import SimulatorBrahma


class BenchBrahma(SimulatorBrahma):
    '''
        Class BenchBrahma
    '''

    def next(self, *, nControl: int = 0) -> None:
        '''
            function next
            Iterates on centroiding and control, with optional parameters

        :param nControl (int): Controller number to use, default 0 (single control configurations)
        '''

        # self.rtc.do_centroids(nControl)
        # self.rtc.do_control(nControl)
        # self.rtc.do_clipping(0)

        if self.rtc is not None:
            self.rtc.publish()
        # if self.tar is not None:
        #     self.tar.publish()

    def loop(self, freq=0., monitoring_freq=100, **kwargs):
        """
        TODO: docstring
        """
        print("----------------------------------------------------")
        print("iter# | S.E. SR | L.E. SR | ETR (s) | Framerate (Hz)")
        print("----------------------------------------------------")
        niter = 0
        t0 = time()
        while True:
            start = time()
            self.next(**kwargs)
            if ((niter + 1) % monitoring_freq == 0):
                framerate = (monitoring_freq) / (time() - t0)
                t0 = time()
                print("%d \t %.1f" % (niter + 1, framerate))
            tmp = time() - start
            if freq > 0:
                sleep(1 / freq - tmp)
            niter += 1
