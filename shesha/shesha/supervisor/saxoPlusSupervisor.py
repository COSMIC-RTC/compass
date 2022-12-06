## @package   shesha.supervisor.canapassSupervisor
## @brief     Initialization and execution of a CANAPASS supervisor
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   5.2.1
## @date      2022/01/24
## @copyright GNU Lesser General Public License
#
#  This file is part of COMPASS <https://anr-compass.github.io/compass/>
#
#  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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
"""
Initialization and execution of a saxo+Supervisor manager. 
It instanciates in one process two compass simulations: 
1 for first stage and 1 for the second stage (as defined in their relative .par files)

IMPORTANT:
The next method of this manager --superseeds-- the compass next method so that the loop is fully handled by the saxoPlus manager. 

Usage:
  saxo+Supervisor.py <saxoparameters_filename> <saxoPlusparameters_filename> [options]

with 'saxoparameters_filename' the path to the parameters file for SAXO+ First stage (I.e current SAXO system)
with 'saxoPlusparameters_filename' the path to the parameters file for SAXO+ Second stage

Options:
  -a, --adopt       used to connect ADOPT client to the manager (via pyro + shm cacao)

Example: 
    ipython -i saxoPlusSupervisor.py ../../data/par/SPHERE+/sphere.py ../../data/par/SPHERE+/sphere+.py
    ipython -i saxoPlusSupervisor.py ../../data/par/SPHERE+/sphere.py ../../data/par/SPHERE+/sphere+.py -- --adopt
"""

import os, sys
import numpy as np
import time

from tqdm import tqdm
import astropy.io.fits as pfits
from threading import Thread
from subprocess import Popen, PIPE

import shesha.ao as ao
import shesha.constants as scons
from shesha.constants import CentroiderType, WFSType

from typing import Any, Dict, Tuple, Callable, List
from shesha.supervisor.compassSupervisor import CompassSupervisor



class SaxoPlusManager():
    """
        DOC To be done..
    """
    def __init__(self, first_stage, second_stage):
        """ 
        Init of the saxoPlusManager object
        """

        self.first_stage = first_stage
        self.second_stage = second_stage
        self.second_stage.atmos.enable_atmos(False) # second stage atmos is not used 
        self.iterations = 0
        mpup_shape = self.second_stage.config.p_geom._mpupil.shape
        self.second_stage_input = np.zeros((mpup_shape[0], mpup_shape[1], 1))
        residual_shape = self.first_stage.config.p_geom._spupil.shape
        self.mpup_offset = (mpup_shape[0] - residual_shape[0]) // 2
        self.frequency_ratio = round(1e-3 / self.second_stage.config.p_loop.ittime)

    def next(self, seeAtmos=True):
        """
        MAIN method that allows to manage properly the 2 AO stages of SAXO+ system. 
        The phase residuals (including turbulence + AO loop residuals) of the first stage simulation is sent to second stage simulation
        at each iteration of the manager. 
        The saxo+ manager disable the seconds stage turbulence simulation (as it is propageated through the first stage residals if any). 

        This next method sould ALWAYS be called to perform a regular SAXO+ simulation 
        instead of the individuals COMPASS next methods to ensure the correct synchronisation of the 2 systems. 
        """
        # Iteration time of the first stage is set as the same as the second stage to allow
        # correct atmosphere movement for second stage integration. Then, first stage as to compute
        # every 3 iteration to be 3 times slower than the second stage
        # self.first_stage.atmos.enable_atmos(seeAtmos) #Enabling (or not) Turbulence
        self.second_stage.atmos.enable_atmos(False) # Turbulence always disabled on 2nd instance of COMPASS
        if not (self.iterations % self.frequency_ratio): # Time for first stage full computation: We update the first stage command every frequency_ratio iterations. 
            self.first_stage.next()
        else: # Only raytracing current tubulence phase (if any) and current DMs phase (no command updates). 
            self.first_stage.next(do_control=False, apply_control=False, compute_tar_psf=False)
        # Get residual of first stage to put it into second stage
        # For now, involves GPU-CPU memory copies, can be improved later if speed is a limiting factor here... 
        first_stage_residual = self.first_stage.target.get_tar_phase(0)
        self.second_stage_input[self.mpup_offset:-self.mpup_offset,self.mpup_offset:-self.mpup_offset,:] = first_stage_residual[:,:,None]
        self.second_stage.tel.set_input_phase(self.second_stage_input) # 1st stage residuals sent to seconds stage simulation. 
        # Second stage computation
        self.second_stage.next(move_atmos=False) #"Updates the second stage siulation accordingly". 
        self.iterations += 1


class loopHandler:

    def __init__(self):
        pass

    def start(self):
        pass

    def stop(self):
        pass

    def alive(self):
        return "alive"

if __name__ == '__main__':
    from docopt import docopt
    from shesha.config import ParamConfig
    arguments = docopt(__doc__)
    adopt = arguments["--adopt"]

    config1 = ParamConfig(arguments["<saxoparameters_filename>"])
    config2 = ParamConfig(arguments["<saxoPlusparameters_filename>"])

    """
    if (arguments["--freq"]):
        print("Warning changed frequency loop to: ", arguments["--freq"])
        config.p_loop.set_ittime(1 / float(arguments["--freq"]))
    if (arguments["--delay"]):
        print("Warning changed delay loop to: ", arguments["--delay"])
        config.p_controllers[0].set_delay(float(arguments["--delay"]))
    if (arguments["--spiders"]):
        print("Warning changed spiders size to: ", arguments["--spiders"])
        config.p_tel.set_t_spiders(float(arguments["--spiders"]))
    if (arguments["--nxsub"]):
        print("Warning changed number of pixels per subaperture to: ", arguments["--nxsub"])
        config.p_wfss[0].set_nxsub(int(arguments["--nxsub"]))
    if (arguments["--pupsep"]):
        print("Warning changed distance between subaperture center and frame center to: ", arguments["--pupsep"])
        config.p_wfss[0].set_pyr_pup_sep(int(arguments["--pupsep"]))
    if (arguments["--gsmag"]):
        print("Warning changed guide star magnitude to: ", arguments["--gsmag"])
        config.p_wfss[0].set_gsmag(float(arguments["--gsmag"]))
    if (arguments["--setr0"]):
        print("Warning changed r0 to: ", arguments["--setr0"])
        config.p_atmos.set_r0(float(arguments["--setr0"]))
    if (arguments["--rmod"]):
        print("Warning changed modulation radius to: ", arguments["--rmod"])
        rMod = int(arguments["--rmod"])
        nbPtMod = int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.) * 4)
        config.p_wfss[0].set_pyr_npts(nbPtMod)
        config.p_wfss[0].set_pyr_ampl(rMod)
    if (arguments["--offaxis"]):
        print("Warning changed target x position: ", arguments["--offaxis"])
        config.p_targets[0].set_xpos(float(arguments["--offaxis"]))
        config.p_targets[1].set_xpos(float(arguments["--offaxis"]))
        config.p_targets[2].set_xpos(float(arguments["--offaxis"]))
    """


    first_stage = CompassSupervisor(config1, cacao=adopt)
    second_stage = CompassSupervisor(config2, cacao=adopt)
    manager = SaxoPlusManager(first_stage, second_stage)

    if(adopt): 
        
        supervisor1 = manager.first_stage
        supervisor2 = manager.second_stage

        
        try:
            from subprocess import Popen, PIPE
            from hraa.server.pyroServer import PyroServer
            import Pyro4
            Pyro4.config.REQUIRE_EXPOSE = False
            p = Popen("whoami", shell=True, stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            if (err != b''):
                print(err)
                raise Exception("ERROR CANNOT RECOGNIZE USER")
            else:
                user = out.split(b"\n")[0].decode("utf-8")
                print("User is " + user)


            devices1 = [
                    supervisor1, supervisor1.rtc, supervisor1.wfs, supervisor1.target,
                    supervisor1.tel, supervisor1.basis, supervisor1.calibration,
                    supervisor1.atmos, supervisor1.dms, supervisor1.config, supervisor1.modalgains
            ]
            devices2 = [
                    supervisor2, supervisor2.rtc, supervisor2.wfs, supervisor2.target,
                    supervisor2.tel, supervisor2.basis, supervisor2.calibration,
                    supervisor2.atmos, supervisor2.dms, supervisor2.config, supervisor2.modalgains
            ]
            names = [
                    "supervisor", "supervisor_rtc", "supervisor_wfs", "supervisor_target",
                    "supervisor_tel", "supervisor_basis", "supervisor_calibration",
                    "supervisor_atmos", "supervisor_dms", "supervisor_config", "supervisor_modalgains"
            ]

            label = "firstStage"
            nname = []
            for name in names:
                nname.append(name + "_" + user + "_" +label)

            label = "secondStage"
            for name in names:
                nname.append(name + "_" + user + "_" +label)

            nname.append('supervisorSAXOPlus'+ "_" + user ) # Adding master next dedicated to trigger SAXO+ hybrid loop
            devices = devices1 + devices2 + [manager]
            server = PyroServer(listDevices=devices, listNames=nname)
            #server.add_device(supervisor, "waoconfig_" + user)
            server.start()
        except:
            raise EnvironmentError(
                    "Missing dependencies (code HRAA or Pyro4 or Dill Serializer)")
        