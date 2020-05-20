## @package   shesha.sim.bench
## @brief     Benchmark class for COMPASS simulation timing (Not used, incomplete)
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   4.3.0
## @date      2020/05/18
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

import time

from typing import Callable, TypeVar

from .simulator import Simulator

_O = TypeVar('_O')


def timeit(function: Callable[..., _O]) -> Callable[..., _O]:
    '''
        Function timing decorator
    '''

    def new_func(*args, **kwargs) -> _O:
        print('** Timed call to function {}'.format(function.__name__))
        t1 = time.time()
        ret = function(*args, **kwargs)
        t2 = time.time()
        print('** Execution time of {}: {} seconds'.format(function.__name__, t2 - t1))
        return ret

    return new_func


class Bench(Simulator):
    '''
        Class Bench

        Timed version of the simulator class using decorated overloads
    '''

    @timeit
    def __init__(self, filepath: str = None) -> None:
        Simulator.__init__(self, filepath)

    @timeit
    def load_from_file(self, filepath: str) -> None:
        Simulator.load_from_file(self, filepath)

    @timeit
    def init_sim(self) -> None:
        Simulator.init_sim(self)

    @timeit
    def timed_next(self) -> None:
        Simulator.next(self)

    @timeit
    def loop(self, n: int = 1, monitoring_freq: int = 100, **kwargs) -> None:
        Simulator.loop(self, n=n, monitoring_freq=monitoring_freq, **kwargs)
