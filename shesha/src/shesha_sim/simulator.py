import sys

from naga import naga_context

import shesha_init as init
import shesha_constants as scons

import Atmos, Telescope, Target, Rtc, Dms, Sensors
import time

from typing import Iterable, Callable, TypeVar, Any


class Simulator:

    def __init__(self, filepath: str=None) -> None:
        """
        TODO: docstring
        """
        self.is_init = False  # type: bool
        self.loaded = False  # type: bool
        self.config = None  # type: Any # types.ModuleType ?
        self.iter = 0  # type: int

        self.c = None  # type: naga_context
        self.atm = None  # type: Atmos.Atmos
        self.tel = None  # type: Telescope.Telescope
        self.tar = None  # type: Target.Target
        self.rtc = None  # type: Rtc.Rtc
        self.wfs = None  # type: Sensors.Sensors
        self.dms = None  # type: Dms.Dms

        if filepath is not None:
            self.load_from_file(filepath)

    def __str__(self) -> str:
        """
        TODO: docstring
        """
        s = ""
        if self.is_init:
            s += "====================\n"
            s += "Objects initialzed on GPU:\n"
            s += "--------------------------------------------------------\n"

            if self.atm is not None:
                s += self.atm.__str__() + '\n'
            if self.wfs is not None:
                s += self.wfs.__str__() + '\n'
            if self.dms is not None:
                s += self.dms.__str__() + '\n'
            if self.tar is not None:
                s += self.tar.__str__() + '\n'
            if self.rtc is not None:
                s += self.rtc.__str__() + '\n'
        else:
            s += "Simulator is not initialized."

        return s

    def load_from_file(self, filepath: str) -> None:
        """
        TODO: docstring
        """
        self.loaded = False
        self.is_init = False
        filename = filepath.split('/')[-1]
        if (len(filepath.split('.')) > 1 and filepath.split('.')[-1] != "py"):
            raise ValueError("Config file must be .py")

        pathfile = filepath.split(filename)[0]
        if (pathfile not in sys.path):
            sys.path.insert(0, pathfile)

        if self.config is not None:
            name = self.config.__name__
            print("Removing previous config")
            self.config = None
            try:
                del sys.modules[name]
            except:
                pass

        print(("loading ", filename.split(".py")[0]))
        self.config = __import__(filename.split(".py")[0])
        # exec("import %s as wao_config" % filename.split(".py")[0])
        sys.path.remove(pathfile)

        # Set missing config attributes to None
        if not hasattr(self.config, 'p_loop'):
            self.config.p_loop = None
        if not hasattr(self.config, 'p_geom'):
            self.config.p_geom = None
        if not hasattr(self.config, 'p_tel'):
            self.config.p_tel = None
        if not hasattr(self.config, 'p_atmos'):
            self.config.p_atmos = None
        if not hasattr(self.config, 'p_dms'):
            self.config.p_dms = None
        if not hasattr(self.config, 'p_target'):
            self.config.p_target = None
        if not hasattr(self.config, 'p_wfss'):
            self.config.p_wfss = None
        if not hasattr(self.config, 'p_centroiders'):
            self.config.p_tel = None
        if not hasattr(self.config, 'p_controllers'):
            self.config.p_tel = None

        self.loaded = True

    def clear_init(self) -> None:
        if self.loaded and self.is_init:
            self.iter = 0

            self.c = None
            self.atm = None
            self.tel = None
            self.tar = None
            self.rtc = None
            self.wfs = None
            self.dms = None
            self.is_init = False

    def init_sim(self) -> None:
        """
        TODO: docstring
        """
        if not self.loaded:
            raise ValueError("Config must be loaded before call to init_sim")
        self.c = naga_context(devices=self.config.p_loop.devices)

        if self.config.p_tel is None or self.config.p_geom is None:
            raise ValueError(
                    "Telescope geometry must be defined (p_geom and p_tel)")

        if self.config.p_atmos is not None:
            r0 = self.config.p_atmos.r0
        else:
            raise ValueError('A r0 value through a Param_atmos is required.')

        if self.config.p_loop is not None:
            ittime = self.config.p_loop.ittime
        else:
            raise ValueError(
                    'An ittime (iteration time in seconds) value through a Param_loop is required.'
            )

        print("->tel")
        self.tel = init.tel_init(
                self.c, self.config.p_geom, self.config.p_tel, r0, ittime,
                self.config.p_wfss)

        if self.config.p_atmos is not None:
            #   atmos
            print("->atmos")
            self.atm = init.atmos_init(
                    self.c, self.config.p_atmos, self.config.p_tel,
                    self.config.p_geom, ittime)
        else:
            self.atm = None

        if self.config.p_dms is not None:
            #   dm
            print("->dm")
            self.dms = init.dm_init(
                    self.c, self.config.p_dms, self.config.p_tel,
                    self.config.p_geom, self.config.p_wfss)
        else:
            self.dms = None

        self._tar_init()

        if self.config.p_wfss is not None:
            print("->wfs")
            self.wfs = init.wfs_init(
                    self.c, self.tel, self.config.p_wfss, self.config.p_tel,
                    self.config.p_geom, self.config.p_dms, self.config.p_atmos)
        else:
            self.wfs = None

        self._rtc_init(ittime)

        self.is_init = True

    def _tar_init(self) -> None:
        if self.config.p_target is not None:
            print("->target")
            self.tar = init.target_init(
                    self.c,
                    self.tel,
                    self.config.p_target,
                    self.config.p_atmos,
                    self.config.p_tel,
                    self.config.p_geom,
                    self.config.p_dms,
                    brama=False)
        else:
            self.tar = None

    def _rtc_init(self, ittime: float) -> None:
        if self.config.p_controllers is not None or self.config.p_centroiders is not None:
            print("->rtc")
            #   rtc
            self.rtc = init.rtc_init(
                    self.c,
                    self.tel,
                    self.wfs,
                    self.dms,
                    self.atm,
                    self.config.p_wfss,
                    self.config.p_tel,
                    self.config.p_geom,
                    self.config.p_atmos,
                    ittime,
                    self.config.p_centroiders,
                    self.config.p_controllers,
                    self.config.p_dms,
                    brama=False)
        else:
            self.rtc = None

    def next(
            self,
            *,
            move_atmos: bool=True,
            nControl: int=0,
            tar_trace: Iterable[int]=None,
            wfs_trace: Iterable[int]=None,
            apply_control: bool=True) -> None:
        '''
            function next
            Iterates the AO loop, with optional parameters

        :param move_atmos (bool): move the atmosphere for this iteration, default: True
        :param nControl (int): Controller number to use, default 0 (single control configurations)
        :param tar_trace (None or list[int]): list of targets to trace. None equivalent to all.
        :param wfs_trace (None or list[int]): list of WFS to trace. None equivalent to all.
        :param apply_control (bool): (optional) if True (default), apply control on DMs
        '''
        if tar_trace is None:
            tar_trace = range(self.config.p_target.ntargets)
        if wfs_trace is None:
            wfs_trace = range(len(self.config.p_wfss))

        if move_atmos:
            self.atm.move_atmos()
        if (
                self.config.p_controllers[nControl].type_control ==
                scons.ControllerType.GEO):
            for t in tar_trace:
                self.tar.raytrace(t, b"atmos", self.tel, self.atm)
                self.rtc.do_control_geo(nControl, self.dms, self.tar, t)
                self.rtc.apply_control(nControl, self.dms)
                self.tar.raytrace(t, b"dm", self.tel, dms=self.dms)
        else:
            for t in tar_trace:
                self.tar.raytrace(t, b"all", self.tel, self.atm, self.dms)
            for w in wfs_trace:
                self.wfs.raytrace(w, b"all", self.tel, self.atm, self.dms)
                self.wfs.comp_img(w)
            self.rtc.do_centroids(nControl)
            self.rtc.do_control(nControl)
            self.rtc.do_clipping(0, -1e5, 1e5)
            if apply_control:
                self.rtc.apply_control(nControl, self.dms)
        self.iter += 1

    def loop(self, n=1, monitoring_freq=100, **kwargs):
        """
        TODO: docstring
        """
        print("----------------------------------------------------")
        print("iter# | S.E. SR | L.E. SR | ETR (s) | Framerate (Hz)")
        print("----------------------------------------------------")
        t0 = time.time()
        for i in range(n):
            self.next(**kwargs)
            if ((i + 1) % monitoring_freq == 0):
                framerate = (i + 1) / (time.time() - t0)
                strehltmp = self.tar.get_strehl(0)
                etr = (n - i) / framerate
                print(
                        "%d \t %.3f \t  %.3f\t     %.1f \t %.1f" %
                        (i + 1, strehltmp[0], strehltmp[1], etr, framerate))
        t1 = time.time()
        print(
                " loop execution time:", t1 - t0, "  (", n, "iterations), ",
                (t1 - t0) / n, "(mean)  ", n / (t1 - t0), "Hz")


class SimulatorBrama(Simulator):
    """
        Class SimulatorBrama: Brama overloaded simulator
        _tar_init and _rtc_init to instantiate Brama classes instead of regular classes
        next() to call rtc/tar.publish()
    """

    def _tar_init(self) -> None:
        '''
            TODO
        '''
        if self.config.p_target is not None:
            print("->target")
            self.tar = init.target_init(
                    self.c,
                    self.tel,
                    self.config.p_target,
                    self.config.p_atmos,
                    self.config.p_tel,
                    self.config.p_geom,
                    self.config.p_dms,
                    brama=True)
        else:
            self.tar = None

    def _rtc_init(self, ittime) -> None:
        '''
            TODO
        '''
        if self.config.p_controllers is not None or self.config.p_centroiders is not None:
            print("->rtc")
            #   rtc
            self.rtc = init.rtc_init(
                    self.c,
                    self.tel,
                    self.wfs,
                    self.dms,
                    self.atm,
                    self.config.p_wfss,
                    self.config.p_tel,
                    self.config.p_geom,
                    self.config.p_atmos,
                    ittime,
                    self.config.p_centroiders,
                    self.config.p_controllers,
                    self.config.p_dms,
                    brama=True)
        else:
            self.rtc = None

    def next(self, **kwargs) -> None:
        Simulator.next(self, **kwargs)
        if self.rtc is not None:
            self.rtc.publish()
        if self.tar is not None:
            self.tar.publish()


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
        print(
                '** Execution time of {}: {} seconds'.format(
                        function.__name__, t2 - t1))
        return ret

    return new_func


class Bench(Simulator):
    '''
        Class Bench

        Timed version of the simulator class using decorated overloads
    '''

    @timeit
    def __init__(self, filepath: str=None) -> None:
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
    def loop(self, n: int=1, monitoring_freq: int=100, **kwargs) -> None:
        Simulator.loop(self, n=n, monitoring_freq=monitoring_freq, **kwargs)
