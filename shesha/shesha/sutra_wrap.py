## @package   shesha
## @brief     Shesha package
## @author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
## @date      2022/01/24
## @copyright 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>
#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>

# COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
# General Public License as published by the Free Software Foundation, either version 3 of the 
# License, or any later version.

# COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
# If not, see <https://www.gnu.org/licenses/>

# Copyright (C) 2011-2024 COSMIC Team <https//://github.com/COSMIC-RTC/compass>
import importlib


def smart_import(mod, cls, verbose=False, silent=False):
    try:
        if verbose:
            print("trying from " + mod + " import " + cls)
        # my_module = __import__(mod, globals(), locals(), [cls], 0)
        my_module = importlib.import_module(mod)
        return getattr(my_module, cls)

    except ImportError as err:
        if not silent:
            import warnings
            warnings.warn(
                    "Error importing %s, it will be simulated due to: %s" %
                    (cls, err.msg), Warning)

        class tmp_cls:

            def __init__(self, *args, **kwargs):
                raise RuntimeError("Can not initilize the simulation with fake objects")

        return tmp_cls

    except AttributeError as err:
        if not silent:
            import warnings
            warnings.warn(
                    "Error importing %s, it will be simulated due to: %s" %
                    (cls, err.args), Warning)

        class tmp_cls:

            def __init__(self, *args, **kwargs):
                raise RuntimeError("Can not initialize the simulation with fake objects")

        return tmp_cls


# The import of carma MUST be done first
# since MAGMA >= 2.5.0
# Otherwise, it causes huge memory leak
# plus not working code
# Why ? We don't know... TB check with further version of MAGMA
carma_context = smart_import("carma", "context")

Dms = smart_import("sutra", "Dms")
Rtc_FFF = smart_import("sutra", "Rtc_FFF")
Rtc_FHF = smart_import("sutra", "Rtc_FHF", silent=True)
Rtc_UFF = smart_import("sutra", "Rtc_UFF", silent=True)
Rtc_UHF = smart_import("sutra", "Rtc_UHF", silent=True)
Rtc_FFU = smart_import("sutra", "Rtc_FFU", silent=True)
Rtc_FHU = smart_import("sutra", "Rtc_FHU", silent=True)
Rtc_UFU = smart_import("sutra", "Rtc_UFU", silent=True)
Rtc_UHU = smart_import("sutra", "Rtc_UHU", silent=True)
Sensors = smart_import("sutra", "Sensors")
Atmos = smart_import("sutra", "Atmos")
Telescope = smart_import("sutra", "Telescope")
Target = smart_import("sutra", "Target")
Gamora = smart_import("sutra", "Gamora")
Groot = smart_import("sutra", "Groot")
