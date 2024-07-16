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
"""@package shesha
Documentation for shesha.

More details.
"""

import rich
from packaging import version
import carma
import sutra

__version__ = "6.0.0"


def check_shesha_compass_versions():
    if carma.__version__ != sutra.__version__:
        rich.print(
            f"[bold red]WARNING: [/bold red] CARMA version ({carma.__version__}) is different from SUTRA version ({sutra.__version__})."
        )
        rich.print("[bold red]WARNING: [/bold red] Please recompile COMPASS")

    compass_version = version.parse(sutra.__version__)
    shesha_version = version.parse(__version__)

    if compass_version != shesha_version:
        rich.print(
            f"[bold red]WARNING: [/bold red] SHESHA version ({shesha_version}) is different from COMPASS version ({compass_version})."
        )
        rich.print("[bold red]WARNING: [/bold red] Please recompile COMPASS")

    assert (
        shesha_version == compass_version
    ), f"SHESHA and COMPASS versions are not matching : {shesha_version} != {compass_version}"


check_shesha_compass_versions()
