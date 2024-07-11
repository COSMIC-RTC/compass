## @package   shesha
## @brief     Shesha package
## @author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
## @version   5.5.0
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
''' @package shesha
Documentation for shesha.

More details.
'''

import sys
import subprocess
import rich
from packaging import version

__version__ = "5.5.0"

__api_version__ = "5.5.0"

def check_shesha_compass_versions():
    compass_package = subprocess.check_output('conda list compass | tail -n1',shell=True).decode(
            sys.stdout.encoding)
    if(compass_package.startswith('compass')):
        # using conda package version
        compass_version = version.parse(compass_package.split()[1])
        
        # using shesha package version
        shesha_version = version.parse(__api_version__)
        
        if compass_version < shesha_version:
            rich.print(f'[bold red]WARNING: [/bold red] SHESHA version ({shesha_version}) is higher than COMPASS version ({compass_version}).')
            rich.print('[bold red]WARNING: [/bold red] Please update COMPASS: please run "conda update -c compass compass"')
        elif compass_version > shesha_version:
            rich.print(f'[bold red]WARNING: [/bold red] COMPASS version ({compass_version}) is higher than SHESHA version ({shesha_version}).')
            rich.print('[bold red]WARNING: [/bold red] Please update SHESHA: please run "git pull origin main"')

        assert(shesha_version == compass_version), f'SHESHA and COMPASS versions are not matching : {shesha_version} != {compass_version}'

check_shesha_compass_versions()
