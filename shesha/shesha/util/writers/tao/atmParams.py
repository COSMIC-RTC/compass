#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
#
# COMPASS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# COMPASS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with COMPASS. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2011-2024 COSMIC Team
import json
from shesha.util.writers import common


def write_json_atm_param(sup, *, file_name="./atm-params.json"):
    """Return a json representation of the atmospheric parameters

    Args:
        sup : (CompassSupervisor) : supervisor to get the json representation from
    """
    atm_json = {
        "notice": common.atmos_json_notice(),
        "profiles": [common.atmos_to_json(sup.config.p_atmos)],
    }
    f = open(file_name, "w")
    f.write(json.dumps(atm_json, indent=4))
    f.close()
