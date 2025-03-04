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


from shesha.supervisor.components.rtc.rtcCompass import RtcCompass
from shesha.supervisor.components.rtc.rtcStandalone import RtcStandalone

__all__ = ["RtcCompass", "RtcStandalone"]

try:
    from shesha.supervisor.components.rtc.rtcCosmic import RtcCosmic

    __all__ += ["RtcCosmic"]
except ImportError:
    pass
