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
from widget_roket import Bokeh_roket
from widget_gamora import Bokeh_gamora
from widget_groot import Bokeh_groot

from bokeh.models.widgets import Tabs


class Bokeh_guardian:
    """
    Class that defines a bokeh layout for all the guardians package
    Usage: see bokeh_roket.py which is the executable
    """

    def __init__(self):
        self.roket = Bokeh_roket()
        self.gamora = Bokeh_gamora()
        self.groot = Bokeh_groot()

        self.tab = Tabs(tabs=[self.roket.tab, self.gamora.tab, self.groot.tab])
