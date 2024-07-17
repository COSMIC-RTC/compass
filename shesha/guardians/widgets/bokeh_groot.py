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
"""
To launch it :

    - locally :
        bokeh serve --show bokeh_display.py
    - as a server :
        bokeh serve --port 8081 --allow-websocket-origin hippo6.obspm.fr:8081 bokeh_roket.py
        then, open a web browser and connect to http://hippo6.obspm.fr:8081/bokeh_roket
"""

from widget_groot import Bokeh_groot
from bokeh.io import curdoc
import glob
import os
import atexit


def remove_files():
    files = glob.glob("/home/fferreira/public_html/roket_display*")
    for f in files:
        os.remove(f)


widget = Bokeh_groot()
curdoc().clear()
# widget.update()
# output_file("roket.html")
# show(widget.tab)
curdoc().add_root(widget.tab)

atexit.register(remove_files)
