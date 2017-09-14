"""
To launch it :

    - locally :
        bokeh serve --show bokeh_display.py
    - as a server :
        bokeh serve --port 8081 --allow-websocket-origin hippo6.obspm.fr:8081 bokeh_roket.py
        then, open a web browser and connect to http://hippo6.obspm.fr:8081/bokeh_roket
"""
from widget_roket import Bokeh_roket
from bokeh.io import curdoc, output_file, show

widget = Bokeh_roket()
curdoc().clear()
widget.update()
#output_file("roket.html")
#show(widget.tab)
curdoc().add_root(widget.tab)
