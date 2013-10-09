#!/bin/bash

/home/brujo/yorick-2.2/relocate/bin/yorick -batch ao_scripts.i relax0Gauss.par
/home/brujo/yorick-2.2/relocate/bin/yorick -batch ao_scripts.i relax10Gauss.par
/home/brujo/yorick-2.2/relocate/bin/yorick -batch ao_scripts.i relax20Gauss.par
/home/brujo/yorick-2.2/relocate/bin/yorick -batch ao_scripts.i relax30Gauss.par
/home/brujo/yorick-2.2/relocate/bin/yorick -batch ao_scripts.i relax40Gauss.par

/home/brujo/yorick-2.2/relocate/bin/yorick -batch ao_scripts.i relax0Exp.par
/home/brujo/yorick-2.2/relocate/bin/yorick -batch ao_scripts.i relax10Exp.par
/home/brujo/yorick-2.2/relocate/bin/yorick -batch ao_scripts.i relax20Exp.par
/home/brujo/yorick-2.2/relocate/bin/yorick -batch ao_scripts.i relax30Exp.par
/home/brujo/yorick-2.2/relocate/bin/yorick -batch ao_scripts.i relax40Exp.par
