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


import shesha.config as conf

from shesha.constants import CONST

import numpy as np
from shesha.sutra_wrap import carma_context, Target, Telescope


def target_init(
    ctxt: carma_context,
    telescope: Telescope,
    p_targets: list,
    p_atmos: conf.ParamAtmos,
    p_tel: conf.ParamTel,
    p_geom: conf.ParamGeom,
    dm=None,
):
    """Create a cython target from parametres structures

    Args:
        ctxt: (carma_context) :
        telescope: (Telescope): Telescope object
        p_targets: (lis of ParamTarget) : target_settings
        p_atmos: (ParamAtmos) : atmos settings
        p_tel: (ParamTel) : telescope settings
        p_geom: (ParamGeom) : geom settings
        dm: (ParamDm) : (optional) dm settings
    :return:
        tar: (Target): Target object
    """
    type_target = "atmos"

    if p_targets is not None:
        for p_target in p_targets:
            if p_target.dms_seen is None and dm is not None:
                p_target.dms_seen = np.arange(len(dm))

    sizes = np.ones(len(p_targets), dtype=np.int64) * p_geom.pupdiam

    ceiled_pupil = np.ceil(p_geom._spupil)

    ceiled_pupil[np.where(ceiled_pupil > 1)] = 1

    if p_target.apod == 1:
        Npts = 0
        # TODO apodizer, Npts=nb element of apodizer>0
        ceiled_apodizer = np.ceil(p_geom._apodizer * p_geom._spupil)
        ceiled_apodizer[np.where(ceiled_apodizer > 1)] = 1
        Npts = int(np.sum(ceiled_apodizer))
    else:
        Npts = int(np.sum(ceiled_pupil))

    xpos = np.array([p_target.xpos for p_target in p_targets], dtype=np.float32)
    ypos = np.array([p_target.ypos for p_target in p_targets], dtype=np.float32)
    Lambda = np.array([p_target.Lambda for p_target in p_targets], dtype=np.float32)
    mag = np.array([p_target.mag for p_target in p_targets], dtype=np.float32)
    zerop = p_targets[0].zerop

    target = Target(
        ctxt,
        telescope,
        len(p_targets),
        xpos,
        ypos,
        Lambda,
        mag,
        zerop,
        sizes,
        Npts,
        ctxt.active_device,
    )

    # cc=i
    for i in range(len(p_targets)):
        p_target = p_targets[i]
        if p_atmos.nscreens > 0:
            for j in range(p_atmos.nscreens):
                xoff = p_target.xpos * CONST.ARCSEC2RAD * p_atmos.alt[j] / p_atmos.pupixsize
                yoff = p_target.ypos * CONST.ARCSEC2RAD * p_atmos.alt[j] / p_atmos.pupixsize
                xoff += float((p_atmos.dim_screens[j] - p_geom._n) / 2)
                yoff += float((p_atmos.dim_screens[j] - p_geom._n) / 2)
                pupdiff = (p_geom._n - p_geom.pupdiam) / 2
                xoff += pupdiff
                yoff += pupdiff
                target.d_targets[i].add_layer(type_target, j, xoff, yoff)

        # if (y_dm != []) {
        if dm is not None:
            # j=ddd
            # for (ddd=1;ddd<=numberof(*y_target(cc).dms_seen);ddd++) {
            for j in range(p_target.dms_seen.size):
                # k=dd
                # dd = (*y_target(cc).dms_seen)(ddd)
                k = p_target.dms_seen[j]
                dims = dm[k]._n2 - dm[k]._n1 + 1
                dim = p_geom._mpupil[2].size
                dim_dm = max(dim, dims)
                xoff = p_target.xpos * CONST.ARCSEC2RAD * dm[k].alt / p_tel.diam * p_geom.pupdiam
                yoff = p_target.ypos * CONST.ARCSEC2RAD * dm[k].alt / p_tel.diam * p_geom.pupdiam

                xoff += float((dim_dm - p_geom._n) / 2)
                yoff += float((dim_dm - p_geom._n) / 2)

                pupdiff = (p_geom._n - p_geom.pupdiam) / 2
                xoff += pupdiff
                yoff += pupdiff

                # if (dm[k].type == scons.DmType.KL):
                #     xoff -= 2
                #     yoff -= 2
                target.d_targets[i].add_layer(dm[k].type, k, xoff, yoff)

        target.d_targets[i].init_strehlmeter()

    return target
