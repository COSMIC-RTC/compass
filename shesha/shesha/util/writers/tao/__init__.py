from shesha.util.writers.tao.sysParams import *
from shesha.util.writers.tao.atmParams import *
from shesha.util.writers import common

def write_parfiles(sup, *, file_name_sys="./sysParams.json",
    file_name_atm="./atmParams.json", file_name_data="sys-input.fits", 
    wfss_indices=None, ts=False, dms_indices=None, imat_type="controller",
    controller_id=-1,influ_index=0):

    write_json_sys_param(sup, wfss_indices=wfss_indices, ts=ts,
        dms_indices=dms_indices,file_name=file_name_sys)

    write_json_atm_param(sup, file_name=file_name_atm)

    common.write_data(file_name_data, sup, wfss_indices=wfss_indices,
        dms_indices=dms_indices, controller_id=controller_id,
        influ=influ_index, compose_type="controller")