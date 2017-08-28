import h5py
import pandas
import os
import numpy as np
from subprocess import check_output

#import shesha


def params_dictionary(config):
    """ Create and returns a dictionary of all the config parameters with the
    corresponding keys for further creation of database and save files

    :param config: (module) : simulation parameters
    :return param_dict: (dictionary) : dictionary of parameters
    """

    commit = check_output(["git", "rev-parse", "--short", "HEAD"]).decode('utf8')

    param_dict = {
            "simul_name": config.simul_name.encode('utf8'),
            "commit": version.encode('utf8')
    }

    param_loop = [
            i for i in dir(config.p_loop)
            if (not i.startswith('_') and not i.startswith('set_') and
                not i.startswith('get_'))
    ]
    for k in param_loop:
        param_dict.update({
                "_Param_loop__" + k: config.p_loop.__dict__["_Param_loop__" + k]
        })
    param_geom = [
            i for i in dir(config.p_geom)
            if (not i.startswith('_') and not i.startswith('set_') and
                not i.startswith('get_'))
    ]
    for k in param_geom:
        param_dict.update({
                "_Param_geom__" + k: config.p_geom.__dict__["_Param_geom__" + k]
        })
    param_tel = [
            i for i in dir(config.p_tel)
            if (not i.startswith('_') and not i.startswith('set_') and
                not i.startswith('get_'))
    ]
    for k in param_tel:
        param_dict.update({
                "_Param_tel__" + k: config.p_tel.__dict__["_Param_tel__" + k]
        })
    if config.p_atmos is not None:
        param_atmos = [
                i for i in dir(config.p_atmos)
                if (not i.startswith('_') and not i.startswith('set_') and
                    not i.startswith('get_'))
        ]
        for k in param_atmos:
            param_dict.update({
                    "_Param_atmos__" + k: config.p_atmos.__dict__["_Param_atmos__" + k]
            })

    if config.p_target is not None:
        param_target = [
                i for i in dir(config.p_target)
                if (not i.startswith('_') and not i.startswith('set_') and
                    not i.startswith('get_'))
        ]
        for k in param_target:
            param_dict.update({
                    "_Param_target__" + k:
                            config.p_target.__dict__["_Param_target__" + k]
            })
    if config.p_wfss is not None:
        param_wfs = [
                i for i in dir(config.p_wfss[0])
                if (not i.startswith('_') and not i.startswith('set_') and
                    not i.startswith('get_'))
        ]
        for k in param_wfs:
            param_dict.update({
                    "_Param_wfs__" + k: [
                            w.__dict__["_Param_wfs__" + k] for w in config.p_wfss
                    ]
            })
    if config.p_dms is not None:
        param_dm = [
                i for i in dir(config.p_dms[0])
                if (not i.startswith('_') and not i.startswith('set_') and
                    not i.startswith('get_'))
        ]
        for k in param_dm:
            param_dict.update({
                    "_Param_dm__" + k: [
                            w.__dict__["_Param_dm__" + k] for w in config.p_dms
                    ]
            })
    if config.p_controllers is not None:
        param_controller = [
                i for i in dir(config.p_controllers[0])
                if (not i.startswith('_') and not i.startswith('set_') and
                    not i.startswith('get_'))
        ]
        for k in param_controller:
            param_dict.update({
                    "_Param_controller__" + k: [
                            w.__dict__["_Param_controller__" + k]
                            for w in config.p_controllers
                    ]
            })
    if config.p_centroiders is not None:
        param_centroider = [
                i for i in dir(config.p_centroiders[0])
                if (not i.startswith('_') and not i.startswith('set_') and
                    not i.startswith('get_'))
        ]
        for k in param_centroider:
            param_dict.update({
                    "_Param_centroider__" + k: [
                            w.__dict__["_Param_centroider__" + k]
                            for w in config.p_centroiders
                    ]
            })

    for k in param_dict.keys():
        if type(param_dict[k]) is list:
            param_dict[k] = [d if d is not None else -10 for d in param_dict[k]]
        elif param_dict[k] is None:
            param_dict[k] = -10
    return param_dict


def create_file_attributes(filename, param_dict):
    """ create_file_attributes(filename,config)
    Create an hdf5 file wtih attributes corresponding to all simulation parameters

    :param:
        filename : (str) : full path + filename to create
        config : () : simulation parameters
    """
    f = h5py.File(filename, "w")

    for i in list(param_dict.keys()):
        f.attrs.create(i, param_dict[i])
    f.attrs.create("validity", False)
    print(filename, "initialized")
    f.close()


def init_hdf5_files(savepath, param_dict, matricesToLoad):
    version = check_output(["git", "rev-parse", "--short", "HEAD"]).decode('utf8')
    # if not(matricesToLoad.has_key("A")):
    if "A" not in matricesToLoad:
        df = pandas.read_hdf(savepath + "matricesDataBase.h5", "A")
        ind = len(df.index)
        filename = savepath + "turbu/A_" + version + "_" + str(ind) + ".h5"
        create_file_attributes(filename, param_dict)
        updateDataBase(filename, savepath, "A")
        df = pandas.read_hdf(savepath + "matricesDataBase.h5", "B")
        ind = len(df.index)
        filename = savepath + "turbu/B_" + version + "_" + str(ind) + ".h5"
        create_file_attributes(filename, param_dict)
        updateDataBase(filename, savepath, "B")
        df = pandas.read_hdf(savepath + "matricesDataBase.h5", "istx")
        ind = len(df.index)
        filename = savepath + "turbu/istx_" + version + "_" + str(ind) + ".h5"
        create_file_attributes(filename, param_dict)
        updateDataBase(filename, savepath, "istx")
        df = pandas.read_hdf(savepath + "matricesDataBase.h5", "isty")
        ind = len(df.index)
        filename = savepath + "turbu/isty_" + version + "_" + str(ind) + ".h5"
        create_file_attributes(filename, param_dict)
        updateDataBase(filename, savepath, "isty")

    if not ("pztok" in matricesToLoad):
        df = pandas.read_hdf(savepath + "matricesDataBase.h5", "pztok")
        ind = len(df.index)
        filename = savepath + "mat/pztok_" + version + "_" + str(ind) + ".h5"
        create_file_attributes(filename, param_dict)
        updateDataBase(filename, savepath, "pztok")
        df = pandas.read_hdf(savepath + "matricesDataBase.h5", "pztnok")
        ind = len(df.index)
        filename = savepath + "mat/pztnok_" + version + "_" + str(ind) + ".h5"
        create_file_attributes(filename, param_dict)
        updateDataBase(filename, savepath, "pztnok")
    if not ("imat" in matricesToLoad):
        df = pandas.read_hdf(savepath + "matricesDataBase.h5", "imat")
        ind = len(df.index)
        filename = savepath + "mat/imat_" + version + "_" + str(ind) + ".h5"
        create_file_attributes(filename, param_dict)
        updateDataBase(filename, savepath, "imat")
        df = pandas.read_hdf(savepath + "matricesDataBase.h5", "eigenv")
        ind = len(df.index)
        filename = savepath + "mat/eigenv_" + version + "_" + str(ind) + ".h5"
        create_file_attributes(filename, param_dict)
        updateDataBase(filename, savepath, "eigenv")
        df = pandas.read_hdf(savepath + "matricesDataBase.h5", "U")
        ind = len(df.index)
        filename = savepath + "mat/U_" + version + "_" + str(ind) + ".h5"
        create_file_attributes(filename, param_dict)
        updateDataBase(filename, savepath, "U")


def initDataBase(savepath, param_dict):
    """ Initialize and create the database for all the saved matrices. This database
    will be placed on the top of the savepath and be named matricesDataBase.h5.

    :parameters:
        savepath : (str) : path to the data repertory
        param_dict : (dictionary) : parameters dictionary
    """
    keys = list(param_dict.keys())
    keys.append("path2file")
    keys.append("validity")
    df = pandas.DataFrame(columns=keys)
    store = pandas.HDFStore(savepath + "matricesDataBase.h5")
    store.put("A", df)
    store.put("B", df)
    store.put("istx", df)
    store.put("isty", df)
    store.put("eigenv", df)
    store.put("imat", df)
    store.put("pztok", df)
    store.put("pztnok", df)
    store.put("U", df)

    store.close()
    print("Matrices database created")


def updateDataBase(h5file, savepath, matrix_type):
    """ Update the database adding a new row to the matrix_type database.

    :parameters:
        h5file : (str) : path to the new h5 file to add
        savepath : (str) : path to the data directory
        matrix_type : (str) : type of matrix to store ("A","B","istx","isty"
                                                         "istx","eigenv","imat","U"
                                                         "pztok" or "pztnok")
    """
    if (matrix_type == b"A" or matrix_type == b"B" or matrix_type == b"istx" or
                matrix_type == b"isty" or matrix_type == b"eigenv" or
                matrix_type == b"imat" or matrix_type == b"U" or
                matrix_type == b"pztok" or matrix_type == b"pztnok"):
        f = h5py.File(h5file, "r")
        store = pandas.HDFStore(savepath + "matricesDataBase.h5")
        df = pandas.read_hdf(savepath + "matricesDataBase.h5", matrix_type)
        ind = len(df.index)
        for i in list(f.attrs.keys()):
            df.loc[ind, i] = f.attrs[i]
        df.loc[ind, "path2file"] = h5file
        df.loc[ind, "validity"] = False
        store.put(matrix_type, df)
        store.close()
        f.close()
    else:
        raise ValueError("Wrong matrix_type specified. See documentation")


def save_hdf5(filename, dataname, data):
    """ save_hdf5(filename, dataname, data)
    Create a dataset in an existing hdf5 file filename and store data in it

    :param:
        filename: (str) : full path to the file
        dataname : (str) : name of the data (imat, cmat...)
        data : np.array : data to save
    """
    f = h5py.File(filename, "r+")
    f.create_dataset(dataname, data=data)
    f.close()


def save_h5(filename, dataname, config, data):
    """ save_hdf5(filename, dataname, config, data)
    Create a hdf5 file and store data in it with full header from config parameters
    Usefull to backtrace data origins

    :param:
        filename: (str) : full path to the file
        dataname : (str) : name of the data (imat, cmat...)
        config : (module) : config parameters
        data : np.array : data to save
    """
    p_dict = params_dictionary(config)
    create_file_attributes(filename, p_dict)
    save_hdf5(filename, dataname, data)
    print(filename, "has been written")


def checkMatricesDataBase(savepath, config, param_dict):
    """ Check in the database if the current config have been already run. If so,
    return a dictionary containing the matrices to load and their path. Matrices
    which don't appear in the dictionary will be computed, stored and added
    to the database during the simulation.
    If the database doesn't exist, this function creates it.

    :parameters:
        savepath : (str) : path to the data repertory
        config : (module) : simulation parameters
        param_dict : (dictionary) : parameters dictionary
    :return:
        matricesToLoad : (dictionary) : matrices that will be load and their path
    """

    matricesToLoad = {}
    if (os.path.exists(savepath + "matricesDataBase.h5")):
        checkTurbuParams(savepath, config, param_dict, matricesToLoad)
        checkDmsParams(savepath, config, param_dict, matricesToLoad)
        #        if(matricesToLoad.has_key("pztok")):
        if "pztok" in matricesToLoad:
            checkControlParams(savepath, config, param_dict, matricesToLoad)

    else:
        initDataBase(savepath, param_dict)
    init_hdf5_files(savepath, param_dict, matricesToLoad)
    return matricesToLoad


def checkTurbuParams(savepath, config, pdict, matricesToLoad):
    """ Compare the current turbulence parameters to the database. If similar parameters
    are found, the matricesToLoad dictionary is completed.
    Since all the turbulence matrices are computed together, we only check the parameters
    for the A matrix : if we load A, we load B, istx and isty too.

    :parameters:
        config : (module) : simulation parameters
        matricesToLoad : (dictionary) :  matrices that will be load and their path
    """
    dataBase = pandas.read_hdf(savepath + "matricesDataBase.h5", "A")
    param2test = [
            "r0", "seeds", "L0", "atm.alt", "tel_diam", "cobs", "pupdiam", "zenithangle",
            "target.xpos", "target.ypos", "wfs.xpos", "wfs.ypos"
    ]

    for i in dataBase.index:
        cc = 0
        version = check_output(["git", "rev-parse", "--short", "HEAD"]).decode('utf8')
        if (dataBase.loc[i, "validity"] and (dataBase.loc[i, "revision"] == version)):
            cond = ((dataBase.loc[i, param2test[cc]] == pdict[param2test[cc]]).all())
            while (cond):
                if (cc >= len(param2test)):
                    break
                else:
                    cond = ((dataBase.loc[i, param2test[cc]] == pdict[param2test[cc]]
                             ).all())
                    cc += 1
            # For debug
            #############################
            # if not cond:
            #    cc -= 1
            #    print(param2test[cc]+" has changed from ",dataBase.loc[i,param2test[cc]], " to ",pdict[param2test[cc]])
            ###############################
        else:
            cond = False

        if (cond):
            matricesToLoad["index_turbu"] = i
            matricesToLoad["A"] = dataBase.loc[i, "path2file"]
            dataBase = pandas.read_hdf(savepath + "matricesDataBase.h5", "B")
            matricesToLoad["B"] = dataBase.loc[i, "path2file"]
            dataBase = pandas.read_hdf(savepath + "matricesDataBase.h5", "istx")
            matricesToLoad["istx"] = dataBase.loc[i, "path2file"]
            dataBase = pandas.read_hdf(savepath + "matricesDataBase.h5", "isty")
            matricesToLoad["isty"] = dataBase.loc[i, "path2file"]
            return


def checkControlParams(savepath, config, pdict, matricesToLoad):
    """ Compare the current controller parameters to the database. If similar parameters
    are found, matricesToLoad dictionary is completed.
    Since all the controller matrices are computed together, we only check the parameters
    for the imat matrix : if we load imat, we load eigenv and U too.

    :parameters:
        config : (module) : simulation parameters
        matricesToLoad : (dictionary) :  matrices that will be load and their path
    """
    dataBase = pandas.read_hdf(savepath + "matricesDataBase.h5", "imat")

    param2test = [
            "tel_diam", "t_spiders", "spiders_type", "pupangle", "referr", "std_piston",
            "std_tt", "type_ap", "nbrmissing", "cobs", "pupdiam", "nwfs", "type_wfs",
            "nxsub", "npix", "pixsize", "fracsub", "wfs.xpos", "wfs.ypos", "wfs.Lambda",
            "dms_seen", "fssize", "fstop", "pyr_ampl", "pyr_loc", "pyr_npts",
            "pyr_pup_sep", "pyrtype", "ndms", "type_dm", "dm.alt", "coupling",
            "margin_in", "margin_out", "nact", "nkl", "type_kl", "push4imat",
            "dm.thresh", "unitpervolt", "ncentroiders", "type_centro", "nmax",
            "centro.nwfs", "sizex", "sizey", "centroider.thresh", "type_fct", "weights",
            "width"
    ]

    for i in dataBase.index:
        cc = 0
        version = shesha.__version__
        if (dataBase.loc[i, "validity"] and (dataBase.loc[i, "revision"] == version)):
            cond = ((dataBase.loc[i, param2test[cc]] == pdict[param2test[cc]]).all())
            while (cond):
                if (cc >= len(param2test)):
                    break
                else:
                    cond = ((dataBase.loc[i, param2test[cc]] == pdict[param2test[cc]]
                             ).all())
                    cc += 1
            # For debug
            #############################
            # if not cond:
            #    cc -= 1
            #    print(param2test[cc]+" has changed from ",dataBase.loc[i,param2test[cc]], " to ",pdict[param2test[cc]])
            ###############################
        else:
            cond = False

        if (cond):
            matricesToLoad["index_control"] = i
            matricesToLoad["imat"] = dataBase.loc[i, "path2file"]
            dataBase = pandas.read_hdf(savepath + "matricesDataBase.h5", "eigenv")
            matricesToLoad["eigenv"] = dataBase.loc[i, "path2file"]
            dataBase = pandas.read_hdf(savepath + "matricesDataBase.h5", "U")
            matricesToLoad["U"] = dataBase.loc[i, "path2file"]
            return


def checkDmsParams(savepath, config, pdict, matricesToLoad):
    """ Compare the current controller parameters to the database. If similar parameters
    are found, matricesToLoad dictionary is completed.
    Since all the dms matrices are computed together, we only check the parameters
    for the pztok matrix : if we load pztok, we load pztnok too.

    :parameters:
        config : (module) : simulation parameters
        matricesToLoad : (dictionary) :  matrices that will be load and their path
    """
    dataBase = pandas.read_hdf(savepath + "matricesDataBase.h5", "pztok")

    param2test = [
            "tel_diam", "t_spiders", "spiders_type", "pupangle", "referr", "std_piston",
            "std_tt", "type_ap", "nbrmissing", "cobs", "pupdiam", "nwfs", "type_wfs",
            "nxsub", "npix", "pixsize", "fracsub", "wfs.xpos", "wfs.ypos", "wfs.Lambda",
            "dms_seen", "fssize", "fstop", "pyr_ampl", "pyr_loc", "pyr_npts", "pyrtype",
            "pyr_pup_sep", "ndms", "type_dm", "dm.alt", "coupling", "margin_in",
            "margin_out", "nkl", "nact", "type_kl", "push4imat", "dm.thresh",
            "unitpervolt"
    ]

    for i in dataBase.index:
        cc = 0
        version = check_output(["git", "rev-parse", "--short", "HEAD"]).decode('utf8')
        if (dataBase.loc[i, "validity"] and (dataBase.loc[i, "revision"] == version)):
            cond = ((dataBase.loc[i, param2test[cc]] == pdict[param2test[cc]]).all())
            while (cond):
                if (cc >= len(param2test)):
                    break
                else:
                    cond = ((dataBase.loc[i, param2test[cc]] == pdict[param2test[cc]]
                             ).all())
                    cc += 1
            # For debug
            #############################
            # if not cond:
            #    cc -= 1
            #    print((param2test[cc]+" has changed from ",dataBase.loc[i,param2test[cc]], " to ",pdict[param2test[cc]]))
            ###############################
        else:
            cond = False

        if (cond):
            matricesToLoad["index_dms"] = i
            dataBase = pandas.read_hdf(savepath + "matricesDataBase.h5", "pztnok")
            matricesToLoad["pztnok"] = dataBase.loc[i, "path2file"]
            dataBase = pandas.read_hdf(savepath + "matricesDataBase.h5", "pztok")
            matricesToLoad["pztok"] = dataBase.loc[i, "path2file"]
            return


def validDataBase(savepath, matricesToLoad):
    store = pandas.HDFStore(savepath + "matricesDataBase.h5")
    if not ("A" in matricesToLoad):
        validInStore(store, savepath, "A")
        validInStore(store, savepath, "B")
        validInStore(store, savepath, "istx")
        validInStore(store, savepath, "isty")
    if not ("pztok" in matricesToLoad):
        validInStore(store, savepath, "pztok")
        validInStore(store, savepath, "pztnok")
    if not ("imat" in matricesToLoad):
        validInStore(store, savepath, "imat")
    if not ("eigenv" in matricesToLoad):
        validInStore(store, savepath, "eigenv")
        validInStore(store, savepath, "U")
    store.close()


def validFile(filename):
    f = h5py.File(filename, "r+")
    f.attrs["validity"] = True
    f.close()


def validInStore(store, savepath, matricetype):
    df = pandas.read_hdf(savepath + "matricesDataBase.h5", matricetype)
    ind = len(df.index) - 1
    df.loc[ind, "validity"] = True
    store[matricetype] = df
    validFile(df.loc[ind, "path2file"])


def configFromH5(filename, config):
    #import shesha as ao

    f = h5py.File(filename, "r")

    config.simul_name = str(f.attrs.get("simulname"))
    # Loop
    config.p_loop.set_niter(f.attrs.get("niter"))
    config.p_loop.set_ittime(f.attrs.get("ittime"))

    # geom
    config.p_geom.set_zenithangle(f.attrs.get("zenithangle"))
    config.p_geom.set_pupdiam(f.attrs.get("pupdiam"))

    # Tel
    config.p_tel.set_diam(f.attrs.get("tel_diam"))
    config.p_tel.set_cobs(f.attrs.get("cobs"))
    config.p_tel.set_nbrmissing(f.attrs.get("nbrmissing"))
    config.p_tel.set_t_spiders(f.attrs.get("t_spiders"))
    config.p_tel.set_type_ap(str(f.attrs.get("type_ap")))
    config.p_tel.set_spiders_type(str(f.attrs.get("spiders_type")))
    config.p_tel.set_pupangle(f.attrs.get("pupangle"))
    config.p_tel.set_referr(f.attrs.get("referr"))
    config.p_tel.set_std_piston(f.attrs.get("std_piston"))
    config.p_tel.set_std_tt(f.attrs.get("std_tt"))

    # Atmos
    config.p_atmos.set_r0(f.attrs.get("r0"))
    config.p_atmos.set_nscreens(f.attrs.get("nscreens"))
    config.p_atmos.set_frac(f.attrs.get("frac"))
    config.p_atmos.set_alt(f.attrs.get("atm.alt"))
    config.p_atmos.set_windspeed(f.attrs.get("windspeed"))
    config.p_atmos.set_winddir(f.attrs.get("winddir"))
    config.p_atmos.set_L0(f.attrs.get("L0"))
    config.p_atmos.set_seeds(f.attrs.get("seeds"))

    # Target
    config.p_target.set_nTargets(f.attrs.get("ntargets"))
    config.p_target.set_xpos(f.attrs.get("target.xpos"))
    config.p_target.set_ypos(f.attrs.get("target.ypos"))
    config.p_target.set_Lambda(f.attrs.get("target.Lambda"))
    config.p_target.set_mag(f.attrs.get("target.mag"))
    if (f.attrs.get("target.dms_seen") > -1):
        config.p_target.set_dms_seen(f.attrs.get("target.dms_seen"))

    # WFS
    config.p_wfss = []
    for i in range(f.attrs.get("nwfs")):
        config.p_wfss.append(ao.Param_wfs())
        config.p_wfss[i].set_type(str(f.attrs.get("type_wfs")[i]))
        config.p_wfss[i].set_nxsub(f.attrs.get("nxsub")[i])
        config.p_wfss[i].set_npix(f.attrs.get("npix")[i])
        config.p_wfss[i].set_pixsize(f.attrs.get("pixsize")[i])
        config.p_wfss[i].set_fracsub(f.attrs.get("fracsub")[i])
        config.p_wfss[i].set_xpos(f.attrs.get("wfs.xpos")[i])
        config.p_wfss[i].set_ypos(f.attrs.get("wfs.ypos")[i])
        config.p_wfss[i].set_Lambda(f.attrs.get("wfs.Lambda")[i])
        config.p_wfss[i].set_gsmag(f.attrs.get("gsmag")[i])
        config.p_wfss[i].set_optthroughput(f.attrs.get("optthroughput")[i])
        config.p_wfss[i].set_zerop(f.attrs.get("zerop")[i])
        config.p_wfss[i].set_noise(f.attrs.get("noise")[i])
        config.p_wfss[i].set_atmos_seen(f.attrs.get("atmos_seen")[i])
        config.p_wfss[i].set_fstop(str(f.attrs.get("fstop")[i]))
        config.p_wfss[i].set_pyr_npts(f.attrs.get("pyr_npts")[i])
        config.p_wfss[i].set_pyr_ampl(f.attrs.get("pyr_ampl")[i])
        config.p_wfss[i].set_pyrtype(str(f.attrs.get("pyrtype")[i]))
        config.p_wfss[i].set_pyr_loc(str(f.attrs.get("pyr_loc")[i]))
        config.p_wfss[i].set_fssize(f.attrs.get("fssize")[i])
        if ((f.attrs.get("dms_seen")[i] > -1).all()):
            config.p_wfss[i].set_dms_seen(f.attrs.get("dms_seen")[i])

        # LGS
        config.p_wfss[i].set_gsalt(f.attrs.get("gsalt")[i])
        config.p_wfss[i].set_lltx(f.attrs.get("lltx")[i])
        config.p_wfss[i].set_llty(f.attrs.get("llty")[i])
        config.p_wfss[i].set_laserpower(f.attrs.get("laserpower")[i])
        config.p_wfss[i].set_lgsreturnperwatt(f.attrs.get("lgsreturnperwatt")[i])
        config.p_wfss[i].set_proftype(str(f.attrs.get("proftype")[i]))
        config.p_wfss[i].set_beamsize(f.attrs.get("beamsize")[i])

    # DMs
    config.p_dms = []
    if (f.attrs.get("ndms")):
        for i in range(f.attrs.get("ndms")):
            config.p_dms.append(ao.Param_dm())
            config.p_dms[i].set_type(str(f.attrs.get("type_dm")[i]))
            config.p_dms[i].set_nact(f.attrs.get("nact")[i])
            config.p_dms[i].set_alt(f.attrs.get("dm.alt")[i])
            config.p_dms[i].set_thresh(f.attrs.get("dm.thresh")[i])
            config.p_dms[i].set_coupling(f.attrs.get("coupling")[i])
            config.p_dms[i].set_unitpervolt(f.attrs.get("unitpervolt")[i])
            config.p_dms[i].set_push4imat(f.attrs.get("push4imat")[i])

    # Centroiders
    config.p_centroiders = []
    if (f.attrs.get("ncentroiders")):
        for i in range(f.attrs.get("ncentroiders")):
            config.p_centroiders.append(ao.Param_centroider())
            config.p_centroiders[i].set_nwfs(f.attrs.get("centro.nwfs")[i])
            config.p_centroiders[i].set_type(str(f.attrs.get("type_centro")[i]))
            config.p_centroiders[i].set_type_fct(str(f.attrs.get("type_fct")[i]))
            config.p_centroiders[i].set_nmax(f.attrs.get("nmax")[i])
            config.p_centroiders[i].set_thresh(f.attrs.get("centroider.thresh")[i])
            if (f.attrs.get("weights")[i]):
                config.p_centroiders[i].set_weights(f.attrs.get("weights")[i])
            config.p_centroiders[i].set_width(f.attrs.get("width")[i])
        config.p_rtc.set_centroiders(config.p_centroiders)

    # Controllers
    config.p_controllers = []
    if (f.attrs.get("ncontrollers")):
        for i in range(f.attrs.get("ncontrollers")):
            config.p_controllers.append(ao.Param_controller())
            config.p_controllers[i].set_type(str(f.attrs.get("type_control")[i]))
            config.p_controllers[i].set_nwfs(f.attrs.get("control.nwfs")[i])
            config.p_controllers[i].set_ndm(f.attrs.get("ndm")[i])
            config.p_controllers[i].set_maxcond(f.attrs.get("maxcond")[i])
            config.p_controllers[i].set_delay(f.attrs.get("delay")[i])
            config.p_controllers[i].set_gain(f.attrs.get("gain")[i])
            config.p_controllers[i].set_modopti(f.attrs.get("modopti")[i])
            config.p_controllers[i].set_nrec(f.attrs.get("nrec")[i])
            config.p_controllers[i].set_nmodes(f.attrs.get("nmodes")[i])
            config.p_controllers[i].set_gmin(f.attrs.get("gmin")[i])
            config.p_controllers[i].set_gmax(f.attrs.get("gmax")[i])
            config.p_controllers[i].set_ngain(f.attrs.get("ngain")[i])
            config.p_controllers[i].set_TTcond(f.attrs.get("TTcond")[i])
            config.p_controllers[i].set_cured_ndivs(f.attrs.get("cured_ndivs")[i])
        config.p_rtc.set_controllers(config.p_controllers)

    config.p_rtc.set_nwfs(f.attrs.get("nwfs"))

    print("Parameters have been read from ", filename, "header")


def writeHdf5SingleDataset(filename, data, datasetName="dataset"):
    """Write a hdf5 file containig a single field

    If the file already exists, it will be overwritten
    :parametres:
        filename: (str) : name of the file to write

        data: (np.ndarray) : content of the file

        datasetName: (str) : name of the dataset to write (default="dataset")
    """

    f = h5py.File(filename, "w")
    f.create_dataset(datasetName, data=data)
    f.close()


def readHdf5SingleDataset(filename, datasetName="dataset"):
    """Read a single dataset from an hdf5 file

    :parameters:
        filename: (str) : name of the file to read from

        datasetName: (str) : name of the dataset to read (default="dataset")
    """

    f = h5py.File(filename, "r")
    data = f[datasetName][:]
    f.close()
    return data
