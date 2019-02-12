import importlib


def smart_import(mod, cls, verbose=False, silent=False):
    try:
        if verbose:
            print("trying from " + mod + " import " + cls)
        # my_module = __import__(mod, globals(), locals(), [cls], 0)
        my_module = importlib.import_module(mod)
        return getattr(my_module, cls)

    except ImportError as err:
        if not silent:
            import warnings
            warnings.warn(
                    "Error importing %s, it will be simulated due to: %s" %
                    (cls, err.msg), Warning)

        class tmp_cls:

            def __init__(self, *args, **kwargs):
                raise RuntimeError("Can not initilize the simulation with fake objects")

        return tmp_cls

    except AttributeError as err:
        if not silent:
            import warnings
            warnings.warn(
                    "Error importing %s, it will be simulated due to: %s" %
                    (cls, err.args), Warning)

        class tmp_cls:

            def __init__(self, *args, **kwargs):
                raise RuntimeError("Can not initilize the simulation with fake objects")

        return tmp_cls


Dms = smart_import("sutraWrap", "Dms")
Rtc_FFF = smart_import("sutraWrap", "Rtc_FFF")
Rtc_FHF = smart_import("sutraWrap", "Rtc_FHF")
Rtc_UFF = smart_import("sutraWrap", "Rtc_UFF")
Rtc_UHF = smart_import("sutraWrap", "Rtc_UHF")
Rtc_FFU = smart_import("sutraWrap", "Rtc_FFU")
Rtc_FHU = smart_import("sutraWrap", "Rtc_FHU")
Rtc_UFU = smart_import("sutraWrap", "Rtc_UFU")
Rtc_UHU = smart_import("sutraWrap", "Rtc_UHU")
Rtc_brahma = smart_import("sutraWrap", "Rtc_brahma", silent=True)
Rtc_cacao = smart_import("sutraWrap", "Rtc_cacao", silent=True)
Sensors = smart_import("sutraWrap", "Sensors")
Atmos = smart_import("sutraWrap", "Atmos")
Telescope = smart_import("sutraWrap", "Telescope")
Target = smart_import("sutraWrap", "Target")
Target_brahma = smart_import("sutraWrap", "Target_brahma", silent=True)
Gamora = smart_import("sutraWrap", "Gamora")
Groot = smart_import("sutraWrap", "Groot")

carmaWrap_context = smart_import("carmaWrap", "context")
