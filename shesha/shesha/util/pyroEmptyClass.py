"""
This file is here only to avoid pyro server crash when a coronograph class is not in the pramam file. 
(Pyro does not manage non callable in server...)

"""

class PyroEmptyClass():
    """ 
    Just a useless empty class for pyro server. 
    """
    def __init__(self):
        pass