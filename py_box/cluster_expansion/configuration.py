class Configuration(object):
    """
    Holds the configuration.
        Attributes:
            sigma - list of int
                Vector that holds the state of the site. Usually either +1 or -1.
            E_DFT - float
                Energy of the configuration obtained by DFT
            E_CE - float
                Enregy of the configuration obtained by cluster expansion
    """

    def __init__(self, name = None, sigma = None, E_DFT = None, E_CE = None):
        self.name = name
        self.sigma = sigma
        self.E_DFT = E_DFT
        self.E_CE = E_CE

default_dict = {'name': str,
                'sigma': list,
                'E_DFT': float,
                'E_CE': float,
                'E_DFT_Raw': float}
