
import uonidtoolbox as unit

class Struct:
    def __init__(self, existing: dict):
        if existing:
            self.__dict__ = existing
            for k in self.__dict__.keys():
                if isinstance(self.__dict__[k], dict):
                    self.__dict__[k] = Struct(self.__dict__[k])
                #endif
            #endfor
        #endif
    #enddef
#endclass

