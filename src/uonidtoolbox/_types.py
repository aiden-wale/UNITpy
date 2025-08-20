
import uonidtoolbox as unit
import numpy as np


class struct:
    def __init__(self, existing: dict = {}):

        if not isinstance(existing, (dict, struct)):
            raise Exception("unit.struct can only be created from 'dict' or another unit.struct")
        #endif

        if not existing:
            return
        #endif

        if isinstance(existing, dict):
            self.__dict__ = existing.copy()
            for k in self.__dict__.keys():
                if isinstance(self.__dict__[k], dict):
                    self.__dict__[k] = struct(self.__dict__[k])
                elif isinstance(self.__dict__[k], (list, np.ndarray)):
                    for i in range(0, len(self.__dict__[k])):
                        if isinstance(self.__dict__[k][i], dict):
                            self.__dict__[k][i] = struct(self.__dict__[k][i])
                        #endif
                    #endfor
                #endif
            #endfor
        elif isinstance(existing, struct):
            tmp = struct(existing.__dict__)
            self.__dict__ = tmp.__dict__
        #endif
    #enddef

    def __iter__(self):
        for item in self.__dict__:
            yield item
        #endfor
    #enddef

    def __len__(self):
        return len(self.__dict__)
    #enddef

    def __getitem__(self, key):
        return self.__dict__[key]
    #enddef

    def __setitem__(self, key, value):
        self.__dict__[key] = value
    #enddef

    def __delitem__(self, key):
        del self.__dict__[key]
    #enddef

    def __repr__(self):
        return repr(self.__dict__)
    #enddef



    def asdict(self):
        d = self.__dict__.copy()
        for k in d.keys():
            if isinstance(d[k], struct):
                d[k] = d[k].asdict()
            elif isinstance(d[k], (list, np.ndarray)):
                for i in range(0, len(d[k])):
                    if isinstance(d[k][i], unit.struct):
                        d[k][i] = d[k][i].asdict()
                    #endif
                #endfor
            #endif
        #endfor
        return d
    #enddef

    def keys(self):
        return self.__dict__.keys()
    #enddef

    def items(self):
        return self.__dict__.items()
    #enddef

    def copy(self):
        return struct(self.__dict__.copy())
    #enddef

#endclass


