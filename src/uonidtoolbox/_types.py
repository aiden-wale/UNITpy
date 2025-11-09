
from numpy import ndarray
from copy import deepcopy


class _struct:
    def __init__(self, given_data = {}):

        if not isinstance(given_data, (dict, _struct)):
            raise Exception("uonidtoolbox._struct can only be created from 'dict' or another uonidtoolbox._struct")

        if not given_data:
            return

        _data = deepcopy(given_data)

        if isinstance(_data, dict):
            self.__dict__ = _data.copy()
            for k in self.__dict__.keys():
                if isinstance(self.__dict__[k], dict):
                    self.__dict__[k] = _struct(self.__dict__[k])
                elif isinstance(self.__dict__[k], (list, ndarray)):
                    for i in range(0, len(self.__dict__[k])):
                        if isinstance(self.__dict__[k][i], dict):
                            self.__dict__[k][i] = _struct(self.__dict__[k][i])
                        #endif
                    #endfor
                #endif
            #endfor
        elif isinstance(_data, _struct):
            tmp = _struct(_data.__dict__)
            self.__dict__ = tmp.__dict__
        #endif
    #enddef __init__

    def __len__(self): return len(self.__dict__)
    def __getitem__(self, key): return self.__dict__[key]
    def __setitem__(self, key, value): self.__dict__[key] = value
    def __delitem__(self, key): del self.__dict__[key]
    def __repr__(self): return repr(self.__dict__)

    def __iter__(self):
        for item in self.__dict__:
            yield item
        #endfor

    def keys(self): return self.__dict__.keys()
    def items(self): return self.__dict__.items()

    def copy(self, deep=False):
        if deep:
            return deepcopy(self)
        else:
            return self.__dict__.copy()


    def asdict(self):
        d = self.__dict__.copy()
        for k in d.keys():
            if isinstance(d[k], _struct):
                d[k] = d[k].asdict()
            elif isinstance(d[k], (list, ndarray)):
                for i in range(0, len(d[k])):
                    if isinstance(d[k][i], _struct):
                        d[k][i] = d[k][i].asdict()
                    #endif
                #endfor
            #endif
        #endfor
        return d
#endclass _struct


# class _Z_structure(_struct):
#     def __init__(self, _data):
#         super().__init__(_data)

        


