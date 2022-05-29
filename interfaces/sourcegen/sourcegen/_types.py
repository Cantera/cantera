from abc import ABCMeta, abstractmethod
from dataclasses import dataclass, fields
import os


def _unpack(obj):
    return (getattr(obj, field.name) for field in fields(obj))   


@dataclass
class Param:
    """ Represents a function parameter """

    p_type: str
    name: str


    @staticmethod
    def parse(c_param: str):
        parts = c_param.strip().rsplit(' ', 1)
        if len(parts) == 2:
            return Param(*parts)


    def __iter__(self):
        return _unpack(self)


@dataclass
class Func:
    """ Represents a function parsed from a C header file """

    ret_type: str
    name: str
    params: list[Param]


    @staticmethod
    def parse(c_func: str):
        lparen = c_func.index('(')
        rparen = c_func.index(')')
        front = c_func[0:lparen].split()
        
        params = (Param.parse(p) for p in c_func[lparen+1:rparen].split(','))
        params = [p for p in params if p]
        
        ret_type = front[-2]
        name = front[-1]
        return Func(ret_type, name, params)


    def __iter__(self):
        return _unpack(self)


class SourceGeneratorBase(metaclass=ABCMeta):
    """ Specifies the interface of a language-specific SourceGenerator  """

    @abstractmethod
    def __init__(self, out_dir: str, config: dict):
        pass

    @abstractmethod
    def generate_source(self, incl_file: os.DirEntry, funcs: list[Func]):
        pass

    @abstractmethod
    def finalize(self):
        pass
