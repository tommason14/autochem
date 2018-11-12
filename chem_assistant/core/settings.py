from os.path import (join, dirname)
import json

__all__ = ['Settings', 'read_template', 'dict_to_settings']

class Settings(dict):
    """Provides a means of updating settings for input and job files. Inherits from a python dictionary, and allows for nesting of dictionaries as values."""
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        for k, v in self.items():
            if isinstance(v, dict): # nesting
                self[k] = Settings(v)

    def as_dict(self):
        """Return a copy as a regular python dictionary"""
        d = {}
        for k, v in self.items():
            if isinstance(v, Settings):
                d[k] = v.as_dict() #required for multi-level
            else:
                d[k] = v
        return d

    def __iter__(self):
        """Iterate through in alphabetical order"""
        return iter(sorted(self.keys()))
    
    def _str(self, indent = 0):
        """Print dict with 2 space indenting"""
        ret = ''
        for name in self:
            value = self[name]
            ret += ' '*indent + str(name) + ':    '
            if isinstance(value, Settings):
                ret += '\n' + value._str(indent+len(str(name))+1)
            else:
                ret += str(value) + '\n'
        return ret

    def __str__(self):
        return self._str()

    def __setitem__(self, key, value):
        """Like regular ``__setitem__``, but needs adjusting for the nested possibility"""
        if isinstance(value, dict):
            value = Settings(value)
        dict.__setitem__(self, key, value)
    
    def __getattr__(self, key):
        """If key is not a magic method, redirect it to ``__getitem__``."""
        if (key.startswith('__') and key.endswith('__')):
            return dict.__getattr__(self, key)
        return self[key]


    def __setattr__(self, key, value):
        """If key is not a magic method, redirect it to ``__setitem__``."""
        if key.startswith('__') and key.endswith('__'):
            dict.__setattr__(self, key, value)
        self[key] = value


    def __delattr__(self, key):
        """If key is not a magic method, redirect it to ``__delitem__``."""
        if key.startswith('__') and key.endswith('__'):
            dict.__delattr__(self, key)
        del self[key]
    
    def __missing__(self, name):
        """When requested key is not present, add it with an empty |Settings| instance as a value.
        This method is essential for automatic insertions in deeper levels. Without it things like::
            >>> s = Settings()
            >>> s.a.b.c = 12
        will not work.
        """
        self[name] = Settings()
        return self[name]

    def copy(self):
        ret = Settings()
        for name in self:
            if isinstance(self[name], Settings):
                ret[name] = self[name].copy()
            else:
                ret[name] = self[name]
        return ret

    def merge(self, other):
        """
        Return new instance of |Settings| that is a copy of this instance
        updated with *other*.
        """
        ret = self.copy()
        ret.update(other)
        return ret

    def update(self, other):
        """Updates the current |Settings| object with the object passed in"""
        for name in other:
            if isinstance(other[name], Settings):
                if name not in self or not isinstance(self[name], Settings):
                    self[name] = other[name].copy()
                else:
                    self[name].update(other[name])
            else:
                self[name] = other[name]

    __repr__ = __str__


def read_template(template):
    """Obtains default parameters for input files of different packages, and returns them as a |Settings| object. Currently GAMESS and PSI4 are supported"""
    path = join(dirname(__file__), '..', 'templates')
    file = join(path, template)
    with open(file, "r") as f:
        tmp = json.load(f)
    return dict_to_settings(tmp)


def dict_to_settings(d):
    """Transform a python dictionary into a |Settings| object. This function works recursively, checking if nested dictionaries are present"""
    s = Settings()
    for k, v in d.items():
        if isinstance(v, dict):
            s[k] = dict_to_settings(v)
        else:
            s[k] = v
    return s
