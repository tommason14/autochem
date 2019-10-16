#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ['GaussJob']

"""
File: gaussian.py 
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Interface between Python and GAUSSIAN input files 
"""

from ..core.atom import Atom
from ..core.molecule import Molecule
from ..core.settings import (Settings, read_template)
from ..core.job import Job
from os import (mkdir, chdir, getcwd)
from os.path import (exists, join)
from shutil import (copyfile, move)

__all__ = ['GaussJob']

class GaussJob(Job):
    """Class for creating Gaussian input files and job scripts. 
    
    The names of files created default to the type of calculation: 
    optimisation (opt), single point energy (spec) or hessian matrix 
    calculation for thermochemical data and vibrational frequencies (hess). 
    If a different name is desired, pass a string with the ``filename`` 
    parameter, with no extension. The name will be used for both input and job
    files.
    
        >>> job = GaussJob(using = 'file.xyz', filename = 'benzene')
    
    This command produces 'benzene.job', containing both input data and 
    job scheduler information.

    """

    def __init__(self, using=None, frags_in_subdir=False, settings=None, filename=None):
        super().__init__(using)
        self.filename = filename
        self.defaults = read_template('gaussian.json')
        if settings is not None:
            self.user_settings = settings.as_dict()
            self.merged = self.defaults.merge(settings)
            self.meta = self.merged.meta
            self.input = self.merged.input
        else:
            self.input = self.defaults.input
        self.input = self.input.remove_none_values()
        if '/' in using:
            self.title = using.split('/')[-1][:-4]
        else:
            self.title = using[:-4]
        self.xyz = using

        self.file_basename()
        if frags_in_subdir:
            self.create_inputs_for_fragments()
        self.write_file(self.inp, filetype='job')

    def file_basename(self):
        """
        If no filename is passed when the class is instantiated, the name of the
        file defaults to the run type: a geometry optimisation (opt), single
        point energy calculation (spec), or a hessian matrix calculation for
        vibrational frequencies (freq). This method creates an attribute
        ``base_name``, used in creating the input file."""

        if self.filename is not None:
            self.base_name = self.filename
        else:
            if self.runtype == '':
                self.base_name = 'spec'
            elif 'opt' in self.runtype and 'freq' in self.runtype:
                self.base_name = 'opt-freq'
            elif 'opt' in self.runtype and 'freq' not in self.runtype:
                self.base_name = 'opt'
            else:
                self.base_name = 'freq'

    # def move_to_subdir(self):
    #     """
    #     Makes a subdirectory based on either the filename passed in, 
    #     or the xyz file used to create the input file.
    #     The xyz file is then copied over and the input file is copied
    #     to the new directory.
    #     As a result, this method must be called after making the input file.
    #     """ 
    #     if self.filename is not None:
    #         newdir = join(getcwd(), self.filename)
    #     else:
    #         newdir = join(getcwd(), self.title)
    #     if not exists(newdir):
    #         mkdir(newdir)
    #     copyfile(self.xyz, newdir)
    #     move(f'{self.base_name}.job', newdir)
        

    
    @property
    def job_data(self):
        return self.get_job_template().replace('name', self.base_name)
            
    @property
    def inp(self):
        """
        Creates Gaussian input file of the form:
        SLURM/PBS data
        <blank>
        name
        <blank>
        charge/mult
        xyzdata
        <blank>
        END
        """
        inp = [
            self.job_data,
            self.metadata,
            self.run_info,
            self.title,
            self.coord_info,
            'END'
        ]
        return '\n\n'.join(inp)

    @property
    def metadata(self):
        """
        Include data such as memory and number of cpus in the Gaussian file.
        """
        meta = []
        for k, v in self.meta.items():
            meta.append(f'%{k}={v}')
        return '\n'.join(meta).replace('name', self.base_name)

    @property
    def runtype(self):
        """
        Decides if the job should be an optimisation, single point,
        frequency job, or an optimisation followed by a frequency job.
        Also adds parameters like ts,eigentest to the opt call, for example.
        """
        params = self.input.keys()
        runtype = ''
        if 'opt' in params:
            runtype += gauss_print(self.input, 'opt')
        if 'freq' in params:
            if len(runtype) > 0:
                runtype += ' '
            runtype += 'freq'
        if 'opt' not in params and 'freq' not in params:
            runtype = ''
        return runtype

    @property
    def additional_params(self):
        """
        Add in parameters to the `run_info` that do not involve 
        a basis set, method or run type.
        """
        addn = ''
        ignore=('opt', 'freq', 'method', 'basis', 'meta')
        for arg in self.input:
            if arg not in ignore:
                addn += gauss_print(self.input, arg) + ' '
        return addn

    @property
    def formatted_run(self):
        """
        If a single point calculation is desired, a runtype of '' leaves a
        double space, if '{...} {self.runtype} {...}' is used in `run_info`.
        Instead, use this function and include with no spaces in `run_info`.
        i.e. '{...}{self.formatted_run}{...}'
        """
        return ' ' if self.runtype == '' else f' {self.runtype} '

    @property
    def run_info(self):
        """
        Creating a line like this:
        #P wB97XD/cc-pVDZ opt=(ts,noeigentest,calcfc) freq SCF=tight SCRF=(SMD,solvent=water) INT=(grid=ultrafine)
        from data stored in self.input
        """
        return (f'#P {self.input.method}/{self.input.basis}'
                f'{self.formatted_run}{self.additional_params}')

    def find_charge_and_mult(self):
        """
        Changes charge and multiplicity unless user defines values. 
        In that case, the user-defined charge and multiplicity are used.
        """
        user_assigned_charge = hasattr(self.user_settings, 'input.charge')
        user_assigned_mult = hasattr(self.user_settings, 'input.mult')
        if not user_assigned_charge:
            self.input.charge = self.mol.overall_charge
        if not user_assigned_mult:
            self.input.mult = self.mol.overall_mult        
    
    @property
    def coord_info(self):
        self.find_charge_and_mult()
        info = [f'{self.input.charge} {self.input.mult}']
        info += [f'{atom.symbol:5s} {atom.x:>10.5f} {atom.y:>10.5f} {atom.z:>10.5f}' for atom in self.mol.coords]
        return '\n'.join(info)

def gauss_print(d, value):
    """
    Decides how to print a parameter.
    For example, opt, opt=ts or opt=(ts,eigentest,calcfc). 
    Checks a |Settings| object, `d`, for a key, `value`.
    For example, sett.input.freq=True in the script would
    produce 'freq', sett.input.scf='tight' would produce
    'scf=tight', and sett.input.opt='ts,noeigentest,calcfc' 
    produces 'opt=(ts,noeigentest,calcfc)'. Note: doesn't
    work with lists or dict values, but unlikely that they would
    be passed in as settings values.
    """
    if isinstance(d[f'{value}'], bool) or isinstance(d[f'{value}'], int):
        return value
    elif len(d[f'{value}'].split(',')) > 1 or '=' in d[f'{value}']:
        return f"{value}=({d[f'{value}']})"
    else:
        return f"{value}={d[f'{value}']}"
