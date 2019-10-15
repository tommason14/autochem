from chem_assistant import Settings, GaussJob
import glob

xyz=glob.glob('.xyz')[0]

sett=Settings()
sett.input.method='wB97XD'
sett.input.basis='aug-cc-pVDZ' 
sett.input.opt=True
sett.input.td='NSTATES=10,root=5'
sett.input.scrf='smd,solvent=water'

GaussJob(using=xyz, settings=sett)
