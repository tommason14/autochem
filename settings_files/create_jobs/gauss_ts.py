from chem_assistant import GaussJob, Settings
import glob

xyz=glob.glob('*xyz')[0]

sett=Settings()
sett.input.opt='ts,noeigentest,calcfc'
sett.input.freq=True
sett.input.scrf='smd,solvent=water'

GaussJob(using=xyz, settings=sett)
