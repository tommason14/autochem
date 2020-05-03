from autochem import GaussJob, Settings
import glob

xyz=glob.glob('*xyz')[0]

sett=Settings()
sett.input.opt=True
sett.input.freq=True
sett.input.scrf='smd,solvent=water'

GaussJob(using=xyz, settings=sett)
