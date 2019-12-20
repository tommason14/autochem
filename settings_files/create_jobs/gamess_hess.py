from chem_assistant import Settings, GamessJob
import glob

xyz=glob.glob('*xyz')[0]

sett=Settings()
sett.input.contrl.runtyp='hessian'
# sett.input.force.method='seminum' # big systems

GamessJob(using=xyz, fmo=True, frags_in_subdir=True, settings=sett)
