from chem_assistant import Settings, GamessJob
import glob

xyz = glob.glob('*xyz')[0]

sett = Settings()
sett.input.basis.gbasis = 'cct'
sett.input.contrl.runtyp = 'energy'

GamessJob(using=xyz, fmo=True, frags_in_subdir=True, settings=sett)
