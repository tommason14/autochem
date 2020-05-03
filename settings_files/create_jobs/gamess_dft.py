from autochem import Settings, GamessJob
import glob

xyz=glob.glob('*xyz')[0]

sett=Settings()
sett.input.mp2=None
sett.input.contrl.mplevl=None
sett.input.contrl.dfttyp='m062x'
sett.input.dft.method='grid'

GamessJob(using=xyz, fmo=True, frags_in_subdir=True, settings=sett)
