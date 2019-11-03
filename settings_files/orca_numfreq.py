from chem_assistant import Settings

sett = Settings()
sett.input.method='M062X'
sett.input.density_fitting=None
sett.input.NumFreq = True
sett.input.basis='6-311+G(d,p)'
sett.input.solvent.model='cpcm'
sett.input.solvent.SMD='true'
sett.input.solvent.Solvent='"Water"'
