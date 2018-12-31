from chem_assistant import GamessJob, Settings

s = Settings()
s.input.contrl.runtyp = 'energy'
s.input.basis.gbasis = 'cct'
GamessJob(fmo = True, frags_in_subdir = True, using = 'equil.xyz', settings = s)