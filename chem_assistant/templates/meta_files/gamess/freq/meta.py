from chem_assistant import GamessJob, Settings

s = Settings()
s.input.contrl.runtyp = 'hessian'
# s.input.force.method = 'seminum'
GamessJob(fmo = True, frags_in_subdir = True, using = 'equil.xyz', settings = s)
