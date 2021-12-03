from autochem import Settings

sett = Settings()
sett.input.method = "MP2"
sett.input.basis = "cc-pVTZ"
sett.input.pop = "chelpg"
sett.meta.ncpus = 24
sett.meta.mem = "96gb"
sett.meta.nodemem = "96gb"
sett.meta.time = "2:00:00"
