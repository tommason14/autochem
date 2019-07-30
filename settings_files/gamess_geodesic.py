from chem_assistant import Settings

sett = Settings()
sett.input.basis.gbasis='cct'
sett.input.contrl.runtyp='energy'
sett.input.elpot.iepot=1
sett.input.elpot.where='pdc'
sett.input.pdc.ptsel='geodesic'

# rm defaults
sett.input.contrl.mplevl=None
sett.input.mp2=None
