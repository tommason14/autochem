from autochem import Settings, GamessJob

sett = Settings()
sett.input.basis.gbasis = "cct"
sett.input.contrl.runtyp = "energy"
sett.input.elpot.iepot = 1
sett.input.elpot.where = "pdc"
sett.input.pdc.ptsel = "geodesic"

# rm defaults
sett.input.contrl.mplevl = None
sett.input.mp2 = None

g = GamessJob("methane.xyz", settings=sett)

with open("methane_geodesic.inp") as f:
    # strip newlines from end of line, keep indentation at start
    reference = [l.rstrip() for l in f.readlines()]

assert g.input.split("\n") == reference
