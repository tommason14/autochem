from chem_assistant import Settings

sett=Settings()
sett.input.run='wB97X-D3 aug-cc-pVTZ RIJCOSX CPCM'
sett.input.meta.tddft="""\
  maxdim 5
  nroots 10"""
sett.input.meta.cpcm="""\
  SMD true
  SMDSolvent \"DMSO\" """
