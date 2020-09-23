from autochem import Settings, PsiJob
from glob import glob
import os

##############################################################
#  Currently no way of separating fragments without          #
#  asking for counterpoise correction, so use cp=True        #
#  and remove unneccesary files and clean up file structure  #
##############################################################

sett = Settings()
sett.input.memory = '190gb'
sett.input.globals.basis = 'aug-cc-pvdz'

#####################
#  Density fitting  #
#####################

sett.input.globals.df_basis_scf = 'aug-cc-pvdz-jkfit'
sett.input.globals.df_basis_sapt = 'aug-cc-pvdz-ri'

# MP2 defaults need to be removed
sett.input.globals.S_ORTHOGONALIZATION = None
sett.input.globals.freeze_core = None

sett.input.unbound.set = "sapt {\n    print 1\n}"

##############
#  Job info  #
##############

sett.meta.ncpus=48
sett.meta.mem='192gb'
sett.meta.time='24:00:00'
sett.meta.jobfs='200gb'

sett.remove_none_values()

# need the cp=True in order to split the monomers up (counterpoise correction),
# but it places the right input in a subdir, so clean up the file structure
for xyz in glob('*xyz'):
    newdir = f"sapt2plus_{xyz.replace('.xyz', '')}"
    if not os.path.isdir(newdir):
        os.mkdir(newdir)
        os.chdir(newdir)
        PsiJob(f'../{xyz}', cp = True, settings=sett)
        os.system('rm spec.inp spec.job')
        os.system('mv cp-hf/* . && rm -r cp-hf')
        os.system('sed -i "s/^energy.*/energy(\'sapt2+\')/" spec.inp')
        os.chdir('..')
