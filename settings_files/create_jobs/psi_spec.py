from chem_assistant import PsiJob
import glob

xyz = glob.glob('*xyz')[0]

PsiJob(using=xyz, frags_in_subdir=True)
