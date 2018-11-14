import re

def read(filename):
    with open(filename, "r") as f:
        for line in f.readlines(): #one at a time
            yield line

def get_runtype():
    for line in read('anion_3.ok'):
        #regex for energy('mp2') or optimize('scf', dertype='hess') (any number of additional args)
        if re.search("[A-z]*\('[A-z0-9]*'(.?\s*[A-z]*='[A-z]*')*\)", line):
            if re.search("[A-z]*\('[A-z0-9]*'\)", line): #energy('mp2')
                return line.split('(')[0]
            else: #optimize('scf', dertype='hess'......)
                return line.split('(')[0] #add to this later, using the collect additional data
