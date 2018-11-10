def read(filename):
    with open(filename, "r") as f:
        for line in f.readlines(): #one at a time
            yield line

def completed():
    found = False
    for line in read('c4mim-ac-p4.log'):
        if 'EXECUTION OF GAMESS TERMINATED NORMALLY' in line:
            found = True
    return found

def get_error():
    pass

if not completed():
    print('Nope')
