import chem_assistant as ca
import glob

files = glob.glob('*xyz')

# water properties for each xyz in a directory

for f in sorted(files):
    print(f)
    mol = ca.Molecule(using=f)
    mol.separate()
    for frag in mol.fragments.values():
        angle_done = False
        if frag['name'] == 'water':
            for atom1 in frag['atoms']:
                for atom2 in frag['atoms']:
                    if atom1 != atom2 and atom1.symbol == 'O':
                        length = atom1.distance_to(atom2)
                        print(f'    {atom1.symbol} -- {length:.3f} Å -- {atom2.symbol}')
                        if not angle_done:
                            o = atom1
                            h1, h2 = atom1.connected_atoms
                            angle = o.angle_between(h1.coords,h2.coords)
                            print(f'    {angle:.3f}°')
                            angle_done = True
            print('-'*50) # if more than one water in system, need to discern the output
