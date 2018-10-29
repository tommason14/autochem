from periodic_table import PeriodicTable as PT

class Molecule:
	"""A class representing a molecule,  comprised of atoms

	An instance of this class has the following attributes:
	
	* ``atoms`` -- a list of |Atom| objects
	* ``bonds`` -- a list of bonds between the ``atoms`` in a the molecule
	...

	"""


	def __init__(self, filename = None, inputformat = None, geometry = 1):
		self.atoms = []
		self.bonds = []

		if filename is not None:
		self.read(filename, inputformat, geometry)


    def readxyz(self, f, frame):
  
        def newatom(line):
            lst = line.split()
            shift = 1 if (len(lst) > 4 and lst[0] == str(i)) else 0
            num = lst[0+shift]
            if isinstance(num, str):
                num = PT.get_atomic_number(num)
            self.add_atom(Atom(atnum=num, coords=(lst[1+shift],lst[2+shift],lst[3+shift])))

        def newlatticevec(line):
            lst = line.split()
            self.lattice.append((float(lst[1]),float(lst[2]),float(lst[3])))

        fr = frame
        begin, first, nohead = True, True, False
        for line in f:
            if first:
                if line.strip() == '' : continue
                first = False
                try:
                    n = int(line.strip())
                    fr -= 1
                except ValueError:
                    nohead = True
                    newatom(line)
            elif nohead:
                if line.strip() == '' : break
                if 'VEC' in line.upper():
                    newlatticevec(line)
                else:
                    newatom(line)
            elif fr != 0:
                try:
                    n = int(line.strip())
                    fr -= 1
                except ValueError:
                    continue
            else:
                if begin:
                    begin = False
                    i = 1
                    if line:
                        self.properties['comment'] = line.rstrip()
                else:
                    if i <= n:
                        newatom(line)
                        i += 1
                    elif 'VEC' in line.upper():
                       newlatticevec(line)
                    else:
                        break
        if not nohead and fr > 0:
            raise FileError('readxyz: There are only %i frames in %s' % (frame - fr, f.name))


    def writexyz(self, f):
        f.write(str(len(self)) + '\n')
        if 'comment' in self.properties:
            comment = self.properties['comment']
            if isinstance(comment, list):
                comment = comment[0]
            f.write(comment)
        f.write('\n')
        for at in self.atoms:
            f.write(str(at) + '\n')
        for i,vec in enumerate(self.lattice):
            f.write('VEC'+str(i+1) + '%14.6f %14.6f %14.6f\n'%tuple(vec))


    def readmol(self, f, frame):
        if frame != 1:
            raise FileError('readmol: .mol files do not support multiple geometries')

        comment = []
        for i in range(4):
            line = f.readline().rstrip()
            if line:
                spl = line.split()
                if spl[-1] == 'V2000':
                    if len(line) == 39:
                        natom = int(line[0:3])
                        nbond = int(line[3:6])
                    else:
                        natom = int(spl[0])
                        nbond = int(spl[1])
                    for j in range(natom):
                        atomline = f.readline().rstrip()
                        if len(atomline) == 69:
                            crd = (float(atomline[:10]),float(atomline[10:20]),float(atomline[20:30]))
                            symb = atomline[31:34].strip()
                        else:
                            tmp = atomline.split()
                            crd = tuple(map(float, tmp[0:3]))
                            symb = tmp[3]
                        try:
                            num = PT.get_atomic_number(symb)
                        except PTError:
                            num = 0
                        self.add_atom(Atom(atnum=num, coords=crd))
                    for j in range(nbond):
                        bondline = f.readline().rstrip()
                        if len(bondline) == 21:
                            at1 = int(bondline[0:3])
                            at2 = int(bondline[3:6])
                            ordr = int(bondline[6:9])
                        else:
                            tmp = bondline.split()
                            at1 = int(tmp[0])
                            at2 = int(tmp[1])
                            ordr = int(tmp[2])
                        if ordr == 4:
                            ordr = Bond.AR
                        self.add_bond(Bond(atom1=self[at1], atom2=self[at2], order=ordr))
                    break
                elif spl[-1] == 'V3000':
                    raise FileError('readmol: Molfile V3000 not supported. Please convert')
                else:
                    comment.append(line)
        if comment:
            self.properties['comment'] = comment



    def writemol(self, f):
        commentblock = ['\n']*3
        if 'comment' in self.properties:
            comment = self.properties['comment']
            if isinstance(comment, str):
                commentblock[0] = comment + '\n'
            elif isinstance(comment, list):
                comment = comment[0:3]
                while len(comment) < 3:
                    comment.append('')
                commentblock = [a+b for a,b in zip(comment,commentblock)]
        f.writelines(commentblock)

        self.set_atoms_id()

        f.write('%3i %2i  0  0  0  0  0  0  0  0999 V2000\n' % (len(self.atoms),len(self.bonds)))
        for at in self.atoms:
            f.write('%10.4f %9.4f %9.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0\n' % (at.x,at.y,at.z,at.symbol))
        for bo in self.bonds:
            order = bo.order
            if order == Bond.AR:
                order = 4
            f.write('%3i %2i %2i  0  0  0  0\n' % (bo.atom1.id,bo.atom2.id,order))
        self.unset_atoms_id()
        f.write('M  END\n')



    def readmol2(self, f, frame):
        if frame != 1:
            raise MoleculeError('readmol: .mol2 files do not support multiple geometries')

        bondorders = {'1':1, '2':2, '3':3, 'am':1, 'ar':Bond.AR, 'du':0, 'un':1, 'nc':0}
        mode = ('', 0)
        for i, line in enumerate(f):
            line = line.rstrip()
            if not line:
                continue
            elif line[0] == '#':
                continue
            elif line[0] == '@':
                line = line.partition('>')[2]
                if not line:
                    raise FileError('readmol2: Error in %s line %i: invalid @ record' % (f.name, str(i+1)))
                mode = (line, i)

            elif mode[0] == 'MOLECULE':
                pos = i - mode[1]
                if pos == 1:
                    self.properties['name'] = line
                elif pos == 3:
                    self.properties['type'] = line
                elif pos == 4:
                    self.properties['charge_type'] = line
                elif pos == 5:
                    self.properties['flags'] = line
                elif pos == 6:
                    self.properties['comment'] = line

            elif mode[0] == 'ATOM':
                spl = line.split()
                if len(spl) < 6:
                    raise FileError('readmol2: Error in %s line %i: not enough values in line' % (f.name, str(i+1)))
                symb = spl[5].partition('.')[0]
                try:
                    num = PT.get_atomic_number(symb)
                except PTError:
                    num = 0
                crd = tuple(map(float, spl[2:5]))
                newatom = Atom(atnum=num, coords=crd, name=spl[1], type=spl[5])
                if len(spl) > 6:
                    newatom.properties['subst_id'] = spl[6]
                if len(spl) > 7:
                    newatom.properties['subst_name'] = spl[7]
                if len(spl) > 8:
                    newatom.properties['charge'] = float(spl[8])
                if len(spl) > 9:
                    newatom.properties['flags'] = spl[9]
                self.add_atom(newatom)

            elif mode[0] == 'BOND':
                spl = line.split()
                if len(spl) < 4:
                    raise FileError('readmol2: Error in %s line %i: not enough values in line' % (f.name, str(i+1)))
                try:
                    atom1 = self.atoms[int(spl[1])-1]
                    atom2 = self.atoms[int(spl[2])-1]
                except IndexError:
                    raise FileError('readmol2: Error in %s line %i: wrong atom ID' % (f.name, str(i+1)))
                newbond = Bond(atom1, atom2, order=bondorders[spl[3]])
                if len(spl) > 4:
                    for flag in spl[4].split('|'):
                        newbond.properties[flag] = True
                self.add_bond(newbond)



    def writemol2(self, f):

        def write_prop(name, obj, separator, space=0, replacement=None):
            form_str = '%-' + str(space) + 's'
            if name in obj.properties:
                f.write(form_str % str(obj.properties[name]))
            elif replacement is not None:
                f.write(form_str % str(replacement))
            f.write(separator)

        f.write('@<TRIPOS>MOLECULE\n')
        write_prop('name', self, '\n')
        f.write('%i %i\n' % (len(self.atoms),len(self.bonds)))
        write_prop('type', self, '\n')
        write_prop('charge_type', self, '\n')
        write_prop('flags', self, '\n')
        write_prop('comment', self, '\n')

        f.write('\n@<TRIPOS>ATOM\n')
        for i,at in enumerate(self.atoms):
            f.write('%5i ' % (i+1))
            write_prop('name', at, ' ', 5, at.symbol+str(i+1))
            f.write('%10.4f %10.4f %10.4f ' % at.coords)
            write_prop('type', at, ' ', 5, at.symbol)
            write_prop('subst_id', at, ' ', 5)
            write_prop('subst_name', at, ' ', 7)
            write_prop('charge', at, ' ', 6)
            write_prop('flags', at, '\n')
            at.id = i+1

        f.write('\n@<TRIPOS>BOND\n')
        for i,bo in enumerate(self.bonds):
            f.write('%5i %5i %5i %4s' % (i+1, bo.atom1.id, bo.atom2.id, 'ar' if bo.is_aromatic() else bo.order))
            write_prop('flags', bo, '\n')

        self.unset_atoms_id()



    def readpdb(self, f, frame):
        pdb = PDBHandler(f)
        models = pdb.get_models()
        if frame > len(models):
            raise FileError('readpdb: There are only %i frames in %s' % (len(models), f.name))

        symbol_columns = [70,6,7,8]
        for i in models[frame-1]:
            if i.name in ['ATOM  ','HETATM']:
                x = float(i.value[0][24:32])
                y = float(i.value[0][32:40])
                z = float(i.value[0][40:48])
                for n in symbol_columns:
                    symbol = i.value[0][n:n+2].strip()
                    try:
                        atnum = PT.get_atomic_number(symbol)
                        break
                    except PTError:
                        if n == symbol_columns[-1]:
                            raise FileError('readpdb: Unable to deduce the atomic symbol in the following line:\n%s'%(i.name+i.value[0]))
                self.add_atom(Atom(atnum=atnum,coords=(x,y,z)))

        return pdb



    def writepdb(self, f):
        pdb = PDBHandler()
        pdb.add_record(PDBRecord('HEADER'))
        model = []
        for i,at in enumerate(self.atoms):
            s = 'ATOM  %5i                   %8.3f%8.3f%8.3f                      %2s  ' % (i+1,at.x,at.y,at.z,at.symbol.upper())
            model.append(PDBRecord(s))
        pdb.add_model(model)
        pdb.add_record(pdb.calc_master())
        pdb.add_record(PDBRecord('END'))
        pdb.write(f)


    def read(self, filename, inputformat=None, geometry=1):
        """Read molecular coordinates from a file.
        *filename* should be a string with a path to a file. If *inputformat* is not ``None``, it should be one of supported formats (keys occurring in the class attribute ``_readformat``). Otherwise, the format is deduced from the file extension. For files without an extension the `xyz` format is used.
        If the chosen format allows multiple geometries in a single file, *geometry* can be used to pick one of them.
        """
        if inputformat is None:
            fsplit = filename.rsplit('.',1)
            if len(fsplit) == 2:
                inputformat = fsplit[1]
            else:
                inputformat = 'xyz'
        if inputformat in self.__class__._readformat:
            with open(filename, 'r') as f:
                ret = self._readformat[inputformat](self, f, geometry)
            return ret
        else:
            raise MoleculeError('read: Unsupported file format')



    def write(self, filename, outputformat=None):
        """Write molecular coordinates to a file.
        *filename* should be a string with a path to a file. If *outputformat* is not ``None``, it should be one of supported formats (keys occurring in the class attribute ``_writeformat``). Otherwise, the format is deduced from the file extension. For files without an extension the `xyz` format is used.
        """
        if outputformat is None:
            fsplit = filename.rsplit('.',1)
            if len(fsplit) == 2:
                outputformat = fsplit[1]
            else:
                outputformat = 'xyz'
        if outputformat in self.__class__._writeformat:
            with open(filename, 'w') as f:
                self._writeformat[outputformat](self, f)
        else:
            raise MoleculeError('write: Unsupported file format')

    _readformat = {'xyz':readxyz, 'mol':readmol, 'mol2':readmol2, 'pdb':readpdb}
    _writeformat = {'xyz':writexyz, 'mol':writemol, 'mol2':writemol2, 'pdb': writepdb}



    def as_dict(self):
        """Store all information about the molecule in a dictionary.
        The returned dictionary is, in principle, identical to ``self.__dict__`` of the current instance, apart from the fact that all |Atom| and |Bond| instances in ``atoms`` and ``bonds`` lists are replaced with dictionaries storing corresponing information.
        This method is a counterpart of :meth:`from_dict`.
        """
        mol_dict = copy.copy(self.__dict__)
        atom_indices = {id(a): i for i, a in enumerate(mol_dict['atoms'])}
        bond_indices = {id(b): i for i, b in enumerate(mol_dict['bonds'])}
        atom_dicts = [copy.copy(a.__dict__) for a in mol_dict['atoms']]
        bond_dicts = [copy.copy(b.__dict__) for b in mol_dict['bonds']]
        for a_dict in atom_dicts:
            a_dict['bonds'] = [bond_indices[id(b)] for b in a_dict['bonds']]
            del(a_dict['mol'])
        for b_dict in bond_dicts:
            b_dict['atom1'] = atom_indices[id(b_dict['atom1'])]
            b_dict['atom2'] = atom_indices[id(b_dict['atom2'])]
            del(b_dict['mol'])
        mol_dict['atoms'] = atom_dicts
        mol_dict['bonds'] = bond_dicts
        return mol_dict



    @classmethod
    def from_dict(cls, dictionary):
        """Generate a new |Molecule| instance based on the information stored in a *dictionary*.
        This method is a counterpart of :meth:`as_dict`.
        """
        mol = cls()
        mol.__dict__ = copy.copy(dictionary)
        atom_dicts = mol.atoms
        bond_dicts = mol.bonds
        mol.atoms=[]
        mol.bonds=[]
        for a_dict in atom_dicts:
            a = Atom()
            a.__dict__ = a_dict
            a.mol = mol
            a.bonds=[]
            mol.add_atom(a)
        for b_dict in bond_dicts:
            b = Bond(None, None)
            b_dict['atom1'] = mol.atoms[b_dict['atom1']]
            b_dict['atom2'] = mol.atoms[b_dict['atom2']]
            b.__dict__ = b_dict
            b.mol = mol
            mol.add_bond(b)
        return mol


	def add_atom(self, atom, adjacent = None):
  		"""Add an *atom* to the molecule"""
		self.atoms.append(atom)
		atom.mol = self
		#add bonds
		if adjacent is not None:
  			for adj in adjacent:
  				  	if isinstance(adj, tuple): #two attachments, in a chain 
  						self.add_bond(atom, adj[0], adj[1])
					else:
  						self.add_bond(atom, adj) #one atom attached i.e. O-h
	
	def add_bond(self, arg1, arg2=None, order=1):
  		"""Add a new bond between two atoms in the molecule.
		  
		Can be called in two ways:
		* create a |Bond| object then pass to method:
			>>> b = Bond(mol[1], mol[4], order=2) # create a double bond between 1st and 4th atom
			>>> mol.add_bond(b)
		* pass two atoms and a bond order directly into method:
			>>> mol.add_bond(mol[1], mol[4], order=2)
		
		Both atoms must be contained in the same molecule, otherwise an exception is raised.
		"""
		if isinstance(arg1, Atom) and isinstance(arg2, Atom):
  			newbond = Bond(arg1, arg2, order = order)
		elif isinstance(arg1, Bond):
  			newbond = arg1
		else:
  			raise MoleculeError('Invalid argument passed when trying to add a new bond')
		
		# check that atoms are in the molecule, then add the bond to the same molecule
		if newbond.atom1.mol == self and newbond.atom2.mol == self:
  			newbond.mol = self
			# add bonds to bonds array of the molecule, and to the bonds array for each atom
			self.bonds.append(newbond)
			newbond.atom1.bonds.append(newbond)
			newbond.atom2.bonds.append(newbond)
		else:
  			MoleculeError('add_bond: Bonded atoms must belong to the same molecule')


	def delete_atom(self, atom):
  		"""Delete an *atom* from the molecule.

		*atom* has to be an |Atom| instance; all bonds attached to this atom are removed also.

		Examples::
            #delete all hydrogens
            mol = Molecule('protein.pdb')
            hydrogens = [atom for atom in mol if atom.atnum == 1]
            for i in hydrogens: mol.delete_atom(i)
        ::
            #delete first two atoms
            mol = Molecule('geom.xyz')
            mol.delete_atom(mol[1])
            mol.delete_atom(mol[1]) #since the second atom of original molecule is now the first

		"""
		if atom.mol != self:
  			raise MoleculeError('delete_atom: Atom needs to belong to the molecule')
		try:
  			self.atoms.remove(atom)
		except:
			raise MoleculeError('delete_atom: invalid argument passed as an atom')
		atom.mol = None
		for b in reversed(atom.bonds):
  			self.delete_bond(b)

	def delete_bond(self, arg1, arg2 = None):
  		"""Delete a bond from the molecule.
        Just like :meth:`add_bond`, this method accepts either a single argument that is a |Bond| instance, or two arguments being instances of |Atom|. In both cases objects used as arguments have to belong to the molecule.
        """
        if isinstance(arg1, Atom) and isinstance(arg2, Atom):
            delbond = self.find_bond(arg1, arg2)
        elif isinstance(arg1, Bond):
            delbond = arg1
        else:
            raise MoleculeError('delete_bond: invalid arguments passed')
        if delbond in self.bonds:
            delbond.mol = None
            self.bonds.remove(delbond)
            delbond.atom1.bonds.remove(delbond)
            delbond.atom2.bonds.remove(delbond)

	def find_bond(self, atom1, atom2):
  		"""Find the bond between *atom1* and *atom2*."""
		if atom1.mol != self or atom2.mol != self:
  			raise MoleculeError('find_bond: Atoms have to belong to the molecule')
		for b in atom1.bonds:
  			if atom2 is b.other_end(atom1):
  				return b
		return None

	def neighbours(self, atom):
  		"""Finds the atoms that are connected to the atom passed."""
		
		if atom.mol != self:
  			raise MoleculeError('neigbours: Atom has to belong to the molecule')
		return [b.other_end(atom) for b in atom.bonds]

	def separate():
  		"""Separate the molecule into connected components.
        Returnsa list of new |Molecule| objects. Each element of this list is identical to one connected component of the base molecule- function uses a copy of the original molecule.

        Example::
            >>> mol = Molecule('xyz_dimers/NH3-H2O.xyz')
            >>> mol.guess_bonds()
            >>> print(mol)
              Atoms:
                1         N     -1.395591     -0.021564      0.000037
                2         H     -1.629811      0.961096     -0.106224
                3         H     -1.862767     -0.512544     -0.755974
                4         H     -1.833547     -0.330770      0.862307
                5         O      1.568501      0.105892      0.000005
                6         H      0.606736     -0.033962     -0.000628
                7         H      1.940519     -0.780005      0.000222
              Bonds:
               (5)--1.0--(7)
               (5)--1.0--(6)
               (1)--1.0--(3)
               (1)--1.0--(4)
               (1)--1.0--(2)
            >>> x = mol.separate()
            >>> for i in x: print(i)
              Atoms:
                1         N     -1.395591     -0.021564      0.000037
                2         H     -1.629811      0.961096     -0.106224
                3         H     -1.862767     -0.512544     -0.755974
                4         H     -1.833547     -0.330770      0.862307
              Bonds:
               (1)--1.0--(3)
               (1)--1.0--(4)
               (1)--1.0--(2)
              Atoms:
                1         O      1.568501      0.105892      0.000005
                2         H      0.606736     -0.033962     -0.000628
                3         H      1.940519     -0.780005      0.000222
              Bonds:
               (1)--1.0--(3)
               (1)--1.0--(2)
        """

		frags = []
        clone = self.copy()
        for at in clone:
            at._visited = False

        def dfs(v, mol):
            v._visited = True
            v.mol = mol
            for e in v.bonds:
                e.mol = mol
                u = e.other_end(v)
                if not u._visited:
                    dfs(u, mol)

        for src in clone.atoms:
            if not src._visited:
                m = Molecule()
                dfs(src, m)
                frags.append(m)

        for at in clone.atoms:
            del at._visited
            at.mol.atoms.append(at)
        for b in clone.bonds:
            b.mol.bonds.append(b)

        return frags

    def guess_bonds(self):
        """Try to guess bonds in the molecule based on types and positions of atoms.
        All previously existing bonds are removed. New bonds are generated based on interatomic distances and information about maximal number of bonds for each atom type (``connectors`` property, taken from |PeriodicTable|).
        The problem of finding molecular bonds for a given set of atoms in space does not have a general solution, especially considering the fact the chemical bond in itself is not a precisely defined concept. For every method, no matter how sophisticated, there will always be corner cases for which the method produces disputable results. Moreover, depending on the context (area of application) the desired solution for a particular geometry may vary. Please do not treat this method as an oracle always providing a proper solution. The algorithm used here gives very good results for geometries that are not very far from the optimal geometry, especially consisting of lighter atoms. All kinds of organic molecules, including aromatic ones, usually work very well. Problematic results can emerge for transition metal complexes, transition states, incomplete molecules etc.
        The algorithm used scales as *n log n* where *n* is the number of atoms.
        .. warning::
            This method works reliably only for geometries representing complete molecules. If some atoms are missing (for example, a protein without hydrogens) the resulting set of bonds would usually contain more bonds or bonds with higher order than expected.
        """

        class HeapElement(object):
            def __init__(self, order, ratio, atom1, atom2):
                eff_ord = order
                if order == 1.5: #effective order for aromatic bonds
                    eff_ord = 1.15
                elif order == 1 and {atom1.symbol, atom2.symbol} == {'C', 'N'}:
                    eff_ord = 1.11 #effective order for single C-N bond
                value = (eff_ord + 0.9) * ratio
                self.data = (value, order, ratio)
                self.atoms = (atom1, atom2)
            def unpack(self):
                val, o, r = self.data
                at1, at2 = self.atoms
                return val, o, r, at1, at2
            def __lt__(self, other): return self.data < other.data
            def __le__(self, other): return self.data <= other.data
            def __eq__(self, other): return self.data == other.data
            def __ne__(self, other): return self.data != other.data
            def __gt__(self, other): return self.data > other.data
            def __ge__(self, other): return self.data >= other.data

        self.delete_all_bonds()

        dmax = 1.28

        cubesize = dmax*2.1*max([at.radius for at in self.atoms])

        cubes = {}
        for i,at in enumerate(self.atoms):
            at._id = i+1
            at.free = at.connectors
            at.cube = tuple(map(lambda x: int(math.floor(x/cubesize)), at.coords))
            if at.cube in cubes:
                cubes[at.cube].append(at)
            else:
                cubes[at.cube] = [at]

        neighbors = {}
        for cube in cubes:
            neighbors[cube] = []
            for i in range(cube[0]-1, cube[0]+2):
                for j in range(cube[1]-1, cube[1]+2):
                    for k in range(cube[2]-1, cube[2]+2):
                        if (i,j,k) in cubes:
                            neighbors[cube] += cubes[(i,j,k)]

        heap = []
        for at1 in self.atoms:
            if at1.free > 0:
                for at2 in neighbors[at1.cube]:
                    if (at2.free > 0) and (at1._id < at2._id):
                        ratio = at1.distance_to(at2)/(at1.radius+at2.radius)
                        if (ratio < dmax):
                            heap.append(HeapElement(0, ratio, at1, at2))
                            #I hate to do this, but I guess there's no other way :/ [MH]
                            if (at1.atnum == 16 and at2.atnum == 8):
                                at1.free = 6
                            elif (at2.atnum == 16 and at1.atnum == 8):
                                at2.free = 6
                            elif (at1.atnum == 7):
                                at1.free += 1
                            elif (at2.atnum == 7):
                                at2.free += 1
        heapq.heapify(heap)

        for at in self.atoms:
            if at.atnum == 7:
                if at.free > 6:
                    at.free = 4
                else:
                    at.free = 3

        while heap:
            val, o, r, at1, at2 = heapq.heappop(heap).unpack()
            step = 1 if o in [0,2] else 0.5
            if at1.free >= step and at2.free >= step:
                o += step
                at1.free -= step
                at2.free -= step
                if o < 3:
                    heapq.heappush(heap, HeapElement(o,r,at1,at2))
                else:
                    self.add_bond(at1,at2,o)
            elif o > 0:
                if o == 1.5:
                    o = Bond.AR
                self.add_bond(at1,at2,o)

        def dfs(atom, par):
            atom.arom += 1000
            for b in atom.bonds:
                oe = b.other_end(atom)
                if b.is_aromatic() and oe.arom < 1000:
                    if oe.arom > 2:
                        return False
                    if par and oe.arom == 1:
                        b.order = 2
                        return True
                    if dfs(oe, 1-par):
                        b.order = 1 + par
                        return True

        for at in self.atoms:
            at.arom = len(list(filter(Bond.is_aromatic, at.bonds)))

        for at in self.atoms:
            if at.arom == 1:
                dfs(at, 1)

        for at in self.atoms:
            del at.cube,at.free,at._id,at.arom

#===========================================================================
#==== Geometry operations ==================================================
#===========================================================================



    def translate(self, vector):
        """Move the molecule in space by *vector*, expressed in angstroms.
        *vector* should be an iterable container of length 3 (usually tuple, list or numpy array). *unit* describes unit of values stored in *vector*.
        """
        for at in self.atoms:
            at.translate(vector)


    def rotate_lattice(self, matrix):
        """Rotate **only** lattice vectors of the molecule with given rotation *matrix*.
        *matrix* should be a container with 9 numerical values. It can be a list (tuple, numpy array etc.) listing matrix elements row-wise, either flat (``[1,2,3,4,5,6,7,8,9]``) or in two-level fashion (``[[1,2,3],[4,5,6],[7,8,9]]``).
        .. note::
            This method does not check if *matrix* is a proper rotation matrix.
        """
        self.lattice = [tuple(np.dot(matrix,i)) for i in self.lattice]


    def rotate(self, matrix, lattice=False):
        """Rotate the molecule with given rotation *matrix*. If *lattice* is ``True``, rotate lattice vectors too.
        *matrix* should be a container with 9 numerical values. It can be a list (tuple, numpy array etc.) listing matrix elements row-wise, either flat (``[1,2,3,4,5,6,7,8,9]``) or in two-level fashion (``[[1,2,3],[4,5,6],[7,8,9]]``).
        .. note::
            This method does not check if *matrix* is a proper rotation matrix.
        """
        for at in self.atoms:
            at.rotate(matrix)
        if lattice:
            self.rotate_lattice(matrix)


    def align_lattice(self, convention='x', zero=1e-10):
        """Rotate the molecule in such a way that lattice vectors are aligned with the coordinate system.
        This method is meant to be used with periodic systems only. Using it on a |Molecule| instance with an empty ``lattice`` attribute has no effect.
        Possible values of the *convention* argument are:
        *   ``x`` (default) -- first lattice vector aligned with X axis. Second vector (if present) aligned with XY plane.
        *   ``z`` (convention used by `ReaxFF <https://www.scm.com/product/reaxff>`_) -- second lattice vector (if present) aligned with YZ plane. Third vector (if present) aligned with Z axis.
        *zero* argument can be used to specify the numerical tolerance for zero (used to determine if some vector is already aligned with a particular axis or plane).
        The returned boolean value indicates if any rotation happened.
        """
        dim = len(self.lattice)

        if dim == 0:
            log('NOTE: align_lattice called on a Molecule without any lattice', 5)
            return False

        rotated = False
        if convention == 'x':
            if abs(self.lattice[0][1]) > zero or abs(self.lattice[0][2]) > zero:
                mat = rotation_matrix(self.lattice[0], [1.0, 0.0, 0.0])
                self.rotate(mat, lattice=True)
                rotated = True

            if dim >= 2 and abs(self.lattice[1][2]) > zero:
                mat = rotation_matrix([0.0, self.lattice[1][1], self.lattice[1][2]], [0.0, 1.0, 0.0])
                self.rotate(mat, lattice=True)
                rotated = True

        elif convention == 'z':
            if dim == 3 and (abs(self.lattice[2][0]) > zero or abs(self.lattice[2][1]) > zero):
                mat = rotation_matrix(self.lattice[2], [0.0, 0.0, 1.0])
                self.rotate(mat, lattice=True)
                rotated = True

            if dim >= 2 and abs(self.lattice[1][0]) > zero:
                mat = rotation_matrix([self.lattice[1][0], self.lattice[1][1], 0.0], [0.0, 1.0, 0.0])
                self.rotate(mat, lattice=True)
                rotated = True

        else:
            raise MoleculeError("align_lattice: unknown convention: {}. Possible values are 'x' or 'z'".format(convention))
        return rotated


    def rotate_bond(self, bond, atom, angle):
        """Rotate given *bond* by an *angle* expressed in radians.
        *bond* should be chosen in such a way, that it divides the molecule into two parts (using a bond being part of a ring results in an error). *atom* has to belong to *bond* and is used to pick which "half" of the molecule is rotated. A positive angle denotes counterclockwise rotation (when looking along the bond, from the stationary part of the molecule).
        """
        if atom not in bond:
            raise MoleculeError('rotate_bond: atom has to belong to the bond')

        atoms_to_rotate = {atom}

        def dfs(v):
            for e in v.bonds:
                if e is not bond:
                    u = e.other_end(v)
                    if u not in atoms_to_rotate:
                        atoms_to_rotate.add(u)
                        dfs(u)

        dfs(atom)

        if len(atoms_to_rotate) == len(self):
            raise MoleculeError('rotate_bond: chosen bond does not divide molecule')

        other_end = bond.other_end(atom)
        v = np.array(other_end.vector_to(atom))
        v /= np.linalg.norm(v)

        W = np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

        a1 = math.sin(angle)
        a2 = 2 * math.pow(math.sin(0.5 * angle), 2)

        rotmat = np.identity(3) + a1 * W + a2 * np.dot(W,W)

        trans = np.array(other_end.vector_to((0,0,0)))
        for at in atoms_to_rotate:
            at.translate(trans)
            at.rotate(rotmat)
            at.translate(-trans)


    def closest_atom(self, point):
        """Return the atom of the molecule that is the closest one to some *point* in space.
        *point* should be an iterable container of length 3 (for example: tuple, |Atom|, list, numpy array). *unit* describes unit of values stored in *point*.
        """
        dist = float('inf')
        for at in self.atoms:
            newdist = at.distance_to(point)
            if newdist < dist:
                dist = newdist
                ret = at
        return ret


    def distance_to_point(self, point):
        """Calculate the distance between the molecule and some *point* in space (distance between *point* and :meth:`closest_atom`).
        *point* should be an iterable container of length 3 (for example: tuple, |Atom|, list, numpy array). *unit* describes unit of values stored in *point*. Returned value is expressed in *result_unit*.
        """
        at = self.closest_atom(point)
        return at.distance_to(point)


    def distance_to_mol(self, other, result_unit='angstrom', return_atoms=False):
        """Calculate the distance between the molecule and some *other* molecule.
        The distance is measured as the smallest distance between any atom of this molecule and any atom of *other* molecule. Returned distance is expressed in *result_unit*.
        If *return_atoms* is ``False``, only a single number is returned.  If *return_atoms* is ``True``, the method returns a tuple ``(distance, atom1, atom2)`` where ``atom1`` and ``atom2`` are atoms fulfilling the minimal distance, with atom1 belonging to this molecule and atom2 to *other*.
        """
        dist = float('inf')
        for at1 in self.atoms:
            for at2 in other.atoms:
                newdist = (at1.x-at2.x)**2 + (at1.y-at2.y)**2 + (at1.z-at2.z)**2
                if newdist < dist:
                    dist = newdist
                    atom1 = at1
                    atom2 = at2
        res = math.sqrt(dist)
        if return_atoms:
            return res, atom1, atom2
        return res


    def wrap(self, length, angle=2*math.pi):
        """wrap(self, length, angle=2*pi)
        Transform the molecule wrapping its x-axis around z-axis. This method is useful for building nanotubes or molecular wedding rings.
        Atomic coordinates are transformed in the following way:
        *   z coordinates remain untouched
        *   x axis gets wrapped around the circle centered in the origin of new coordinate system. Each segment of x axis of length *length* ends up as an arc of a circle subtended by an angle *angle*. The radius of this circle is R = *length*/*angle*.
        *   part of the plane between the x axis and the line y=R is transformed into the interior of the circle, with line y=R being squashed into a single point - the center of the circle.
        *   part of the plane above line y=R is dropped
        *   part of the plane below x axis is transformed into outside of the circle
        *   transformation is done in such a way that distances along y axis are preserved
        Before:
        .. image:: ../_static/wrap.*
        After:
        .. image:: ../_static/wrap2.*
        """

        xs = [atom.x for atom in self.atoms]
        if max(xs)-min(xs) > length:
            raise MoleculeError('wrap: x-extension of the molecule is larger than length')

        if angle < 0 or angle > 2*math.pi:
            raise MoleculeError('wrap: angle must be between 0 and 2*pi')

        R = length / angle

        def map_ring(x,y):
            return ((R-y) * math.cos(x/R), (R-y) * math.sin(x/R))

        for at in self.atoms:
            at.x, at.y = map_ring(at.x, at.y)


    def get_center_of_mass(self):
        """Return the center of mass of the molecule (as a tuple). Returned coordinates are expressed in *unit*."""
        center = [0.0,0.0,0.0]
        total_mass = 0.0
        for at in self.atoms:
            total_mass += at.mass
            for i in range(3):
                center[i] += at.mass*at.coords[i]
        for i in range(3):
            center[i] = center[i]/total_mass
        return tuple(center)


    def get_mass(self):
        """Return the mass of the molecule, expressed in atomic units."""
        return sum([at.mass for at in self.atoms])


    def get_formula(self, as_dict=False):
        """Calculate the molecular formula of the molecule.
        Here molecular formula is a dictionary with keys being atomic symbols. The value for each key is the number of atoms of that type. If *as_dict* is ``True``, that dictionary is returned. Otherwise, it is converted into a string::
            >>> mol = Molecule('Ubiquitin.xyz')
            >>> print(m.get_formula(True))
            {'N': 105, 'C': 378, 'O': 118, 'S': 1, 'H': 629}
            >>> print(m.get_formula(False))
            C378H629N105O118S1
        """
        ret = {}
        for atom in self:
            if atom.symbol not in ret:
                ret[atom.symbol] = 0
            ret[atom.symbol] +=1
        if as_dict:
            return ret
        s = ''
        for key in sorted(ret):
            s += '{}{}'.format(key,ret[key])
        return s


#===========================================================================
#==== Magic methods ========================================================
#===========================================================================



    def __len__(self):
        """The length of the molecule is the number of atoms."""
        return len(self.atoms)

    def __str__(self):
        """Return a string representation of the molecule.
        Information about atoms is printed in ``xyz`` format fashion -- each atom in a separate, enumerated line. Then, if the molecule contains any bonds, they are printed. Each bond is printed in a separate line, with information about both atoms and bond order. Example:
        .. code-block:: none
                  Atoms:
                    1         N       0.00000       0.00000       0.38321
                    2         H       0.94218       0.00000      -0.01737
                    3         H      -0.47109       0.81595      -0.01737
                    4         H      -0.47109      -0.81595      -0.01737
                  Bonds:
                    (1)----1----(2)
                    (1)----1----(3)
                    (1)----1----(4)
        """
        s = '  Atoms: \n'
        for i,atom in enumerate(self.atoms):
            s += ('%5i'%(i+1)) + str(atom) + '\n'
        if len(self.bonds) > 0:
            for j,atom in enumerate(self.atoms):
                atom._tmpid = j+1
            s += '  Bonds: \n'
            for bond in self.bonds:
                s += '   (%d)--%1.1f--(%d)\n'%(bond.atom1._tmpid, bond.order, bond.atom2._tmpid)
            for atom in self.atoms:
                del atom._tmpid
        if self.lattice:
            s += '  Lattice:\n'
            for vec in self.lattice:
                s += '    {:16.10f} {:16.10f} {:16.10f}\n'.format(*vec)
        return s


    def __iter__(self):
        """Iterate over atoms."""
        return iter(self.atoms)


    def __getitem__(self, key):
        """The bracket notation can be used to access atoms or bonds directly.
        If *key* is a single int (``mymol[i]``), return i-th atom of the molecule. If *key* is a pair of ints (``mymol[(i,j)]``), return the bond between i-th and j-th atom (``None`` if such a bond does not exist). Negative integers can be used to access atoms enumerated in the reversed order.
        This notation is read only: things like ``mymol[3] = Atom(...)`` are forbidden.
        Numbering of atoms within a molecule starts with 1.
        """
        if isinstance(key, int):
            if key == 0:
                raise MoleculeError('Numbering of atoms starts with 1')
            if key < 0:
                return self.atoms[key]
            return self.atoms[key-1]
        if isinstance(key, tuple) and len(key) == 2:
            return self.find_bond(self[key[0]], self[key[1]])
        raise MoleculeError('Molecule: invalid argument {} inside []'.format(key))


    def __add__(self, other):
        """Create a new molecule that is a sum of this molecule and some *other* molecule::
            newmol = mol1 + mol2
        The new molecule has atoms, bonds and all other elements distinct from both components. The ``properties`` of ``newmol`` are a copy of the ``properties`` of ``mol1`` :meth:`soft_updated<scm.plams.core.settings.Settings.soft_update>` with the ``properties`` of ``mol2``.
        """
        m = self.copy()
        m += other
        return m


    def __iadd__(self, other):
        """Add some *other* molecule to this one::
            protein += water
        All atoms and bonds present in *other* are copied and copies are added to this molecule. The ``properties`` of this molecule are :meth:`soft_updated<scm.plams.core.settings.Settings.soft_update>` with the  ``properties`` of the *other* molecules.
        """
        othercopy = other.copy()
        self.atoms += othercopy.atoms
        self.bonds += othercopy.bonds
        for atom in self.atoms:
            atom.mol = self
        for bond in self.bonds:
            bond.mol = self
        self.properties.soft_update(othercopy.properties)
        return self


    def __copy__(self):
        return self.copy()
