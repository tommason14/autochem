class Gamess_results:
	"""Class for obtaining results from Gamess simulations. This class requires
	a log file to be read.
	Usage:
		>>> results = Gamess_results(using = 'filename.log')

	Note: assumes that FMO has been used.
	May have to change in future, and write methods for both non-FMO and FMO
	calculations.

	Currently, these methods look for FMO3 data primarily, with an FMO2 fall
back if not found.

	Instances of this class have the following attributes:

	* ``log`` -- filename of the log file of the calculation
	* ``basis`` -- basis set of the calculation, as this class is assumed to be
	* used for ab initio calculations. This attribute may be read from the
	* input file i.e. set as gamess.input.basis = 'CCT')
	* ``coords`` -- coordinates of system, in xyz format
	"""
	def __init__(self, log, basis = None, coords = None):
		self.log = log
		if basis is not None:
			self.basis = basis
		if coords is not None:
			self.coords = coords
			# coords needed for hessian calcs


	def read(self):
		with open(self.log, "r") as f:
			for line in f.readlines():
				yield line

	def get_runtype(self):
		"""Returns type of calculation ran"""
		with open(self.log, "r") as f:
			for line in f.readlines():
				if 'RUNTYP=' in line:
					return line.split()[4].split('=')[1]
		# only want to find the first instance of the run type, might not work
		# with a second successful search

################################
#                              #
#      AB INITIO ENERGIES      #
#                              #
################################

	def get_hf(self):
		"""Returns Hartree-Fock energy"""
		for line in self.read():
			if 'Euncorr HF(3)=' in line:
				return float(line.split()[-1])
			elif 'Euncorr HF(2)=' in line:
				return float(line.split()[-1])

	def get_mp2(self):
		"""Returns MP2 energy"""
		for line in self.read():
			if 'E corr MP2(3)=' in line:
				return float(line.split()[-1])
			elif 'E corr MP2(2)=' in line:
				return float(line.split()[-1])

	def get_srs(self):
		"""Returns SRS-MP2 energy. Looks for SCS as this component is set in
the input file using SRS values"""
		for line in self.read():
			if 'E corr SCS(3)=' in line:
				return float(line.split()[-1])
			elif 'E corr SCS(2)=' in line:
				return float(line.split()[-1])

	def get_mp2_corr(self):
		"""Returns MP2 correlation energy"""
		for line in self.read():
			if 'Edelta MP2(3)=' in line:
				return float(line.split()[-1])
			elif 'Edelta MP2(2)=' in line:
				return float(line.split()[-1])

	def get_srs_corr(self):
		"""Returns SRS-MP2 correlation energy"""
		for line in self.read():
			if 'Edelta SCS(3)=' in line:
				return float(line.split()[-1])
			elif 'Edelta SCS(2)=' in line:
				return float(line.split()[-1])

	def get_e_os(self):
		"""Calculates the opposite spin component to the interaction energy,
for both MP2 and SRS-MP2 results.

		MP2: Corr = E_OS + E_SS, regardless of basis set.
		SRS: Corr = 1.752*E_OS for cc-pVDZ, 1.64*E_OS for cc-pVTZ

		Other basis sets have coefficients for the same spin parameter, these
need to be accounted for also when the time comes.
		Currently only returns the opposite spin energy found from SRS
calculations due to the greater accuracy when compared to MP2.
		"""
		os_coeff = {
			"CCD": 1.752,
			"CCT": 1.64,
		}

		return self.get_srs_corr() / os_coeff[self.basis]

	def get_e_ss(self):
		if self.get_e_os() is not None:
			return self.get_mp2() - self.get_e_os()
		else:
			return None

################################
#                              #
#     VIBRATIONAL ANALYSIS     #
#                              #
################################

	def get_geom(self):
		for i in self.input.
