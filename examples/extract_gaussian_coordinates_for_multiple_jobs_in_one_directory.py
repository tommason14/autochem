from autochem import GaussianResults
from glob import glob


class GaussResults(GaussianResults):
    def __init__(self, log):
        super().__init__(log)

    @property
    def title(self):
        """
        Redefine self.title so that the xyzs from example.log are written to spec/example-equil.xyz or
        rerun/example-rerun.xyz
        """
        return self.log.split('.')[0]


for log in glob('*log'):
    res = GaussResults(log)
    res.get_equil_coords()
