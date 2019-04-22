import numpy as np


class PressureVessel:
	def __init__(self, R, t, p_i, p_o):
		self.R = R
		self.t = t
		self.p_i = p_i
		self.p_o = p_o
	
	@property
	def area(self):
		return np.pi*self.R**2
	
	@property
	def sig_r(self):
		return -(self.p_i + self.p_o)/2
	
	@property
	def sig_a(self):
		return (self.R/self.t)*(self.p_i - self.p_o)
	
	@property
	def sig_z(self):
		return -(self.p_o*self.R)/(2*self.t)
	
	def _get_diffs(self):
		return np.array((self.sig_r - self.sig_a,
		                 self.sig_r - self.sig_z,
		                 self.sig_a - self.sig_z))
	
	def get_tresca_stress(self):
		return max(abs(self._get_diffs()))
	
	def get_von_mises_stress(self):
		return np.sqrt(0.5*(self._get_diffs()**2).sum())
	
	def get_max_stress(self):
		return max(self.get_tresca_stress(), self.get_von_mises_stress())


if __name__ == "__main__":
	# test
	pv1 = PressureVessel(2, 0.2, 2, 10)
	print(pv1.get_max_stress())
		