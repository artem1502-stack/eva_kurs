import pandas as pd
import numpy as np
from scipy.optimize import fsolve

R = 8.31441


def BWR(ro, p, T):
	
	N = np.array(N)
	a = np.array(a)
	

N = [0, -1.8439486666 * 1.e-2, 1.0510162064, -16.057820303, 848.44027562, -42738.409106, 7.6565285254 * 1.e-4,\
	-0.48360724197, 85.195473835, -16607.434721, -3.7521074532 * 1.e-5, 0.028616309259,\
	5.7974531455*1.e-6, -2.8685295973, 1.1906973942*1.e-4, -8.5315715699*1.e-3, 3.8365063841, 2.4986828379*1.e-5,\
	5.7974531455*1.e-6, -7.1648329297*1.e-3, 1.2577853784*1.e-4, 22240.102466, -1480051.2328,\
	50.498054887, 1642837.5992, 0.21325387196, 37.791273422, -1.1857016815*1.e-5, -31.630780767, -4.1006782941 *1.e-6,\
	1.4870043284*1.e-3, 3.1512261532*1.e-9, -2.1670774745*1.e-6, 12.4000551079*1.e-5]
N = np.array(N)
y = 0.0096

def bvr(ro, p, T):
	#p, T = args[0], args[1]

	a = [0, R*T, N[1]*T + N[2]*T**(0.5) + N[3] + N[4]/T + N[5]/T**2, N[6]*T + N[7] + N[8]/T + N[9]/T**2, N[10]*T + N[11] + N[12]*T,\
	N[13], N[14]/T + N[15]/T**2, N[16]/T, N[17]/T + N[18]/T**2, N[19]/T**2, N[20]/T**2 + N[21]/T**3, N[22]/T**2 + N[23]/T**4,\
	N[24]/T**2 + N[25]/T**3, N[26]/T**2 + N[27]/T**4, N[28]/T**2 + N[29]/T**3, N[30]/T**2 + N[31]/T**3 + N[32]/T**4]
	a = np.array(a)

	res = -p
	for i in range(1, 10):
		res += a[i] * ro**i
	for i in range(10, 16):
		res += a[i] * ro ** (2 * i - 17) * np.exp(-y * ro**2)
	return res

def find_ro(p, T, r0):
	ro = fsolve(bvr, r0, args = (p, T))
	return ro

class Fluid:

	def __init__(self, rc, pc, Tc, pt, Tt, M, df, name):
		self.rc = rc
		self.pc = pc
		self.Tc = Tc
		self.pt = pt
		self.Tt = Tt
		self.M = M
		self.u_nist = df.Viscosity
		self.PT = df[["Pressure", "Temperature"]]
		self.name = name

	def __str__(self):
		return f"{self.name} :\nrc = {self.rc}\npc = {self.pc}\nTc = {self.Tc}\npt = {self.pt}\nTt = {self.Tt}\nM = {self.M}\nPT :\n{self.PT}"

def main():
	colnames1 = ["Temperature", "Pressure", "Density", "Volume", "Internal_Energy", "Enthalpy", "Entropy",\
				"Cv", "Cp", "Sound_Spd", "Joule-Thomson", "Viscosity", "Therm_Cond", "Surf_Tension" ,"Phase"]
	colnames2 = ["Temperature", "Pressure", "Density", "Volume", "Internal_Energy", "Enthalpy", "Entropy",\
				"Cv", "Cp", "Sound_Spd", "Joule-Thomson", "Viscosity", "Therm_Cond", "Phase"] 
	df_liquid = pd.read_csv("liq_methan.csv", names=colnames1, header=None)
	df_vapor = pd.read_csv("vapor.csv", names=colnames2, header=None)

	Tc = 190.4
	pc = 46
	rc = 163.5 
	metan = Fluid(190.4, 46, 163.5, 0.11696, 90.67, 16.04, df_liquid, "Metan")
	propan = Fluid(5000, 42.512, 369.89, 1.7 * 1.e-9, 85.48, 44.097, df_liquid, "Propane")
	
	print(df_liquid.iloc[300])
	s = find_ro(6.4118, 140, 400)
	print(s)

if __name__ == "__main__":
	main()
