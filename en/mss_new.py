import matplotlib.pyplot as plt
import numpy as np
from math import log, tanh,exp
import H_M as hm
import pandas as pd

R = 8.31441

class Fluid:

	def __init__(self, rc, pc, Tc, pt, Tt, M, df, name):
		self.rc = rc
		self.pc = pc
		self.Tc = Tc
		self.pt = pt
		self.Tt = Tt
		self.M = M
		self.u_nist = df.Viscosity
		self.PT = df[["Pressure", "Temperature", "Density"]]
		self.name = name

#		
		#
	def __str__(self):
		return f"{self.name} :\nrc = {self.rc}\npc = {self.pc}\nTc = {self.Tc}\npt = {self.pt}\nTt = {self.Tt}\nM = {self.M}\nPT :\n{self.PT}"

def Rr(ro,rco):
	return ro/rco
def alfa(Rr,Mo,Mm):
	ao = 1 + 0.007378*(Rr**1.847)*(Mo**0.5173)
	am = 1 + 0.007378*(Rr**1.847)*(Mm**0.5173)
	return ao,am
	
def main():
	colnames1 = ["Temperature", "Pressure", "Density", "Volume", "Internal_Energy", "Enthalpy", "Entropy",\
				"Cv", "Cp", "Sound_Spd", "Joule-Thomson", "Viscosity", "Therm_Cond", "Surf_Tension" ,"Phase"]
	colnames2 = ["Temperature", "Pressure", "Density", "Volume", "Internal_Energy", "Enthalpy", "Entropy",\
				"Cv", "Cp", "Sound_Spd", "Joule-Thomson", "Viscosity", "Therm_Cond", "Phase"]
	df_p_l = pd.read_csv("liquid_propan.csv", names=colnames1, header=None)
	#df_p_v = pd.read_csv("vapor_propan.csv", names=colnames2, header=None)

	df_m_l = pd.read_csv("liquid_metan.csv", names=colnames1, header=None)

	Tc = 190.4
	pc = 46
	rc = 163.5 
	metan = Fluid(0.1628, 45.992, 190.564, 0.11696, 90.67, 16.043, df_m_l, "Metan")
	propan = Fluid(0.220, 42.4766, 369.825, 1.7*10**(-9), 85.48, 44.097, df_p_l, "Propane")
	#To = t*(metan.Tc-metan.Tt)/(propan.Tc-propan.Tt)-(metan.Tc*propan.Tt-metan.Tt*propan.Tc)/(propan.Tc-propan.Tt)
	#print(0.35 *(metan.Tc-metan.Tt)/(propan.Tc-propan.Tt))
	#print(90*(metan.Tc-metan.Tt)/(propan.Tc-propan.Tt)-(metan.Tc*propan.Tt-metan.Tt*propan.Tc)/(propan.Tc-propan.Tt))
	#print(300*(metan.Tc-metan.Tt)/(propan.Tc-propan.Tt)-(metan.Tc*propan.Tt-metan.Tt*propan.Tc)/(propan.Tc-propan.Tt))
	#print(metan)
	#print(propan)
	rcm = 0.220
	Pcm = 42.4766
	#Tcm=400
	Tcm = 369.825
	Mm = 44.097
	
	Tco = 190.564
	Pco = 45.992
	rco = 0.1628
	Mo = 16.043
	
	Tto = 90.6941
	Pto = 0.11696
	
	Ttm = 85.48
	Ptm = 1.7*10**(-9)

	_,_,_,a,b,c,f0,_,_,_,_,A,GV=hm.import_data()

	print(metan)
	print(propan)

	#kp=(Pco-Pto)/(Pcm-Ptm)
	#bp=(Pco*Ptm-Pto*Pcm)/(Pcm-Ptm)
	EEE=((propan.Tc/metan.Tc)**(-1/6.))*((propan.pc/metan.pc)**(2/3.))*((propan.M/metan.M)**(1/2.))
	V=np.zeros(len(metan.PT.Density))
	for j,i in enumerate(metan.PT.iterrows()):
		#print(i[1])
		V[j] = EEE * hm.Viscosity_t(i[1].Density * 1.e-5,metan.rc,GV,i[1].Temperature,a,b,c,f0,A)
		
	plt.plot(propan.PT.Temperature, V, "r", propan.PT.Temperature, df_p_l.Viscosity, "b")
	plt.show()

if __name__ == "__main__":
	main()
