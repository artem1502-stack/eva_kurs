import matplotlib.pyplot as plt
import numpy as np
from math import log, tanh,exp
def import_data(d=None):
	np.set_printoptions(precision=7)
	GV_t=np.zeros(9)
	A=np.zeros(7)

	A[0]=-10.2391
	A[1]=174.2282
	A[2]=17.4605
	A[3]=-2847.6328
	A[4]=0.1337
	A[5]=142.0724
	A[6]=5002.6697

	GV_t[0] =2907741.307
	GV_t[1] =-3312874.033
	GV_t[2] =1608101.838
	GV_t[3] =-433190.4871
	GV_t[4] =70624.81330
	GV_t[5] =-7116.620750
	GV_t[6] =432.51744
	GV_t[7] =-14.45911210
	GV_t[8] =0.2037

	a = 1.696985927
	b = -0.133372146
	c = 1.4
	f0 = 168.
	
	i = 0
	
	t=100.
	
	Tco = 190.564
	Pco = 45.992
	rco_t=0.1628
	
	Po=np.zeros(103)
	ro=np.zeros(103)
	vo=np.zeros(103)
	if (d):
		with open(d,'r') as f:
			for line in f:
				Po[i]=float(line.split()[1])
				ro[i]=float(line.split()[2])*0.001
				vo[i]=float(line.split()[11])
				i+=1
	return Po,ro,vo,a,b,c,f0,rco_t,Tco,Pco,t,A,GV_t

def dMt(ro,t,kj,rco):
	return dM1t(ro,t,kj,rco) * ( exp(dM2t(ro,t,kj,rco) + dM3t(ro,t,kj,rco) ) -1.0 )
def dM1t(ro,t,kj,rco):
	#print(f"t = {t}")
	res =exp(kj[0]+kj[1]/t)
	#print(f"dM1t = {res}")
	return res
def dM2t(ro,t,kj,rco):
	res =(ro**(0.1)) * (kj[2]+kj[3]/(t**(1.5)))
	#print(f"dM2t = {res}")
	return res
def dM3t(ro,t,kj,rco):
	#print(f"ro = {ro}, t = {t}, kj = {kj}, rco = {rco}")
	res = ((ro-rco)/rco)*(ro**0.5)*(kj[4] + kj[5]/t +kj[6]/(t**2))
	#print(f"dM3t = {res}")
	return res
def M1(GV,t):
	s=0
	for  i in range(0,9):
		s += GV[i] * (t**(-1 + (i)/3.))
	return s
def M2(a, b, c, f0,t):
	
	return a+(b*((c - log(t/f0))**2))
def Viscosity_t(ro,rco,GV,t,a,b,c,f0,A):
	v=M1(GV,t)+M2(a,b,c,f0,t)+dMt(ro,t,A,rco)
	return v/10000.

"""def main():
	Po,ro,vo,a,b,c,f0,rco,Tco,Pco,t,A,GV=import_data('data_methane_H_M.txt')	
	Res_t=np.zeros(103)

	for i in range(0,len(Po)):
		Res_t[i]=Viscosity_t(ro[i],rco,GV,t,a,b,c,f0,A)
	fig = plt.figure()
	for i in range(0,len(Po)):
		if i%5==0:
			plt.scatter(Po[i],Res_t[i],color='b')	
	
	plt.title("Вязкость метана по Хенли CH4")
	graph1=plt.plot(Po,Res_t)
	graph2=plt.plot(Po,vo)
	grid1=plt.grid(True)
	plt.ylabel('Вязкость, сП')
	plt.xlabel('Давление, бар')
	plt.show()"""
#main()

