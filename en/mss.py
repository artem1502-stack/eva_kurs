import matplotlib.pyplot as plt
import numpy as np
from math import log, tanh,exp
import H_M as hm
def import_data(d="data_test_mss.txt"):

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

	t=300
	kp=(Pco-Pto)/(Pcm-Ptm)
	bp=(Pco*Ptm-Pto*Pcm)/(Pcm-Ptm)
	To=t*(Tco-Tto)/(Tcm-Ttm)-(Tco*Ttm-Tto*Tcm)/(Tcm-Ttm)
	EEE=((Tcm/Tco)**(-1/6.))*((Pcm/Pco)**(2/3.))*((Mm/Mo)**(1/2.))
	with open(d,"r") as f:
		i=0
		N=102
		Po=np.zeros(N)
		P=np.zeros(N)
		
		ro=np.zeros(N)
		vo=np.zeros(N)
		for line in f:
			Po[i]=float(line.split()[1])
			ro[i]=float(line.split()[2])*0.001
			P[i]=(Po[i]+bp)*kp
			i+=1
	with open("data_test_mss.txt","r") as f:
		i=0
		for line in f:
			vo[i]=float(line.split()[11])
			i+=1
	return rcm,Pcm,Tcm,Mm,Tco,Pco,rco,Mo,t,EEE,To,Po,ro,vo,P
def Rr(ro,rco):
	return ro/rco
def alfa(Rr,Mo,Mm):
	ao = 1 + 0.007378*(Rr**1.847)*(Mo**0.5173)
	am = 1 + 0.007378*(Rr**1.847)*(Mm**0.5173)
	return ao,am
'''
def gen_res():
	
	with open("vsayt.txt","w") as f:
		for i in range(0,N):
			RoR=ro[i]/rco
			ao,am=alfa(RoR,Mo,Mm)
			p0=Po[i]*(Pco*ao)/(Pcm*am)
			T0=t*(Tco*ao)/(Tcm*am)
			print("{} {}".format(p0,T0),file=f)
#gen_res()
	'''
def no_alfa(d="data_methane_mss.txt"):
	rcm,Pcm,Tcm,Mm,Tco,Pco,rco,Mo,t,EEE,t1,Po,ro,vo,P=import_data(d)
	d=0

	V=np.zeros(len(Po))
	V_=np.zeros(len(Po))
	_,_,_,a,b,c,f0,_,_,_,_,A,GV=hm.import_data('data_methane_mss.txt')
	for i in range(0,len(Po)):
		ao,am=alfa(Rr(ro[i],rco),Mo,Mm)
		V[i]=EEE*hm.Viscosity_t(ro[i],rco,GV,t1,a,b,c,f0,A)*am/ao
		V_[i]=EEE*hm.Viscosity_t(ro[i],rco,GV,t1,a,b,c,f0,A)
		d+=(V[i]-vo[i])/V[i]
#		print(f'{i+1}. {V[i]}, {V_[i]}, {vo[i]}')

	fig = plt.figure()
	for i in range(0,len(Po)):
		if i%5==0:
			plt.scatter(P[i],V[i],color='b')	
			plt.scatter(P[i],V_[i],color='g')	
#	P=[i for i in range(1,103)]
	plt.title("Вязкость пропана по МСС C3H8")
	graph1=plt.plot(P,V)
	graph2=plt.plot(P,vo)
	graph3=plt.plot(P,V_)
	grid1=plt.grid(True)
	plt.ylabel('Вязкость, сП')
	plt.xlabel('Давление, бар')
	plt.show()
no_alfa()
