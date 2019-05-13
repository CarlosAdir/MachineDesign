# Material para colocar a teoria e algumas contas do projeto de eixo

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import interpolate
import sys

def delta(t):
	#return sp.Heaviside(t)
	return 0.5*sp.sign(t) + 0.5
def int_delta(t):
	return t*delta(t)
def int2_delta(t):
	return (t**2/2)*delta(t)

class Distribuicao:
	def __init__(self, q, inter):
		self.q = np.array(q)
		self.a = inter[0]
		self.b = inter[1]

class Forca:
	def __init__(self, F, p):
		self.F = np.array(F)
		self.p = p

class Momento:
	def __init__(self, M, p):
		self.M = np.array(M)
		self.p = p

class Eixo:
	def __init__(self, listas, L):
		self.listas = np.array(listas)
	def get_momento_inercia_y(self):
		t = sp.symbols('t')
		fd = 0
		for lista in self.listas:
			a, b, d = lista
			fd += d*(delta(t-a) - delta(t-b))
		I = np.pi*(fd**4)/64
		I = sp.lambdify(t, I, modules = ["numpy"])
		return I
	def get_momento_inercia_z(self):
		# Porque temos um eixo, que tem simetria radial
		return self.get_momento_inercia_y()
	def get_momento_polar(self):
		t = sp.symbols('t')
		fd = 0
		for lista in self.listas:
			a, b, d = lista
			fd += d*(delta(t-a) - delta(t-b))
		J = np.pi*(fd**4)/32
		J = sp.lambdify(t, J, modules = ["numpy"])
		return J


def calcula_diagramas_esforcos(esforcos, pm, Ls):
	ds, forcas, momentos = esforcos
	newForces = []
	for dist in ds:
		a, b = dist.a, dist.b
		p    = (b+a)/2
		F    = dist.q * (b-a)
		newForces.append(Forca(F, p))

	a, b = pm
	Fa, Fb = calcula_mancais(forcas + newForces, momentos, [a, b])
	
	forcas.append(Forca(Fa, a))
	forcas.append(Forca(Fb, b))

	Fs, Ms = get_funcoes(ds, forcas, momentos)

	return Fs, Ms, [Fa, Fb]

def plota_diagramas(Fs, Ms, Ls, nome = '1'):
	tamanho = 1200
	fator_sobra = 0.1 # 10% para cada borda
	titulo_superior = "DiagramaEsforcos"

	L0, Lf = Ls
	X = np.linspace(L0, Lf, tamanho)
	nomesF = ["N", "Vy", "Vz"]
	nomesM = ["T", "My", "Mz"]

	Xlabel = "$L \ (\mathrm{mm})$"
	Ylabel = [	"$N \ (\mathrm{N})$", \
				"$V \ (\mathrm{N})$", \
				"$T \ (\mathrm{N \cdot mm})$", \
				"$M \ (\mathrm{N \cdot mm})$"]
	titles = [	"Força normal", \
				"Força cortante", \
				"Momento torsor", \
				"Momento fletor"]

	N, Vy, Vz = Fs
	T, My, Mz = Ms

	YN, YVy, YVz = N(X), Vy(X), Vz(X)
	YT, YMy, YMz = T(X), My(X), Mz(X)

	minimos = [	min(YN), \
				min(min(YVy), min(YVz)), \
				min(YT), \
				min(min(YMy), min(YMz))]
	maximos = [	max(YN), \
				max(max(YVy), max(YVz)), \
				max(YT), \
				max(max(YMy), max(YMz))]
	fatores_multiplicadores = np.zeros(4)
	


	for i in range(4):
		if maximos[i] < 0:
			maximos[i] = 0
		if minimos[i] > 0:
			minimos[i] = 0
		dif = maximos[i] - minimos[i]
		while 1:
			if dif >= 1000:
				dif /= 1e+3
				maximos[i] /= 1e+3
				minimos[i] /= 1e+3
				fatores_multiplicadores[i] += 3
			else:
				break
		if dif == 0:
			dif = 1/fator_sobra
		minimos[i] -= fator_sobra * dif
		maximos[i] += fator_sobra * dif


	fig, axs = plt.subplots(2, 2)
	#fig.suptitle(titulo_superior, fontsize = 18)
	axs = axs.flatten()

	axs[0].set_title(titles[0])
	axs[0].set_xlim(L0, Lf)
	axs[0].set_ylim(minimos[0], maximos[0])
	axs[0].plot(X, YN/(10**fatores_multiplicadores[0]), label = nomesF[0], color = 'r')
	axs[0].grid()
	if fatores_multiplicadores[0] == 0:
		adicional = ""
	else:
		adicional = "$\cdot 10^{" + str(int(fatores_multiplicadores[0])) + "}$"
	axs[0].set(xlabel = Xlabel, ylabel = Ylabel[0] + adicional)

	axs[1].set_title(titles[1])
	axs[1].set_xlim(L0, Lf)
	axs[1].set_ylim(minimos[1], maximos[1])
	axs[1].plot(X, YVy/10**fatores_multiplicadores[1], label = nomesF[1], color = 'g')
	axs[1].plot(X, YVz/10**fatores_multiplicadores[1], label = nomesF[2], color = 'b')
	axs[1].legend()
	axs[1].grid()
	if fatores_multiplicadores[1] == 0:
		adicional = ""
	else:
		adicional = "$\cdot 10^{" + str(int(fatores_multiplicadores[1])) + "}$"
	axs[1].set(xlabel = Xlabel, ylabel = Ylabel[1] + adicional)

	axs[2].set_title(titles[2])
	axs[2].set_xlim(L0, Lf)
	axs[2].set_ylim(minimos[2], maximos[2])
	axs[2].plot(X, YT/10**fatores_multiplicadores[2], label = nomesM[0], color = 'r')
	axs[2].grid()
	if fatores_multiplicadores[2] == 0:
		adicional = ""
	else:
		adicional = "$\cdot 10^{" + str(int(fatores_multiplicadores[2])) + "}$"
	axs[2].set(xlabel = Xlabel, ylabel = Ylabel[2] + adicional)

	axs[3].set_title(titles[3])
	axs[3].set_xlim(L0, Lf)
	axs[3].set_ylim(minimos[3], maximos[3])
	axs[3].plot(X, YMy/10**fatores_multiplicadores[3], label = nomesM[1], color = 'b')
	axs[3].plot(X, YMz/10**fatores_multiplicadores[3], label = nomesM[2], color = 'g')
	axs[3].legend()
	axs[3].grid()
	if fatores_multiplicadores[3] == 0:
		adicional = ""
	else:
		adicional = "$\cdot 10^{" + str(int(fatores_multiplicadores[3])) + "}$"
	axs[3].set(xlabel = Xlabel, ylabel = Ylabel[3] + adicional)

	plt.subplots_adjust(top = 0.9, bottom=0.15, left = 0.15, right=0.90, wspace = 0.6, hspace = 0.6)
	plt.savefig("img/" + titulo_superior + nome + ".png")
	#plt.show()

	'''
	for i, F in zip(nomesF, Fs):
		Y = F(X)
		plt.plot(X, Y, label = str(i))
	plt.legend()
	plt.title("Diagrama de forcas")
	plt.show()
	plt.clf()
	for i, M in zip(nomesM, Ms):
		Y = M(X)
		plt.plot(X, Y, label = str(i))
	plt.legend()
	plt.title("Diagrama de momentos")
	plt.show()
	'''
	

def get_funcoes(ds, forcas, momentos):
	t = sp.symbols('t')
	N,  T  = 0, 0
	Vy, Vz = 0, 0
	My, Mz = 0, 0

	
	for dist in ds:
		N  -= dist.q[0]*(int_delta(t - dist.a)  - int_delta(t - dist.b))
		Vy += dist.q[1]*(int_delta(t - dist.a)  - int_delta(t - dist.b))
		Vz += dist.q[2]*(int_delta(t - dist.a)  - int_delta(t - dist.b))
		T  -= dist.q[0]*(int2_delta(t - dist.a) - int2_delta(t - dist.b))
		My += dist.q[2]*(int2_delta(t - dist.a) - int2_delta(t - dist.b))
		Mz += dist.q[1]*(int2_delta(t - dist.a) - int2_delta(t - dist.b))
	
	for forca in forcas:
		N  -= forca.F[0]*delta(t - forca.p)
		Vy += forca.F[1]*delta(t - forca.p)
		Vz += forca.F[2]*delta(t - forca.p)
		T  -= forca.F[0]*int_delta(t - forca.p)
		My += forca.F[2]*int_delta(t - forca.p)
		Mz += forca.F[1]*int_delta(t - forca.p)
	
	for momento in momentos:
		T  -= momento.M[0]*delta(t - momento.p)
		My += momento.M[1]*delta(t - momento.p)
		Mz += momento.M[2]*delta(t - momento.p)
	
	

	if N == 0:
		N = lambda t: np.zeros(len(t))
	else:
		N  = sp.lambdify(t, N,  modules = ["numpy"])
	if Vy == 0:
		Vy = lambda t: np.zeros(len(t))
	else:
		Vy = sp.lambdify(t, Vy, modules = ["numpy"])
	if Vz == 0:
		Vz = lambda t: np.zeros(len(t))
	else:
		Vz = sp.lambdify(t, Vz, modules = ["numpy"])
	if T == 0:
		T = lambda t: np.zeros(len(t))
	else:
		T = sp.lambdify(t, T, modules = ["numpy"])
	if My == 0:
		My = lambda t: np.zeros(len(t))
	else:
		My = sp.lambdify(t, My, modules = ["numpy"])
	if Mz == 0:
		Mz = lambda t: np.zeros(len(t))
	else:
		Mz = sp.lambdify(t, Mz, modules = ["numpy"])
	return [N, Vy, Vz], [T, My, Mz]

def calcula_mancais(forcas, momentos, pm):
	'''
	pm = posicao dos mancais = [a, b]
	'''
	a, b = pm
	#print("pm = " + str(pm))
	F, M = np.zeros(3), np.zeros(3)
	Fa, Fb = np.zeros(3), np.zeros(3)
	for forca in forcas:
		F += forca.F
	for forca in forcas:
		p = forca.p 
		M[1] += forca.F[2]*p
		M[2] += forca.F[1]*p
	for momento in momentos:
		M += momento.M
	Fy = F[1]
	Fz = F[2]
	My = M[1]
	Mz = M[2]


	Fa[1] = (+b*Fy - Mz)/(a-b)
	Fb[1] = (-a*Fy + Mz)/(a-b)
	Fa[2] = (+b*Fz - My)/(a-b)
	Fb[2] = (-a*Fz + My)/(a-b)

	return Fa, Fb

def integ(x, tck, constant=0):
	x = np.atleast_1d(x)
	Y = tck(x)
	out = np.zeros(x.shape, dtype=x.dtype)
	for n in range(1, len(out)):
		out[n] = out[n-1] + (x[n]-x[n-1])*(Y[n]+Y[n-1])/2
	out += constant
	return out

def get_deslocamentos(Ms, E, I, Ls, pm):
	Iy, Iz = I
	My, Mz = Ms
	L0, Lf = Ls
	a, b   = pm

	I = Iz # Simplificação feita por assumirmos que nao tem outra direcao

	tamanho = 1200
	# Mr = lambda t: np.sqrt(My(t)**2 + Mz(t)**2)
	Mr = Mz # COnsideracao que o momento em My é nulo
	X  = np.linspace(L0, Lf, 1200)
	Y  = Mr(X)/(E*I(X)) # Vemos que Mr está em N*mm

	f      = interpolate.interp1d(X, Y, kind='cubic')
	Y_int  = integ(X, f)
	int_f  = interpolate.interp1d(X, Y_int, kind='cubic')
	Y_int2 = integ(X, int_f)
	int_f2 = interpolate.interp1d(X, Y_int2, kind='cubic')

	C1 = -(int_f2(b) - int_f2(a))/(b-a) # As condicoes de contorno
	C2 = +(a*int_f2(b) - b*int_f2(a))/(b-a) # Condicoes de contorno

	Y_int  += C1
	Y_int2 += C1*X + C2


	return Y_int2, Y_int

def plota_deslocamentos(Ls, y, theta, nome = '1'):
	fator_sobra 	= 0.1 # 10% para cada borda
	titulo_superior = "Deslocamentos"
	Xlabel 			= "L ($mm$)"
	Ylabel 			= ["$\\theta$ (rad)", "y (mm)"]
	titles 			= ["Deslocamento angular", "Deslocamento da linha elastica",]
	L0, Lf 			= Ls
	
	minimos 		= [min(theta), min(y)]
	maximos 		= [max(theta), max(y)]

	# Para colocar as linhas horizontais que limitam as inclinacoes e deflexão maxima
	#print('max, min theta = ' + str(maximos[0]) + " " + str(minimos[0]))
	'''
	if maximos[0] > np.abs(inc_max) or minimos[0] < -np.abs(inc_max):
		color1 = 'r'
	else:
		color1 = 'g'
	if maximos[1] > np.abs(defl_max) or minimos[1] < -np.abs(defl_max):
		color2 = 'r'
	else:
		color2 = 'b'
	'''
	color1, color2 = 'g', 'b'

	for i in range(2):
		difs = maximos[i] - minimos[i]
		if difs == 0:
			difs = 10
		minimos[i] -= fator_sobra*difs
		maximos[i] += fator_sobra*difs

	X = np.linspace(L0, Lf, len(y))

	fig, axs = plt.subplots(2, 1)
	fig.suptitle(titulo_superior, fontsize = 18)
	axs = axs.flatten()

	axs[0].set_title(titles[0])
	axs[0].set_xlim(L0, Lf)
	axs[0].set_ylim(minimos[0], maximos[0])
	axs[0].plot(X, theta, color = color1)
	#shape1 = patches.Rectangle((L0, -inc_max), Lf-L0, 2*inc_max, color = '.75')
	#axs[0].add_patch(shape1)
	axs[0].grid()
	axs[0].set(xlabel = Xlabel, ylabel = Ylabel[0])

	axs[1].set_title(titles[1])
	axs[1].set_xlim(L0, Lf)
	axs[1].set_ylim(minimos[1], maximos[1])
	axs[1].plot(X, y, color = color2)
	#shape2 = patches.Rectangle((L0, -defl_max), Lf-L0, 2*defl_max, color = '.75')
	#axs[1].add_patch(shape2)
	axs[1].grid()
	axs[1].set(xlabel = Xlabel, ylabel = Ylabel[1])

	plt.subplots_adjust(top = 0.85, bottom=0.15, left = 0.15, right=0.90, wspace = 0.4, hspace = 0.45)
	plt.savefig("img/" + titulo_superior + nome + ".png")
	#plt.show()

if __name__ == "__main__":
	try:
		args = sys.argv
		teste = int(args[1])
	except:
		raise ValueError("Hey, coloque algo do tipo: python3 arquivo 3")

	n 				= 2.84 # Fator de segurança
	g 				= 9.85
	rho  			= 7.85e-6 # kg/mm^3
	E    			= 206*1e+3 # Em MPa, ou também N/mm²
	G 				= 79.3*1e+3
	phi_n 			= 20*(np.pi/180) # Angulo em radianos
	Lq 				= 298 # mm
	La 				= 162 # mm
	q1 				= 726*g/(2*Lq) # Em N/mm
	q2 				= 160*g/(2*La) # Em N/mm
	L1, L2, L3		= 310, 100, 100		# mm
	L4, L5, L6		= 119.5, 113, 123 	# mm
	L7, L8, L9		= 123, 240, 113	# mm
	L10, L11, L12	= 113, 113, 50 	# mm
	P1 				= 4*49.5*g	# São 4 massas planetarias, cada uma com 49.5 quilos, vezes a gravidade
	P2 				= 22*g		# O peso da solar, que tem 22 quilogramas
	P3 				= 279.7*g
	P4 				= 69.9*g
	P5 				= 150.6*g
	P6 				= 16.7*g
	Ft2 			= 53476.06
	Ft4 			= 33422.54
	T1, T2, T3, T4  = 80214.1e+3, 16042.82e+3, 8021.4e+3, 2673.8e+3 # N*mm
	Fn2 			= Ft2*np.tan(phi_n)
	Fn4 			= Ft4*np.tan(phi_n)
	Sy 				= 450  # MPa, limite de escoamento do eixo
	Sfratura 		= 515  # MPa, limite de resistencia à fratura do eixo
	Se 				= (Sfratura/2) # Se' vezes todos os vatores de correcao, concentrador de tensoes, temperatura, tamanho, etc


	# Inclinação no mancal
	#   Rolos cilindricos deve ser menor que 0.001 rad
	#   Rolos conicos deve ser menor que 0.0005 rad
	#   Rolos de sulco profundo e pista profunda deve ser menor que 0.004
	#   0.0087 para mancais de esfera
	#inc_max = 0.0005 # Radianos, para mancal de rolo
	inc_max = 0.0087 # Radianos
	defl_max = 0.120 # As deflexoes nas engrenagens nao devem exceder 120 um

	distibuicoes = []
	forcas 		 = []
	momentos 	 = []

	print("teste = " + str(teste))
	if teste == 1:
		d    	= 172.77 # Diametro
		P1 		= [0, -P1, 0]# Os pesos das planetarias
		qeixo  	= [0, -np.pi*rho*(d**2)*g/4, 0] # A distribuicao de forca que o peso do proprio eixo tem
		qacopl  = [0, -q1, 0] # Distribuicao do 
		TE1		= [T1, 0, 0]

		L    	= L1+L2+L3  # O comprimento total da viga/eixo
		a, b 	= L1, L1+L2 # A posicao dos mancais
		

		distibuicoes.append(Distribuicao(qeixo,  [0, L]))
		distibuicoes.append(Distribuicao(qacopl, [0, Lq]))
		forcas.append(Forca(P1, L))
		momentos.append(Momento(TE1, 0))

	elif teste == 2:
		d 		= 101.14
		
		L 		= L4+L5+L6
		a, b 	= L4, L4+L5+L6

		P2 		= [0, -P2, 0] 	 # Peso da solar
		FE2 	= [0, -Fn2, -Ft2] # A reação que a engrenagem 3 faz na 2
		P3 		= [0, -P3, 0] 		 # O peso da engrenagem 3
		qeixo  	= [0, 	-np.pi*rho*(d**2)*g/4,		0] # A distribuicao de forca que o peso do proprio eixo tem
		TE2		= [T2, 0, 0]
		TE3		= [-T2, 0, 0]


		distibuicoes.append(Distribuicao(qeixo,  [0, L]))
		forcas.append(Forca(P2, 0))
		forcas.append(Forca(FE2, L4+L5))
		forcas.append(Forca(P3, L4+L5))
		momentos.append(Momento(TE2, 0))
		momentos.append(Momento(TE3, L4+L5))



	elif teste == 3:
		d = 102.72
		
		L 		= L7+L8+L9
		a, b 	= 0, L

		FE2 	= [0, +Fn2, +Ft2] # Peso da solar
		FE4 	= [0, +Fn4, -Ft4] # A reação que a engrenagem 3 faz na 2
		P4 		= [0, -P4, 0]
		P5 		= [0, -P5, 0]
		qeixo  	= [0, 	-np.pi*rho*(d**2)*g/4,		0] # A distribuicao de forca que o peso do proprio eixo tem
		TE4		= [T3, 0, 0]
		TE5		= [-T3, 0, 0]

		distibuicoes.append(Distribuicao(qeixo,  [0, L]))
		forcas.append(Forca(FE2, L7))
		forcas.append(Forca(FE4, L7+L8))
		forcas.append(Forca(P4,  L7))
		forcas.append(Forca(P5,  L7+L8))
		momentos.append(Momento(TE4, L7))
		momentos.append(Momento(TE5, L7+L8))

		



	elif teste == 4:
		d = 56.26

		L = L10+L11+L12
		a, b = 0, L10+L11

		FE4 	= [0, -Fn4, +Ft4]
		P6 		= [0, -P6, 0]
		qeixo  	= [0, 	-np.pi*rho*(d**2)*g/4,		0] # A distribuicao de forca que o peso do proprio eixo tem
		qacopl 	= [0, -q2, 0]
		TE6		= [T4, 0, 0]
		TEF		= [-T4, 0, 0]

		distibuicoes.append(Distribuicao(qeixo,  [0, L]))
		distibuicoes.append(Distribuicao(qacopl, [L-Lq, L]))
		forcas.append(Forca(P6, L10))
		forcas.append(Forca(FE4, L10))
		momentos.append(Momento(TE6, L10))
		momentos.append(Momento(TEF, L))

	elif teste == 5:
		d = 70
		C0, C1, C2, C3, C4 = 10, 40, 139-2*40, 40, 150-139
		L = C0 + C2 + C2 + C3 + C4
		a, b = 0, 10

		Fs = 100.5*1e+3 # N
		q  = Fs/(C1+C3)
		q1 	= [0, -q, 0]
		q2  = [0, -q, 0]

		distibuicoes.append(Distribuicao(q1, [C0, C1]))
		distibuicoes.append(Distribuicao(q2, [C0+C1+C2, C0+C1+C2+C3]))
	else:
		print("Coloque um valor de teste valido!")
		sys.exit(0)


	print("L = %.2f mm" % L)

	Fs, Ms = calcula_diagramas_esforcos(distibuicoes, forcas, momentos, [a, b], [0, L])
	T, My, Mz = Ms

	N 		= 1200
	X 		= np.linspace(0, L, N)
	Torqus 	= T(X)
	Momens 	= Mz(X)

	d_tor	= ((32*180*20)/(np.pi**2)) * (np.abs(Torqus)/G)
	d_tor 	= (d_tor)**(1/3) 

	momentosss 	= (np.abs(Momens**2 + Torqus**2))**(1/2)
	d_est	= (32*n*momentosss/(np.pi*Sy))**(1/3)


	Kf  	= 1.0 # Fator para flexao
	Kfs 	= 1.0 # Fator para torcao
	Ma, Ta 	= 0, 0 # Fletor e Torsor alternado
	A 		= np.sqrt(4*(Kf*Ma)**2 		+ 3*(Kfs*Ta)**2)
	B 		= np.sqrt(4*(Kf*Momens)**2 	+ 3*(Kfs*Torqus)**2)
	d_fad 	= (16*n*np.sqrt(np.abs((A/Se)**2 + (B/Sy)**2))/np.pi)**(1/3)

	'''
	plt.plot(X, d_fad, color = 'b')
	plt.show()
	'''

	print('max tor = %.2f' % max(d_tor))
	print('max est = %.2f' % max(d_est))
	print('max fad = %.2f' % max(d_fad))



	dmin, dmax = 10, 500
	N, Nmax = 1, 15
	while 1:
		d = (dmin+dmax)/2
		eixo_util = Eixo([[0, L, d]], L)
		Iz = eixo_util.get_momento_inercia_z()
		Iy = eixo_util.get_momento_inercia_y()
		J = eixo_util.get_momento_polar()
		y, theta 	= get_deslocamentos([My, Mz], E, Iz, [0, L])
		#print('y   = ' + str(y))
		#print('t   = ' + str(theta))
		#print('|y| = ' + str(k))
		#print('max(k) = ' + str(max(k)))
		theta_max 	= max(np.abs(theta))
		y_max 		= max(np.abs(y))
		#print("[t_max, y_max] = [%.1f, %.1f]\n" % (theta_max, y_max))
		if (y_max < defl_max) and (theta_max < inc_max):
			#print("[dmin, dmax] = [%.1f, %.1f]" % (dmin, dmax))
			dmax = d
	#		print("ymax = %.5f, dmax = %.2f" % (y_max, dmax))
		else:
			dmin = d
		N += 1
		if N == Nmax:
			break
	d_def = dmax
	print("max def = %.2f" % d_def)

	diametro_minimo = max([max(d_tor), max(d_est), max(d_fad), d_def])
	print("\nminimo necessario = %.2f" % diametro_minimo)


	eixo_util = Eixo([[0, L, diametro_minimo]], L)
	Iz = eixo_util.get_momento_inercia_z()
	#Iy = eixo_util.get_momento_inercia_y()
	#J = eixo_util.get_momento_polar()
	y, theta 	= get_deslocamentos([My, Mz], E, Iz, [0, L])
	plota_graficos(Fs, Ms, [0, L])
	plot_deslocamentos([0, L], y, theta)


# Em resumo, o projeto de eixos deve incluir principalmente os seguintes topicos
# * Analise das Tensões e Resistência à Fadiga
# * Deflexões
# * Velocidades criticas

# As deflexoes nas engrenagens nao devem exceder 120 um
# As deflexoes angulares nos mancais auto alinhantes nao devem superar 0.04º ou 0.000698
# Se há cargas axiais deveremos usar apenas um rolamento axial em cada direção da carga.
# A primeira frequência natural do eixo deve ser no mínimo três vezes maior do que a freqüência forçada esperada em serviço (de preferencia 10 vezes)
# Se houverem ressaltos, momento de inercia muda, o que muda a equação da linha elástica
# Os deslocamentos angulares por torção devem ser considerados com base nos requisitos associados à frequência natural torcional e às limitações das deflexões torcionais

# Restrições geométricas
# # Mancais
# # # Inclinação no mancal
# # # # Rolos cilindricos deve ser menor que 0.001 rad
# # # # Rolos conicos deve ser menor que 0.0005 rad
# # # # Rolos de sulco profundo e pista profunda deve ser menor que 0.004
# # # # 0.0087 para mancais de esfera

# As deflexoes por flexao do eixo podem ser calculadas pelos seguintes procedimentos
# a) Integração da linha elástica
# b) Calculo da area do diagrama de momentos fletores
# c) Integracao por funcoes singulares
# d) Integração gráfica
# e) Integracao numerica
# f) Método da superposição
# g) Método dos elementos finitos
# O calculo das inclinacoes e deflexoes para um eixo com varios diametros sao mais complicados, uma vez que ambos momentos variam ao longo do comprimento do eixo


# d²y/dx² = M/(EI)
# theta = dy/dx = int M/(EI) dx
# Basicamente tem que pegar o diagrama de momento, e calcular o y e o theta resolvendo a equação diferencial. 
# Daí pegamos o valor da deflexão máxima, e das inclinações nos mancais e verificamos o resultado.


# Deflexao angular
# No slide 59 do toninho fala: "Na ausência de norma específica, ..."
# Na literatura inglesa, é comum a indicação da deformação admissivel de 1 grau em 20 x d (em polegada)
# Para árvores longas, a variação admissível é de 0,25 a 0,50 graus por metro linear