# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:08:04 2019

@author: carlos
"""

import numpy as np
from scipy.interpolate import interp1d, interp2d
import matplotlib.pyplot as plt
from Engrenagens_lib import EngExtDentesRetos, EngExtDentesHelicoidais


# Esse import é só usado a __main__
from Materiais_lib import Material

class FatoresDenteFundamental:
	'''
	DESCRIPTION
		Essa classe determina os esforços fundamentais 
		As variáveis comuns para as duas tensões:
		* b = largura de face do membro mais fino, em mm
		* Ft = Força tangencial, em N
		* K0 = Fator de sodbrecarga
		* Kv = Fator dinâmico
		* Ks = Fator de tamanho (Size factor)
		* KH = Fator de carregamento
		  - KHpf
		  - KHpf
		  - KHpf
		  - KHpf
		  - KHpf
	'''
	@staticmethod
	def K0():
		# Porque a fonte de potencia e a maquina acionada sao uniformes
		return 1
	@staticmethod
	def Kv(Qv, v):
		'''
		Esse é o fator dinâmico. Precisa do acabamento da peça e a velocidade no diâmetro primitivo
		Qv é um valor numério de 5 a 11
		v é a velocidade, em m/s
		Pag 15 da norma
		'''
		if not (5 <= Qv <= 11):
			raise ValueError("O valor de Qv deve estar no intervalo [5, 11]")
		B = 0.25*(12-Qv)**(2/3.)
		A = 50+56*(1-B)
		Vmax = ((A+Qv-3)**2)/200 # Vmax é retornado em m/s
		if v > Vmax:
			raise ValueError("No calculo do fator dinamico Kv, a velocidade aplicada é maior que a velocidade maxima permitida")
		return ((A+np.sqrt(200*v))/A)**B
	@staticmethod
	def __Y(N, phi = 20):
		'''
		Precisa dessa funcao para a funcao Ks e para YJ
		É conhecido como fator de Lewis
		'''
		if phi != 20: 
			raise ValueError("Nao foi possivel calcular o Y pois o angulo nao eh 20 graus e eu nao sei calcular quando acontece isso")

		# É obriatorio que o anulo de pressao seja 20° para usar essa equaçao
		tabela_Ys = {	12:0.245, 13:0.261, 14:0.277, 15:0.290, 16:0.296, 
						17:0.303, 18:0.309, 19:0.314, 20:0.322, 21:0.328,
						22:0.331, 24:0.337, 26:0.346, 28:0.353, 30:0.359,
						34:0.371, 38:0.384, 43:0.397, 50:0.409, 60:0.422,
						75:0.435, 100:0.447, 150:0.460, 300:0.472, 400:0.480,
						10000:0.485}
		if N < 12:
			raise ValueError("Nao existe Y para N < 12")
		elif N > 10000:
			return tabela_Ys[10000]
		n = len(tabela_Ys)
		Ns = np.zeros(n)
		Ys = np.zeros(n)
		for i, e in enumerate(tabela_Ys):
			Ns[i] = e
		Ns.sort()
		for i in range(n):
			Ys[i] = tabela_Ys[Ns[i]]
		f = interp1d(Ns, Ys)
		return f(N)
	@staticmethod
	def Ks(N, b, m):
		'''
		N é o numero de dentes
		b é a largura do dente, supostamente
		m é o módulo do dente
		
		'''
		if N == None:
			raise ValueError("Numero de dentes nao pode ser nulo!")
		if b == None:
			raise ValueError("A largura de face nao pode ser nula!")
		if m == None:
			raise ValueError("O modulo nao pode ser nulo")
		Y = FatoresDenteFundamental.__Y(N)
		Ks = 0.8433*(m*b*np.sqrt(Y))**0.0535
		if Ks < 1:
			return 1
		else:
			return Ks
	@staticmethod
	def KH(b, d):
		'''
		b = largura de face
		d = diametro
		'''
		razao = b/(10*d)
		if razao < 0.05:
			razao = 0.05
		# Slide 81 de PROJETO DE ENGRENAGENS DE DENTES RETOS 2018
		nomes = ["Engrenamento aberto", \
				 "Unidades fechadas, comerciais", \
				 "unidades fechadas, de precisão", \
				 "Unidades de engrenagens fechadas, extraprecisas"]
		# Essa tabela abaixo é para shigley
		# abcs = { nomes[0]: [0.247, 0.0167, -0.765e-4],
		#		  nomes[1]: [0.127, 0.0158, -0.930e-4],
		#		  nomes[2]: [0.0675, 0.0128, -0.926e-4],
		#		  nomes[3]: [0.00360, 0.0102, -0.822e-4]}
		# Essa atabela abaixo, pag 22 da AGMA 2101 C95, tabela 2 - constantes empiricas A, B, C
		abcs = { nomes[0]: [2.48e-1, 0.657e-3, -1.186e-7],
				 nomes[1]: [1.27e-1, 0.622e-3, -1.69e-7],
				 nomes[2]: [0.675e-1, 0.504e-3, -1.44e-7],
				 nomes[3]: [0.380e-1, 0.402e-3, -1.27e-7]}
		
		decisao = nomes[3] # Pois escolhemos o QV como 11, entao é preciso
		A, B, C = abcs[decisao]
		
		# Agora sao 5 constantes
		KHmc = 1 # pois escolhemos dentes sem coroamento
		if b <= 25:
			KHpm = razao - 0.025
		elif 25 < b <= 432:
			KHpm = razao - 0.0375 + 0.492*1e-3*b
		elif 432 < b < 1020:
			KHpm = razao - 0.1109 + 0.815*1e-3*b - 0.353*1e-6*b**2
		else:
			raise ValueError("A largura da engrenagem excede 1 metro! Ta errado!")
		KHpf = 1.1 # Pois o pinhão é montado fora do centro do mancal
		KHma = A+B*b+C*b**2
		KHe  = 1 # Nao sabemos bem o porquê ainda
		
		Kh = 1+KHmc*(KHpf*KHpm+KHma*KHe)
		
		return Kh

class FatoresDenteFundamentalFlexao(FatoresDenteFundamental):
	'''
	DESCRIPTION
		As variáveis utilizadas apenas na tensão de flexão no pé do dente
		* KB = fator de espessura da "gengiva"
		* YmJ = fator geométrico para esforço de flexão	 
	'''
	@staticmethod
	def KB():
		# Tem relação com a gengiva do dente de engrenagem
		# Para quando o tamanho da gengiva for 1.2 vezes maior que o tamanho do dente, o valor é 1
		Kb = 1
		return Kb
	@staticmethod
	def mt(mn, gamma):
		'''
		mn = módulo
		gamma = angulo transverssal
		'''
		_mt = mn/np.cos(gamma)
		return _mt
	@staticmethod
	def YJ(N1, N2):
		X = np.log(np.array([17, 25, 35, 50, 85, 170, 1000]))
		Y = np.log(np.array([20, 30, 50, 80, 275]))
		Z = [[0.310, 0.360, 0.394, 0.420, 0.445], \
			 [0.320, 0.366, 0.406, 0.430, 0.460], \
			 [0.327, 0.375, 0.418, 0.442, 0.470], \
			 [0.330, 0.381, 0.429, 0.456, 0.490], \
			 [0.338, 0.390, 0.439, 0.467, 0.500], \
			 [0.342, 0.399, 0.450, 0.480, 0.520], \
			 [0.347, 0.405, 0.460, 0.490, 0.530]]
		f = interp2d(Y, X, Z, kind = 'cubic')
		return float(f(np.log(N2), np.log(N1)))

class FatoresDenteFundamentalContato(FatoresDenteFundamental):
	'''
	DESCRIPTION
		As variáveis utilizadas apenas na tensão de contato
		* ZR = fator de condição de superficie para "pitting resistence"
		* ZI = fator geométrico para "pitting resistence"
		* dw1 = diâmetro de operação do pinhão, mm
	'''
	@staticmethod
	def ZI(phi_t, mG):
		# Se dentes retos
		# mG é a razao de dentes, sempre maior ou igual a 1
		# mG = NG/NP, segundo shigley, pag 678
		# phi_t é o angulo transversal, phi_t = arctan(tan(phi_n)/cos(gamma))
		mN = 1 # É 1 pois é dentes retos
		COS = np.cos(phi_t)
		SIN = np.sin(phi_t)
		zi = (SIN*COS/(2*mN))*(mG/(mG+1))
		return zi
	@staticmethod
	def ZR():
		# Se soubermos o efeito prejudicial causado pelo acabamento superficial, o fator é 1
		return 1
	@staticmethod
	def ZE(Mat1, Mat2):
		# Tem relação com o material, o tanto de deformação
		E1 = Mat1.get_elasticidade()
		v1 = Mat1.get_poisson()
		E2 = Mat2.get_elasticidade()
		v2 = Mat2.get_poisson()
		ZE = np.sqrt(1/(np.pi*(     (1-v1**2)/E1  +  (1-v2**2)/E2     )))
		return ZE
















class FatoresDentePermitido:
	'''
	DESCRIPTION
		Y_theta = fator de temperatura
		Yz = fator de reabilidade
		S  = fator de segurança, utilizaremos o mesmo para ambos casos
	'''
	@staticmethod
	def Yz():
		# Confiabilidade, usaremos como 99% de confiabilidade, de modo que nosso fator é 1
		return 1
	@staticmethod
	def Ytheta():
		# Assumimos que não passa a temperatura de 120 ºC, o que nos garante o fator de 1
		return 1
	@staticmethod
	def S():
		# O fator de segurança escolhido para o projeto
		return 2.84

class FatoresDentePermitidoFlexao(FatoresDentePermitido):
	'''
	DESCRIPTION
		* sigmaFP = tensão de flexão permitida
		* YN = fator de tensão ciclica em flexao
		* SF = fator de segurança para flexão
	'''
	@staticmethod
	def sigmaFP(HB):
		# Aqui assumimos que é grau 2
		return 0.703*HB+113 # retorna em MPa
	@staticmethod
	def YN():
		# Por enquanto que não sabemos o valor de YN
		return 1
	@staticmethod
	def SF():
		# Aqui vem o fator de segurança
		return FatoresDentePermitido.S() # Usaremos o padrão, o mesmo para flexão e contato

class FatoresDentePermitidoContato(FatoresDentePermitido):
	'''
	DESCRIPTION
		As variáveis utilizadas apenas na tensão de contato
		* sigmaHP = tensão de contato permitido
		* ZN = fator de tensão ciclica para "pitting resistence"
		* ZW = fator razão de dureza para "pitting resistence"
		* SH = fator de segurança para "pitting resistence"
	'''
	@staticmethod
	def SH():
		# Aqui vem o fator de segurança
		return (FatoresDentePermitido.S())**0.5
	@staticmethod
	def sigmaHP(HB):
		# Aqui assumimos que é grau 2
		return 2.41*HB+237
	@staticmethod
	def ZN():
		# Tem a ver com fadiga e a quantidade de ciclos da engrenagem.
		# Pagina 755 do shigley
		return 1 # Por enquanto é isso, porque precisamos pegar informações de um gráfico
	@staticmethod
	def ZW():
		# Fator razaão de dureza para "pitting resistence"
		return 1














class EsforcoNoDente:
	@staticmethod
	def tensao_flexao(E1, E2, verbose = False):
		Ft = E1.get_Wt()
		b  = E1.get_largura()

		K0 = FatoresDenteFundamentalFlexao.K0()
		Kv = FatoresDenteFundamentalFlexao.Kv(Qv = 11, v = E1.get_velocidade())
		Ks = FatoresDenteFundamentalFlexao.Ks(N = E1.get_dentes(), b = E1.get_largura(), m = E1.get_modulo())
		KH = FatoresDenteFundamentalFlexao.KH(b = E1.get_largura(), d = E1.get_diametro())
		if type(E1) == EngExtDentesRetos:
			mt = FatoresDenteFundamentalFlexao.mt(E1.get_modulo(), 0) # Poderia ser E2.get_modulo() pois os modulos sao iguais
		elif type(E1) == EngExtDentesHelicoidais:
			mt = FatoresDenteFundamentalFlexao.mt(E1.get_modulo(), E1.get_angulo_helice_degrees())
		else:
			raise ValueError("WHat the hell! Nem deveria entrar aqui")

		KB = FatoresDenteFundamentalFlexao.KB()
		YJ = FatoresDenteFundamentalFlexao.YJ(E2.get_dentes(), E1.get_dentes())

		if verbose:
			print("Tensoes nos dentes entre:")
			print("Flexao:   [K0, Kv, Ks, KH, mt, KB, YJ] = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]" % (K0, Kv, Ks, KH, mt, KB, YJ))
		return (Ft*K0*Kv*Ks*KH*KB)/(b*mt*YJ)
	@staticmethod
	def tensao_flexao_admissivel(E1, E2, verbose = False):

		HB = E1.get_material().get_brinell()

		sigmaFP = FatoresDentePermitidoFlexao.sigmaFP(HB)
		YN = FatoresDentePermitidoFlexao.YN()
		Ytheta = FatoresDentePermitidoFlexao.Ytheta()
		Yz = FatoresDentePermitidoFlexao.Yz()
		SF = FatoresDentePermitidoFlexao.SF()

		if verbose:
			print("Flex max: [sigmaFP, YN, Ytheta, Yz, SF] = [%.2f, %.2f, %.2f, %.2f, %.2f]" % (sigmaFP, YN, Ytheta, Yz, SF))

		return (sigmaFP*YN)/(Ytheta*Yz)  #########

	@staticmethod
	def tensao_contato(E1, E2, verbose = False):
		debug = True
		if E1.get_dentes() > E2.get_dentes():
			Ep = E2
			Eg = E1
		else:
			Ep = E1
			Eg = E2
		Ft = E1.get_Wt()
		b  = E1.get_largura()
		dw1 = Ep.get_diametro()

		if debug:
			print("Ft = " + str(Ft))
			print("b  = " + str(b))

		K0 = FatoresDenteFundamentalFlexao.K0()
		Kv = FatoresDenteFundamentalFlexao.Kv(Qv = 11, v = E1.get_velocidade())
		Ks = FatoresDenteFundamentalFlexao.Ks(N = E1.get_dentes(), b = E1.get_largura(), m = E1.get_modulo())
		KH = FatoresDenteFundamentalFlexao.KH(b = E1.get_largura(), d = E1.get_diametro())

		mG = Eg.get_dentes()/Ep.get_dentes()
		ZI = FatoresDenteFundamentalContato.ZI(E1.get_angulo_pressao_degrees(), mG)
		ZR = FatoresDenteFundamentalContato.ZR()
		ZE = FatoresDenteFundamentalContato.ZE(E1.get_material(), E2.get_material())
		if verbose or debug:
			print("Contato:  [K0, Kv, Ks, KH, ZE, ZR, ZI] = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]" % (K0, Kv, Ks, KH, ZE, ZR, ZI))
		return ZE*np.sqrt(  (Ft*K0*Kv*Ks*KH*ZR)/(dw1*ZI*b)  )
	@staticmethod
	def tensao_contato_admissivel(E1, E2, verbose = False):
		debug = False
		HB = E1.get_material().get_brinell()
		if debug:
			print("HB = " + str(HB))
		sigmaHP = FatoresDentePermitidoContato.sigmaHP(HB)
		ZN = FatoresDentePermitidoContato.ZN()
		ZW = FatoresDentePermitidoContato.ZW()
		SH = FatoresDentePermitidoContato.SH()
		Ytheta = FatoresDentePermitidoContato.Ytheta()
		Yz = FatoresDentePermitidoContato.Yz()
		if verbose or debug:
			print("Cont_max: [sigmaHP, ZN, Ytheta, Yz, ZW, SH] = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f]" % (sigmaHP, ZN, Ytheta, Yz, ZW, SH))
		return (sigmaHP*ZN*ZW)/(Ytheta*Yz)
		
	@staticmethod
	def get_tensoes_dente(E1, E2, verbose = False):
		tensoes = []
		tensoes.append(EsforcoNoDente.tensao_flexao(E1, E2, verbose))
		tensoes.append(EsforcoNoDente.tensao_flexao_admissivel(E1, E2, verbose))
		tensoes.append(EsforcoNoDente.tensao_contato(E1, E2, verbose))
		tensoes.append(EsforcoNoDente.tensao_contato_admissivel(E1, E2, verbose))
		SF = tensoes[1]/tensoes[0]
		SH = (tensoes[3]/tensoes[2])**2
		tensoes.append(SF)
		tensoes.append(SH)
		return tensoes








def calculc(n, b = 5):
	m = modulo = 10 # modulo em metros
	#Ns, Na, Np = 20, 80, 30
	Ns, Na, Np = 16, 64, 24
	M = 400*1000*	60/(2*np.pi*50) # 400 kW, a 50 rpm, torque
	phi_n = 20*np.pi/180
	phi_t = phi_n
	Ft = M*1000/(4*m*(Np+Ns))
	Fr = Ft*np.tan(phi_t)
	mat = Material("steel 1045")

	n_in = 50 # rpm
	w    = n_in * 2*np.pi/60

	ds, da, dp = Ns*m, Na*m, Np*m
	v    = ((w * dp/2)*(ds+dp)/(dp))/1000 # Dividimos por 1000 para transformar em m/s

	HB   = 500
	b *= np.pi*m

	#print('Modulo = %.3f m' % (m/1000))
	#print("[Ns, Na, Np] = " + str([Ns, Na, Np]))
	print("[ds, da, dp] = " + str([ds, da, dp]) + " mm")
	ms, ma, mp = 7.85*np.pi*(ds**2/4)*b/1e+6, 7.85*np.pi*(da**2/4)*b/1e+6, 7.85*np.pi*(dp**2/4)*b/1e+6
	print("[ms, ma, mp] = [%.1f, %.1f, %.1f] kg" % (ms, ma, mp))
	#print("Torque = %.4f N*m" % M)
	#print("Angulo = %.2f o" % (phi_t*180/np.pi))
	#print("Ft = %.1f kN" % (Ft/1000))
	#print("Fr = %.1f kN" % (Fr/1000))
	#print("v  = %.2f m/s" % v)
	print("b  = %.3f m" % (b/1000))

	if n == 1:
		mensagem = "Calculo dos fatores da planeta engrenando com a solar:"
		Npar, Nimpar = Ns, Np
		dpar, dimpar = ds, dp
	elif n == 2:
		mensagem = "Calculo dos fatores da planeta engrenando com a anular:"
		Npar, Nimpar = Na, Np
		dpar, dimpar = da, dp
	elif n == 3:
		mensagem = "Calculo dos fatores da solar engrenando com a planeta:"
		Npar, Nimpar = Np, Ns
		dpar, dimpar = dp, ds
	elif n == 4:
		mensagem = "Calculo dos fatores da anular engrenando com a planeta:"
		Npar, Nimpar = Np, Na
		dpar, dimpar = dp, da
	print(mensagem)

	dpinhao = dpar if dpar > dimpar else dimpar
	dw1 = dpinhao
	#print("\n\n\n")


	#print("Flexao")
	K0 = FatoresDenteFundamentalFlexao.K0()
	Kv = FatoresDenteFundamentalFlexao.Kv(Qv = 11, v = v)
	Ks = FatoresDenteFundamentalFlexao.Ks(N = Nimpar, b = b, m = m)
	KH = FatoresDenteFundamentalFlexao.KH(b = b, d = dimpar)
	mt = FatoresDenteFundamentalFlexao.mt(m, 0) # é necessario que esteja em mm
	KB = FatoresDenteFundamentalFlexao.KB()
	YJ = FatoresDenteFundamentalFlexao.YJ(Npar, Nimpar)
	#print("[K0, Kv, Ks, KH, mt, KB, YJ] = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]" % (K0, Kv, Ks, KH, mt, KB, YJ))
	sigma_flexao = Ft*K0*Kv*Ks*KH*KB/(b*mt*YJ)
	
	sigmaFP = FatoresDentePermitidoFlexao.sigmaFP(HB)
	YN = FatoresDentePermitidoFlexao.YN()
	Ytheta = FatoresDentePermitidoFlexao.Ytheta()
	Yz = FatoresDentePermitidoFlexao.Yz()
	SF = FatoresDentePermitidoFlexao.SF()
	#print("[sigmaFP, YN, Ytheta, Yz, SF] = [%.2f, %.2f, %.2f, %.2f, %.2f]" % (sigmaFP, YN, Ytheta, Yz, SF))
	sigma_max_flexao = (sigmaFP*YN)/(Ytheta*Yz)  #########

	print("sigma_flexao = [%.3f / %.3f] MPa" % (sigma_flexao, sigma_max_flexao))

	#print("Contato")
	K0 = FatoresDenteFundamentalContato.K0()
	Kv = FatoresDenteFundamentalContato.Kv(Qv = 11, v = v)
	Ks = FatoresDenteFundamentalContato.Ks(N = Nimpar, b = b, m = m)
	KH = FatoresDenteFundamentalContato.KH(b = b, d = dimpar)
	ZI = FatoresDenteFundamentalContato.ZI(phi_n, 2)
	ZR = FatoresDenteFundamentalContato.ZR()
	ZE = FatoresDenteFundamentalContato.ZE(mat, mat)
	#print("[K0, Kv, Ks, KH, ZE, ZR, ZI] = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]" % (K0, Kv, Ks, KH, ZE, ZR, ZI))
	sigma_contato = ZE*np.sqrt(Ft*K0*Kv*Ks*KH*ZR/(dw1*ZI*b))

	sigmaHP = FatoresDentePermitidoContato.sigmaHP(HB)
	ZN = FatoresDentePermitidoContato.ZN()
	ZW = FatoresDentePermitidoContato.ZW()
	SH = FatoresDentePermitidoContato.SH()
	Ytheta = FatoresDentePermitidoContato.Ytheta()
	Yz = FatoresDentePermitidoContato.Yz()
	#print("[sigmaHP, ZN, Ytheta, Yz, ZW, SH] = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f]" % (sigmaHP, ZN, Ytheta, Yz, ZW, SF))
	sigma_max_contato = sigmaHP*ZN*ZW/(Ytheta*Yz)
	print("sigma_contato = [%.3f / %.3f] MPa" % (sigma_contato, sigma_max_contato))

	sf = (sigma_max_flexao/sigma_flexao)
	sh = (sigma_max_contato/sigma_contato)**2 
	print("Fator segurança:")
	print("SF = %.2f" % sf)
	print("SH = %.2f" % sh)
	print("\n")

	return sf, sh


if __name__ == "__main__":
	sfs, shs = np.zeros(4), np.zeros(4)
	bi, bf = 3, 5
	N = 0
	fator_seguranca = 2.84
	while N < 15:
		newb = (bi+bf)/2
		for i in range(4):
			sf, sh = calculc(i+1, newb)
			sfs[i] = sf
			shs[i] = sh
		minimo = min([min(sfs), min(shs)])
		if minimo < fator_seguranca:
			bi = newb
		else:
			bf = newb
		N += 1

	

	

	
	

	
	
	

