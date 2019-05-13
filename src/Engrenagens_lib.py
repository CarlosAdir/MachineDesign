# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:08:04 2019

@author: carlos
"""

import numpy as np
from Materiais_lib import *

class Engrenagem:
	@staticmethod
	def __validar_modulo(m):
		if type(m) != int and type(m) != float:
			raise TypeError("O modulo deve ser um numero!")
		if m <= 0:
			raise ValueError("O modulo deve ser positivo!")
	@staticmethod
	def __validar_dentes(N):
		if type(N) != int and type(N) != float:
			raise TypeError("O numero de dentes deve ser um numero!")
		if N <= 0:
			raise ValueError("O numero de dentes deve ser positivo!")
	@staticmethod
	def __validar_angulo_pressao(phi_n):
		if type(phi_n) != int and type(phi_n) != float:
			#print("phi_n = " + str(phi_n))
			raise TypeError("O angulo de pressao deve ser um numero!")
		if not phi_n > 0:
			raise TypeError("O angulo de pressao deve ser maior que 0")
	@staticmethod
	def __validar_largura(b):
		if type(b) != int and type(b) != float:
			raise TypeError("A largura da engrenagem deve ser um valor!")
		if not b > 0:
			raise TypeError("A largura da engrenagem deve ser maior que zero!")
	def __init__(self):
		self.__N = None
		self.__m = None
		self.__phi_n = None
		self.__material = None
		self.__b = None
		self.__mates = []
	def set_modulo(self, m):
		self.__validar_modulo(m)
		self.__m = m
	def get_modulo(self):
		if self.__m == None:
			raise ValueError("Nao foi possivel pegar o modulo porque ainda nao foi colocado")
		return self.__m
	def set_dentes(self, N):
		self.__validar_dentes(N)
		self.__N = N
	def get_dentes(self):
		if self.__N == None:
			raise ValueError("Nao foi possivel pegar o numero de dentes porque ainda nao foi colocado")
		return self.__N
	def get_diametro(self):
		try:
			m = self.get_modulo()
			N = self.get_dentes()
			return m*N
		except:
			raise ValueError("Para pegar o diametro é necessário ter o módulo e o numero de dentes")
	def set_material(self, mat):
		# Material
		self.__material = Material(mat)
	def get_material(self):
		if self.__material == None:
			raise ValueError("Não foi possivel pegar o material porque ainda não o colocaram!")
		return self.__material
	def set_angulo_pressao(self, phi_n):
		self.__validar_angulo_pressao(phi_n)
		self.__phi_n = phi_n
	def get_angulo_pressao_degrees(self):
		if self.__phi_n == None:
			raise ValueError("O valor do ângulo ainda não foi definido! Coloque!")
		return self.__phi_n
	def get_angulo_pressao_radians(self):
		# 
		return self.get_angulo_pressao_degrees() * (np.pi/180)
	def set_largura(self, b):
		self.__validar_largura(b)
		self.__b = b
	def get_largura(self):
		if self.__b == None:
			raise ValueError("Nao da pra pegar a largura se voce nao colocou a largura!")
		return self.__b

class EngInterna(Engrenagem):
	def __init__(self):
		Engrenagem.__init__(self)
	def get_massa(self):
		# Assumiremos que a massa da anular é um anel com diametro como o do circulo primitivo
		# E o diametro externo sera meia altura(metade do dente), somado com 1.2* altura do dente
		e = self.get_largura() # Dado em mm
		h = 2.166*self.get_modulo() # Dado em mm
		dint = self.get_diametro() # Dado emm mm
		dext = dint + (h/2) + 1.2*h
		Area = np.pi*(dext**2 - dint**2)/4
		rho = self.get_material().get_densidade() # Dado em g/cm3
		rho /= 1e+3 # transformamos a densidade em kg/cm3
		rho /= 1e+3 # transformamos a densidade em kg/mm3
		return e*Area*rho

	def __str__(self):
		mensagem  = "Numero de dentes: N     = %d \n" % self.get_dentes()
		mensagem += " Módulo do dente: m     = %d mm\n" % self.get_modulo()
		mensagem += "Diâmetro da engr: d     = %d mm\n" % self.get_diametro()
		mensagem += "        Material: mat   = " + str(self.get_material()) + "\n"
		mensagem += "  Angulo Pressao: phi_n = %.1f o\n" % self.get_angulo_pressao_degrees()
		mensagem += "         Largura: b     = %.1f mm\n" % self.get_largura()
		return mensagem

class EngIntDentesRetos(EngInterna):
	def __init__(self):
		EngInterna.__init__(self)
	def __str__(self):
		mensagem  = "Engrenagem Interna de dentes retos\n"
		mensagem += "\n" 
		mensagem += super().__str__()
		return mensagem

class EngIntDentesHelicoidais(EngInterna):
	def __init__(self):
		EngInterna.__init__(self)
	def __str__(self):
		mensagem  = "Engrenagem Interna de dentes helicoidais\n"
		mensagem += "\n" 
		mensagem += super().__str__()
		return mensagem


class EngExterna(Engrenagem):
	'''
	Engrenagem de dentes externos, sendo dentes retos ou helicoidais
	'''
	
	@staticmethod
	def __validar_razao_transmissao(r):
		if type(r) != int and type(r) != float:
			raise TypeError("A razão de transmissão deve ser um numero!")
		if r < 1:
			raise TypeError("A razão de transmissão não pode ser menor que 1")
	@staticmethod
	def menor_pinhao(r, phi_n = 20, gamma = 30):
		'''
		Calcula o mínimo número de dentes do pinhão para que não haja interferência.
		Esse mínimo pega como parâmetros:
		* r, que é a razão de transmissão.
		* phi_n é o ângulo de de pressão da engrenagem, tem como padrão 20° porque é bem comum
		* gamma é o angulo de hélice da engrenagem helicoidal
		'''
		EngExterna.__validar_razao_transmissao(r)
		gamma *= (np.pi/180)	# Transformar para radianos 
		phi_n *= (np.pi/180)	# Transformar para radianos
		if gamma == 0:
			phi_t = phi_n
		else:
			phi_t = np.arctan(np.tan(phi_n)/np.cos(gamma)) # Equacao 13-19 do shigley
		SIN = np.sin(phi_t)
		COS = np.cos(gamma)
		k   = 1				  # Eu nao sei bem de onde vem esse k, mas na equação (13-11) tem, e é 1 para full-depth teeth e 0.8 para stub teeth
		if r > 30:			  # Se caso for considerado uma cremalheira. O valor de 30(numero magico) é porque dificilmente se tem uma redução/aumento de 1 para 30
			Np  = 2*k/(SIN**2) # Equação 13-13 do shigley
		else:
			A   = 2*k*COS*(r + np.sqrt(r**2+(1+2*r)*SIN**2))
			B   = (1+2*r)*SIN**2
			Np  = A/B # Equação 13-11 do shigley
		return int(np.ceil(Np))
	@staticmethod
	def maior_coroa(Np, phi_n = 20, gamma = 30):
		'''
		Calcula o maior número de dentes da coroa para que não haja interferência.
		Esse mínimo pega como parâmetros:
		* Np, que é o numero de dentes do pinhão.
		* phi_n é o ângulo de de pressão da engrenagem, tem como padrão 20° porque é bem comum
		* gamma é o angulo de hélice da engrenagem helicoidal
		'''
		# Aqui verificamos qual o menor pinha que roda em uma cremalheira, com mesmos angulos
		Np_cremalheira = EngrenagemExterna.menor_pinhao(1000, phi_n, gamma)
		if Np_cremalheira <= Np: # Comparamos se o pinhao Np tem mais dentes que Np_cremalheira, porque se tiver mais, com certeza roda
			return 10000
		gamma *= (np.pi/180)	# Transformar para radianos 
		phi_n *= (np.pi/180)	# Transformar para radianos
		if gamma == 0:
			phi_t = phi_n
		else:
			phi_t = np.arctan(np.tan(phi_n)/np.cos(gamma)) # Equacao 13-19 do shigley
		SIN = np.sin(phi_t)
		COS = np.cos(gamma)
		k   = 1
		A   = Np**2 * SIN**2 - 4*k**2*COS**2
		B   = 4*k*COS - 2*Np*SIN**2
		Ng  = int(np.floor(A/B))
		return Ng
	def __init__(self):
		Engrenagem.__init__(self)
		self.__eixo = None
	def set_eixo(self, eixo):
		# Aqui deveriamos testar o eixo, mas como passaremos sempre o eixo, esperamos não ter problema
		self.__eixo = eixo
	def get_eixo(self):
		if self.__eixo == None:
			raise ValueError("Ainda não foi definido o eixo! Não foi possivel pega-lo!")
		return self.__eixo
	def get_Wt(self):
		T = self.get_eixo().get_torque() # o torque, em N*m
		print("T = %.2f N*m" % T)
		d = self.get_diametro()/1000  # O diametro em metros, para que Wt saia em N
		Wt = 2*T/d
		print("Wt = %.2f N" % Wt)
		return Wt
	def get_velocidade(self):
		n = self.get_eixo().get_rotacao()
		w = 2*np.pi*n/60 # Rotacao em rad/s
		d = self.get_diametro()
		v = d/(2*w)
		return v
	def get_massa(self):
		e = self.get_largura() # Dado em mm
		d  = self.get_diametro() # Dado em mm
		Area = np.pi*d**2/4
		rho = self.get_material().get_densidade() # Dado em g/cm3
		rho /= 1e+3 # transformamos a densidade em kg/cm3
		rho /= 1e+3 # transformamos a densidade em kg/mm3
		return e*Area*rho
	def __str__(self):
		mensagem  = "Numero de dentes: N     = %d \n" % self.get_dentes()
		mensagem += " Módulo do dente: m     = %d mm\n" % self.get_modulo()
		mensagem += "Diâmetro da engr: d     = %d mm\n" % self.get_diametro()
		mensagem += "        Material: mat   = " + str(self.get_material().get_nome()) + "\n"
		mensagem += "  Angulo Pressao: phi_n = %.1f o\n" % self.get_angulo_pressao_degrees()
		mensagem += "         Largura: b     = %.1f mm\n" % self.get_largura()
		mensagem += "           Massa: mass  = %.1f kg\n" % self.get_massa()
		return mensagem

class EngExtDentesRetos(EngExterna):
	@staticmethod
	def menor_pinhao(r, phi_n = 20):
		return EngExterna.menor_pinhao(r, phi_n, 0)
	@staticmethod
	def maior_coroa(Np, phi_n = 20):
		return EngExterna.maior_coroa(Np, phi_n, 0)
	def __init__(self):
		EngExterna.__init__(self)
	def __str__(self):
		mensagem  = "Engrenagem de dentes retos\n"
		mensagem += "\n" 
		mensagem += super().__str__()
		return mensagem

class EngExtDentesHelicoidais(EngExterna):
	@staticmethod
	def __validar_angulo_helice(gamma):
		if type(gamma) != int and type(gamma) != float:
			raise TypeError("O angulo de helice deve ser um numero")
		if not gamma > 0:
			raise ValueError("O angulo de helice deve ser maior que zero!")
	def __init__(self):
		EngExterna.__init__(self)
		self.__gamma = None
	def set_angulo_helice(self, gamma):
		self.__validar_angulo_helice(gamma)
		self.__gamma = gamma
	def get_angulo_helice_degrees(self):
		if self.__gamma == None:
			raise ValueError("Voce precisa colocar o angulo de helice primeiro!")
		return self.__gamma
	def get_angulo_helice_radians(self):
		return self.get_angulo_helice_degrees() * (np.pi/180)
	def __str__(self):
		mensagem  = "Engrenagem de dentes retos\n"
		mensagem += "\n" 
		mensagem += super().__str__()
		mensagem += "   Angulo Helice: gamma = " + str(self.get_angulo_helice_degrees())
		return mensagem


if __name__ == "__main__":
	E1 = EngExtDentesRetos()
	E1.set_modulo(1)
	E1.set_dentes(10)
	E1.set_angulo_pressao(20)
	E1.set_material("steel 1050")
	print(E1)

	print("\n\n\n")

	E2 = EngExtDentesHelicoidais()
	E2.set_modulo(1)
	E2.set_dentes(10)
	E2.set_angulo_pressao(20)
	E2.set_angulo_helice(30)
	E2.set_material("steel 1050")
	print(E2)