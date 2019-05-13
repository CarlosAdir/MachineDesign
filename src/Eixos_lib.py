# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:08:04 2019

@author: carlos
"""

import numpy as np
from Materiais_lib import Material

class Eixo:
	@staticmethod
	def __validar_rotacao(n):
		if type(n) != int and type(n) != float:
			raise TypeError("A rotação deve ser um número!")
		if n < 0:
			raise ValueError("A rotacao deve ser positiva!")
	@staticmethod
	def __validar_potencia(P):
		if type(P) != int and type(P) != float:
			raise TypeError("A potencia deve ser um número!")
		if P < 0:
			raise ValueError("A potencia deve ser positiva!")
	@staticmethod
	def __validar_diametro(d):
		if type(d) != int and type(d) != float:
			raise TypeError("O diametro deve ser um número!")
		if not 0 < d:
			raise ValueError("O diametro deve ser positivo!")
	def __init__(self):
		self.__n = None
		self.__P = None
		self.__L = None
		self.__d = None
		self.__material = None
	def set_diametro(self, d):
		self.__validar_diametro(d)
		self.__d = d
	def set_material(self, mat):
		self.__material = Material(mat)

	def set_rotacao(self, n):
		self.__validar_rotacao(n)
		self.__n = n
	def set_potencia(self, P):
		self.__validar_potencia(P)
		self.__P = P
	def set_comprimento(self, L):
		self.__validar_comprimento(L)
		self.__L = L
	
	def get_material(self):
		return self.__material
	def get_diametro(self):
		if self.__d == None:
			raise ValueError("O diametro ainda nao foi colocado!")
		return self.__d

	def get_peso_por_comprimento(self):
		material = self.get_material()
		d 		 = self.get_diametro()
		rho 	 = material.get_densidade()
		rho 	/= 1e+6
		return (np.pi)*(d**2)*rho/4
	def get_rotacao(self):
		if self.__n == None:
			raise ValueError("Ainda nao foi colocado a rotacao")
		return self.__n
	def get_potencia(self):
		if self.__P == None:
			raise ValueError("Ainda nao foi colocado a potência")
		return self.__P
	def get_comprimento(self):
		if self.__L == None:
			raise ValueError("Ainda nao foi colocado o comprimento do eixo")
		return self.__L
	def get_torque(self):
		# Aqui consideraremos que P esta em kW
		# E que n esteja em rpm
		P = self.get_potencia()
		n = self.get_rotacao()
		w = 2*np.pi*n/60 # Rotacao em rad/s
		T = 1000*P/w # Retornaria em kN * m, então colocando o fator de 1000, o retorno fica na verdade N*m
		return T

	def get_momento_inercia_y(self):
		#t = sp.symbols('t')
		#fd = 0
		#for lista in self.listas:
		#	a, b, d = lista
		#	fd += d*(delta(t-a) - delta(t-b))
		#I = np.pi*(fd**4)/64
		#I = sp.lambdify(t, I, modules = ["numpy"])
		I = lambda t: np.pi*(self.get_diametro()**4)/64
		return I
	def get_momento_inercia_z(self):
		# Porque temos um eixo, que tem simetria radial
		return self.get_momento_inercia_y()
	def get_momento_polar(self):
		'''
		t = sp.symbols('t')
		fd = 0
		for lista in self.listas:
			a, b, d = lista
			fd += d*(delta(t-a) - delta(t-b))
		J = np.pi*(fd**4)/32
		J = sp.lambdify(t, J, modules = ["numpy"])
		'''
		J = lambda t: np.pi*(self.get_diametro()**4)/32
		return J