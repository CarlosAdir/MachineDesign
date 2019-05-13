# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:08:04 2019

@author: carlos
"""

from Engrenagens_lib import Engrenagem
from Eixos_lib import Eixo

class Mate:
	@staticmethod
	def validar_engrenagem(eng):
		if not issubclass(type(eng), Engrenagem):
			raise TypeError("O argumento passado para Engrenagem não é Engrenagem!")
	@staticmethod
	def validar_eixo(eixo):
		if not issubclass(type(eixo), Eixo):
			raise TypeError("O argumento passado para eixo não é um eixo!")
	def __init__(self):
		pass

class Engrenamento(Mate):
	@staticmethod
	def __validar_razao_engrenagem(r):
		if type(r) != int and type(r) != float:
			raise TypeError("O valor passado como razão de engrenagem não é um numero!")
		if r < 1:
			raise ValueError("A razão de engrenagem deve ser maior ou igual a 1")
	@staticmethod
	def __validar_engrenamento(E1, E2):
		if E1 == E2:
			raise ValueError("Uma engrenagem não pode engrenar em si mesma!")
		if type(E1) != type(E2):
			raise ValueError("Uma engrenagem do tipo " + str(type(E1)) + " não pode engrenar em uma do tipo " + str(type(E2)))
	def __init__(self, Ep, Eg):
		Mate.__init__(self)
		self.set(Ep, Eg)
		self.__r = None
	def set(self, Ep, Eg):
		self.validar_engrenagem(Ep)
		self.validar_engrenagem(Eg)
		self.__validar_engrenamento(Ep, Eg)
		self.__Ep = Ep
		self.__Eg = Eg
	def get_pinhao(self):
		return self.__Ep
	def get_coroa(self):
		return self.__Eg
	def set_razao_engrenagem(self, r):
		self.__validar_razao_engrenagem(r)
		pinhao = self.get_pinhao()
		coroa = self.get_coroa()
		try:
			Np = pinhao.get_dentes()
		except:
			Np = None
		try:
			Ng = coroa.get_dentes()
		except:
			Ng = None
		if Np != None:
			coroa.set_dentes(Np*r)
		elif Np == None and Ng != None:
			newNp = Ng/r
			if int(newNp) != newNp:
				raise ValueError("A coroa não é multiplo da razão, gerando dente quebrado, errado!")
			pinhao.set_dentes(newNp)
		self.__r = r
	def get_razao_engrenagem(self):
		pinhao = self.get_pinhao()
		coroa = self.get_coroa()
		return self.__r
	def run(self):
		pass


class EngNoEixo(Mate):

	@staticmethod
	def __validar_posicao(p):
		if type(p) != int and type(p) != float:
			raise TypeError("A posição passada não é um valor! É necessário que seja")
		if not (0 <= p <= 1):
			raise ValueError("É necessário que a engrenagem esteja em uma posição entre 0 e 1")
	def __init__(self, engrenagem, eixo, posicao):
		Mate.__init__(self)
		self.set(engrenagem, eixo, posicao)
	def set(self, engrenagem, eixo, posicao):
		self.validar_engrenagem(engrenagem)
		self.validar_eixo(eixo)
		self.__validar_posicao(posicao)
		self.__engrenagem = engrenagem
		self.__eixo = eixo
		self.__posicao = posicao
		self.__engrenagem.set_eixo(eixo)
	def get_engrenagem(self):
		return self.__engrenagem
	def get_eixo(self):
		return self.__eixo
	def get_posicao(self):
		return self.__posicao
