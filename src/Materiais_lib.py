# -*- coding: utf-8 -*-
"""

@author: carlos
"""

class Material:
	def __init__(self, nome = None):
		if nome != None:
			self.set(nome)
		else:
			self.__nome = None
			self.__HB	= None
			self.__E	= None
			self.__G	= None
			self.__v	= None
			self.__rho	= None
			self.__Sy	= None
			self.__Sfat	= None  
	def __del__(self):
		pass
	def set(self, nome):
		NOME = "../data/materiais.txt"
		arq = open(NOME, "r")
		linhas = arq.readlines()
		arq.close()
		lista_materiais = {}
		for i in range(len(linhas)):
			linhas[i] = linhas[i].strip().split(';')
			material = {}
			for j in range(1, len(linhas[i])):
				unidade, valor = linhas[i][j].split("=")
				material[unidade] = valor
			lista_materiais[linhas[i][0]] = material
		if nome in lista_materiais:
			material = lista_materiais[nome]
			self.__nome = nome
			self.__HB  = float(material["HB"])  # Brinell
			self.__E   = float(material["E"])   # Modulo de elasticidade, em MPa
			self.__G   = float(material["G"])   # Modulo de elasticidade torcional, em MPa
			self.__v   = float(material["v"])   # Coeficiente de poisson
			self.__rho = float(material["rho"]) # densidade, gramas/cm3
			self.__Sy  = float(material["Sy"])   # Coeficiente de poisson
			self.__Sfrat = float(material["Sfrat"]) # densidade, gramas/cm3
		else:
			raise ValueError("O material proposto nao esta cadastrado no sistema")
	def get_nome(self):
		if self.__nome == None:
			raise ValueError("Voce ainda nao definiu o material")
		return self.__nome
	def get_brinell(self):
		if self.__HB == None:
			raise ValueError("Voce ainda nao definiu o material")
		return self.__HB
	def get_elasticidade(self):
		if self.__E == None:
			raise ValueError("Voce ainda nao definiu o material")
		return self.__E
	def get_elasticidade_torcional(self):
		if self.__G == None:
			raise ValueError("Voce ainda nao definiu o material")
		return self.__G
	def get_limite_escoamento(self):
		if self.__Sy == None:
			raise ValueError("Voce ainda nao definiu o material")
		return self.__Sy
	def get_limite_resistencia_tracao(self):
		if self.__Sfrat == None:
			raise ValueError("Voce ainda nao definiu o material")
		return self.__Sfrat
	def get_poisson(self):
		if self.__v == None:
			raise ValueError("Voce ainda nao definiu o material")
		return self.__v
	def get_densidade(self):
		if self.__rho == None:
			raise ValueError("Voce ainda nao definiu o material")
		return self.__rho
	def get_infos(self):
		mensagem  = "Material " + self.get_nome() + "\n"
		mensagem += "Dureza Brinell: " + self.get_brinell() + " HB\n"
		mensagem += "Coef Elasticid: " + self.get_elasticidade() + "MPa\n"
		mensagem += "  Coef Poisson: " + self.get_poisson() + "\n"
		mensagem += "     Densidade: " + self.get_densidade() + "\n"

	def __str__(self):
		return self.get_nome()