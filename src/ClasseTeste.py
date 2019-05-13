from CalculaTensaoEngrenagens import Engrenagem, Material

class TU:
	SUCESSO = True
	FALHA = False
	def __init__(self, verbose = True):
		self._result = TU.SUCESSO # Ou seja, até o momento, tudo de acordo
		self._verbose = verbose
		self.__numerador = 1
		self._denominador = 1
	def print_result(self, mensagem):
		texto = "[" + str(self.__numerador) + "/" + str(self._denominador) + "] - "
		if self._result == TU.SUCESSO:
			texto += "\033[1;32m"
		else:
			texto += "\033[1;31m"
		texto += mensagem
		texto += "\033[0m"
		self.__numerador += 1
		print(texto)


class TUMaterial(TU):
	__MATERIAIS_VALIDOS 	= ["steel 1050"]
	__MATERIAIS_INVALIDOS 	= ["Steel 1050", "steel 1040", "aluminium"]

	def __init__(self, verbose = True):
		TU.__init__(self, verbose)
		testes = [	len(TUMaterial.__MATERIAIS_VALIDOS),	len(TUMaterial.__MATERIAIS_INVALIDOS) ]
		self._denominador = 1 + sum(testes) + 1
	def __setUp(self):
		try:
			self.material = Material()
		except:
			self._result = TU.FALHA
		if self._verbose:
			self.print_result("TUMaterial::setUp")
	def __MateriaisValidos(self):
		for nome_material in TUMaterial.__MATERIAIS_VALIDOS:
			try:
				self.material.set(nome_material)
				# if self.material.get() != material:
				# 	self.result = TU.FALHA
			except:
				self._result = TU.FALHA
			if self._verbose:
				self.print_result("TUMaterial::__MateriaisValidos")
			if self._result == TU.FALHA:
				break
	def __MateriaisInvalidos(self):
		for material in TUMaterial.__MATERIAIS_INVALIDOS:
			try:
				self.material.set(nome_material)
				self._result = TU.FALHA
			except:
				pass
			if self._verbose:
				self.print_result("TUMaterial::__MateriaisInvalidos")
			if self._result == TU.FALHA:
				break
	def __tearDown(self):
		try:
			del self.material
		except:
			self._result = TU.FALHA
		if self._verbose:
			self.print_result("TUMaterial::tearDown")
	def run(self):
		funcoes_teste = [	self.__setUp, \
							self.__MateriaisValidos, self.__MateriaisInvalidos, \
							self.__tearDown]
		for funcao in funcoes_teste:
			if self._result == TU.SUCESSO:
				funcao()
		return self._result

class TUEngrenagem(TU):

	__DENTES_VALIDOS		= [12, 13, 14, 15, 16, 17, 18, 19, 20]
	__DENTES_INVALIDOS		= ["a", -5, 0]
	__MODULOS_VALIDOS		= [1.25, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 30, 40, 50]
	__MODULOS_INVALIDOS		= ["w", -1, 0, 1.1, 3.2, 49, 51, 100]
	__ANGULOS_VALIDOS		= [15, 17.5, 20, 22.5, 25, 27.5, 29]
	__ANGULOS_INVALIDOS		= ["k", 0, 5, 10, 30, 33, 40, 60]
	__LARGURAS_VALIDAS		= [20, 30, 40, 50, 100]
	__LARGURAS_INVALIDAS	= ["KK", "12", -1, 0]
	

	def __init__(self, verbose = True):
		TU.__init__(self, verbose)
		testes = [	len(TUEngrenagem.__DENTES_VALIDOS),		len(TUEngrenagem.__DENTES_INVALIDOS), \
					len(TUEngrenagem.__MODULOS_VALIDOS),	len(TUEngrenagem.__MODULOS_INVALIDOS), \
					len(TUEngrenagem.__ANGULOS_VALIDOS),	len(TUEngrenagem.__ANGULOS_INVALIDOS), \
					len(TUEngrenagem.__LARGURAS_VALIDAS),	len(TUEngrenagem.__LARGURAS_INVALIDAS)	]
		self._denominador = 1 + sum(testes) + 1
	def __setUp(self):
		try:
			self.engrenagem = Engrenagem("a")
		except:
			self._result = TU.FALHA
		if self._verbose:
			self.print_result("TUEngrenagem::setUp")
	def __DentesValidos(self):
		for dentes in TUEngrenagem.__DENTES_VALIDOS:
			try:
				self.engrenagem.set_dentes(dentes)
				if self.engrenagem.get_dentes() != dentes:
					self.result = TU.FALHA
			except:
				self._result = TU.FALHA
			if self._verbose:
				self.print_result("TUEngrenagem::__DentesValidos")
			if self._result == TU.FALHA:
				break
	def __DentesInvalidos(self):
		for dentes in TUEngrenagem.__DENTES_INVALIDOS:
			try:
				self.engrenagem.set_dentes(dentes)
				self._result = TU.FALHA
			except:
				pass
			if self._verbose:
				self.print_result("TUEngrenagem::__DentesInvalidos")
			if self._result == TU.FALHA:
				break
	def __ModulosValidos(self):
		for modulo in TUEngrenagem.__MODULOS_VALIDOS:
			try:
				self.engrenagem.set_modulo(modulo)
				if self.engrenagem.get_modulo() != modulo:
					self.result = TU.FALHA
			except:
				self._result = TU.FALHA
			if self._verbose:
				self.print_result("TUEngrenagem::__ModulosValidos")
			if self._result == TU.FALHA:
				break
	def __ModulosInvalidos(self):
		for modulo in TUEngrenagem.__MODULOS_INVALIDOS:
			try:
				self.engrenagem.set_modulo(modulo)
				self._result = TU.FALHA
			except:
				pass
			if self._verbose:
				self.print_result("TUEngrenagem::__ModulosInvalidos")
			if self._result == TU.FALHA:
				break
	def __AngulosValidos(self):
		for angulo in TUEngrenagem.__ANGULOS_VALIDOS:
			try:
				self.engrenagem.set_angulo_pressao(angulo)
				if self.engrenagem.get_angulo_pressao_degrees() != angulo:
					self.result = TU.FALHA
			except:
				self._result = TU.FALHA
			if self._verbose:
				self.print_result("TUEngrenagem::__AngulosValidos")
			if self._result == TU.FALHA:
				break
	def __AngulosInvalidos(self):
		for angulo in TUEngrenagem.__ANGULOS_INVALIDOS:
			try:
				self.engrenagem.set_angulo_pressao(angulo)
				self._result = TU.FALHA
			except:
				pass
			if self._verbose:
				self.print_result("TUEngrenagem::__AngulosInvalidos")
			if self._result == TU.FALHA:
				break
	def __LargurasValidas(self):
		for largura in TUEngrenagem.__LARGURAS_VALIDAS:
			try:
				self.engrenagem.set_largura(largura)
				if self.engrenagem.get_largura() != largura:
					self.result = TU.FALHA
			except:
				self._result = TU.FALHA
			if self._verbose:
				self.print_result("TUEngrenagem::__LargurasValidas")
			if self._result == TU.FALHA:
				break
	def __LargurasInvalidas(self):
		for largura in TUEngrenagem.__LARGURAS_INVALIDAS:
			try:
				self.engrenagem.set_largura(largura)
				self._result = TU.FALHA
			except:
				pass
			if self._verbose:
				self.print_result("TUEngrenagem::__LargurasInvalidas")
			if self._result == TU.FALHA:
				break
	def __tearDown(self):
		try:
			del self.engrenagem
		except:
			self._result = TU.FALHA
		if self._verbose:
			self.print_result("TUEngrenagem::tearDown")
	def run(self):
		funcoes_teste = [	self.__setUp, \
							self.__DentesValidos, self.__DentesInvalidos, \
							self.__ModulosValidos, self.__ModulosInvalidos, \
							self.__AngulosValidos, self.__AngulosInvalidos, \
							self.__LargurasValidas, self.__LargurasInvalidas, \
							self.__tearDown]
		for funcao in funcoes_teste:
			if self._result == TU.SUCESSO:
				funcao()
		return self._result





if __name__ == "__main__":
	print("*** TESTE UNIDADE MATERIAL ***")
	tumaterial = TUMaterial()
	result = tumaterial.run()
	print(result)

	print("*** TESTE UNIDADE ENGRENAGEM ***")
	tuengrenagem = TUEngrenagem()
	result = tuengrenagem.run()
	print(result)


	pass
	