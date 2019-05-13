import Engrenagens_lib as Engs
import numpy as np
import Tensoes_dente as Tensoes

class TEPsimples:
	def __get_numero_planetas_possivel(self):
		Na = self.__anular.get_dentes()
		Ns = self.__solar.get_dentes()
		Np = self.__planeta.get_dentes()

		k = 1/(1+Ns/Np)
		Nmax = np.pi/np.arcsin(k)
		Nmax = int(np.floor(Nmax))

		Npossiveis = []
		Ntotal = Ns + Na

		for N in range(1, Nmax+1):
			divisao = Ntotal/N
			if divisao == int(divisao):
				Npossiveis.append(N)
		return Npossiveis
	def __init__(self, tipo):
		if tipo == 'DentesRetos':
			self.set_dentesretos()
		elif tipo == 'DentesHelicoidais':
			raise ValueError("Por enquanto ainda nao foi implementado para dentes helicoidais")
			self.set_helicoidais()
		else:
			raise TypeError('Os argumentos validos sao "DentesRetos" e "DentesHelicoidais"')
		self.__n_planetas = None
		self.__r = None
	def set_helicoidais(self):
		self.__anular = Engs.EngIntDentesHelicoidais()
		self.__planeta = Engs.EngExtDentesHelicoidais()
		self.__solar = Engs.EngExtDentesHelicoidais()
	def set_dentesretos(self):
		self.__anular = Engs.EngIntDentesRetos()
		self.__planeta = Engs.EngExtDentesRetos()
		self.__solar = Engs.EngExtDentesRetos()
	def set_modulo(self, m):
		self.__anular.set_modulo(m)
		self.__planeta.set_modulo(m)
		self.__solar.set_modulo(m)
	def set_angulo_pressao(self, phi_n):
		self.__anular.set_angulo_pressao(phi_n)
		self.__planeta.set_angulo_pressao(phi_n)
		self.__solar.set_angulo_pressao(phi_n)
	def set_largura(self, b):
		self.__anular.set_largura(b)
		self.__planeta.set_largura(b)
		self.__solar.set_largura(b)
	def set_material(self, M):
		self.__anular.set_material(M)
		self.__planeta.set_material(M)
		self.__solar.set_material(M)
	def get_modulo(self):
		return self.__anular.get_modulo()
	def get_largura(self):
		return self.__anular.get_largura()
	def get_angulo_pressao_degrees(self):
		return self.__anular.get_angulo_pressao_degrees()
	def get_angulo_pressao_radians(self):
		return self.__anular.get_angulo_pressao_radians()
	def get_dentes_anular(self):
		return self.__anular.get_dentes()
	def get_dentes_planeta(self):
		return self.__planeta.get_dentes()
	def get_dentes_solar(self):
		return self.__solar.get_dentes()
	def get_diametro_anular(self):
		return self.__anular.get_diametro()
	def get_diametro_planeta(self):
		return self.__planeta.get_diametro()
	def get_diametro_solar(self):
		return self.__solar.get_diametro()
	def get_massa_anular(self):
		return self.__anular.get_massa()
	def get_massa_planeta(self):
		return self.__planeta.get_massa()
	def get_massa_solar(self):
		return self.__solar.get_massa()
	def set_razao_transmissao(self, r, Ns):
		# e = -Na/Ns
		# Como Ns < Na e ambos positivos, entao e < -1
		#
		# Em um TEP simples, temos
		# Na = Ns + 2*Np
		#
		# r é a razao de transmissao
		# r é a razão entre a rotaçao entre a solar e o braço que as planetas movem
		if r == 1:
			raise ValueError("Pra que raios voce quer uma razao de 1?")
		elif r == 0:
			raise ValueError("Nao existe possibilidade do valor ser zero!")
		# Nesse caso estamos supondo que a anular e fixa
		if r > 1:
			e = 1-r
		if r < 1:
			e = 1-(1/r)
		Na = -e*Ns
		if int(Na) != Na:
			raise ValueError("Nao é possivel colocar o numero de dentes como %.1f" % Na)
		Np = (Na-Ns)/2
		if int(Np) != Np:
			raise ValueError("Nao é possivel colocar o numero de dentes como %.1f" % Np)
		self.__r = r
		self.__anular.set_dentes(Na)
		self.__planeta.set_dentes(Np)
		self.__solar.set_dentes(Ns)
	def get_razao_transmissao(self):
		if self.__r == None:
			raise ValueError("Voce ainda nao colocou a razao de transmissao!")
		return self.__r
	def set_numero_planetas(self, N):
		possiveis = self.__get_numero_planetas_possivel()
		if not N in possiveis:
			mensagem  = "\nNão é possivel ter " + str(N) + " planetas no TEP.\n" 
			mensagem += "Sao possiveis apenas os valores de " + str(possiveis)
			raise ValueError(mensagem)
		self.__n_planetas = N
	def set_numero_planetas_maximo(self):
		possiveis = self.__get_numero_planetas_possivel()
		self.__n_planetas = max(possiveis)
	def get_numero_planetas(self):
		if self.__n_planetas == None:
			raise ValueError("Voce ainda nao colocou o numero de planetas!")
		return self.__n_planetas
	def get_material(self):
		return self.__planeta.get_material()
	def __str__(self):
		debug = True
		Na, Np, Ns = self.get_dentes_anular(), self.get_dentes_planeta(), self.get_dentes_solar()
		da, dp, ds = self.get_diametro_anular(), self.get_diametro_planeta(), self.get_diametro_solar()
		ma, mp, ms = self.get_massa_anular(), self.get_massa_planeta(), self.get_massa_solar()
		mensagem  = "Trem epicicloidal\n"
		mensagem += "\n"
		mensagem += (" Modulo: m     = %.1f mm\n" % self.get_modulo())
		mensagem += ("Largura: b     = %.1f mm\n" % self.get_largura())
		if debug:
			mensagem += ("Larg. perm. = [%.1f, %.1f]\n" % (3*np.pi*self.get_modulo(), 5*np.pi*self.get_modulo()))
		mensagem += (" Angulo: phi_n = %.1f graus\n" % self.get_angulo_pressao_degrees())
		mensagem += ("Planeta: N     = %d planetas\n" % self.get_numero_planetas())
		mensagem += ("Razao T: r     = %.1f\n" % self.get_razao_transmissao())
		mensagem += "\n"
		mensagem +=  "          Anular | Planeta | Solar | \n"
		mensagem += ("  Dentes:  %3d   |   %3d   |  %3d  | \n" % (Na, Np, Ns) )
		mensagem += ("Diametro:  %4d  |  %4d   | %4d  | mm\n" % (da, dp, ds) )
		mensagem += ("   Massa:  %3d   |   %3d   |  %3d  | kg\n" % (np.ceil(ma), np.ceil(mp), np.ceil(ms)) )

		return mensagem


class TEPsimples_modificado(TEPsimples):
	@staticmethod
	def __validar_rotacao(n):
		if type(n) != int and type(n) != float:
			raise TypeError("A rotacao deve ser um numero em rpm!")
	@staticmethod
	def __validar_potencia(P):
		if type(P) != int and type(P) != float:
			raise TypeError("A rotacao deve ser um numero em kW!")	
	@staticmethod
	def __validar_eficiencia(efi):
		if type(efi) != int and type(efi) != float:
			raise TypeError("A eficiencia deve ser um numero no intervalo (0, 1]")
		if not (0 < efi or efi <= 1):
			raise ValueError("A eficiencia deve ser um numero no intervalo (0, 1]")
	@staticmethod
	def __validar_massa(m):
		if type(m) != int and type(m) != float:
			raise TypeError("A massa deve ser um numero maior ou igual a zero!")
		if m < 0:
			raise ValueError("A massa deve ser um numero maior ou igual a zero!")
	def __init__(self, tipo):
		TEPsimples.__init__(self, tipo)
		self.__n_in = None
		self.__P = None
		self.set_eficiencia(1)
		self.set_massa_braco(0)
	def set_massa_braco(self, m):
		self.__validar_massa(m)
		self.__massa_braco = m
	def get_massa_braco(self):
		return self.__massa_braco
	def set_eficiencia(self, efi):
		self.__validar_eficiencia(efi)
		self.__efi = efi
	def get_eficiencia(self):
		return self.__efi
	def set_rotacao_braco(self, n):
		self.__validar_rotacao(n)
		self.__n_in = n
	def get_rotacao_braco(self):
		if self.__n_in == None:
			raise ValueError("Voce ainda nao colocou o valor da rotacao de entrada")
		return self.__n_in
	def get_rotacao_solar(self):
		#
		return self.get_razao_transmissao()*self.get_rotacao_braco()
	def get_rotacao_planeta(self):
		Ns = self.get_dentes_solar()
		Np = self.get_dentes_planeta()
		nb = self.get_rotacao_braco()
		return (1+Ns/Np)*nb
	def get_rotacao_anular(self):
		#
		return 0 # Porque nossa hiptese é que a anular esteja parada


	def set_potencia_braco(self, P):
		self.__validar_potencia(P)
		self.__P = P
	def get_potencia_braco(self):
		if self.__P == None:
			raise ValueError("Voce ainda nao colocou a potencia no braço!")
		return self.__P
	def get_torque_braco(self):
		P = self.get_potencia_braco() * 1e+3 # Para transformar o valor da potencia de kW para W
		n = self.get_rotacao_braco() # Retornado em rpm
		w = 2*np.pi*n/60
		return P/w # Vai ser retornado em N*m
	def get_torque_solar(self):
		efi = self.get_eficiencia()
		torque_braco = self.get_torque_braco()
		torque_solar = torque_braco*efi/self.get_razao_transmissao()
		return torque_solar # Vai ser retornado em N*m
	def get_potencia_solar(self):
		torque_solar = self.get_torque_solar()
		n = self.get_rotacao_solar()
		w = 2*np.pi*n/60
		P = torque_solar * w
		return P/1e+3 # Iria retornar em W, mas a potencia colocaremos como kW
	def get_forca_na_planeta(self):
		gravidade = 9.85

		torque_braco = self.get_torque_braco()
		m = self.get_modulo()
		Ns = self.get_dentes_solar()
		Np = self.get_dentes_planeta()
		distancia_centro = m*(Ns+Np)/2
		Forca_motora = 1000*torque_braco/(distancia_centro * self.get_numero_planetas()) # A distancia esta em mm, para isso que multiplicamos o 1000 no inicio, para poder passar para metros e a forca voltar em N
		Forca_peso   = gravidade*self.get_massa_planeta()
		return Forca_peso + Forca_motora
	def get_Ft(self):
		forca	= self.get_forca_na_planeta()
		return forca/2
	def get_Fn(self):
		Ft		= self.get_Ft()
		phi_n	= self.get_angulo_pressao_radians()
		return Ft*np.tan(phi_n)
	def get_velocidade(self):
		# Assumiremos que sera a velocidade da planeta
		d = self.get_diametro_planeta() / 1000 # Para passar para metros
		n_p = self.get_rotacao_planeta()
		w = 2*np.pi*n_p/60 # Para passar de rpm para rad/s
		return w*d/2
	
	def get_fatores(self, verbose = False):
		b = self.get_largura()
		m = self.get_modulo()
		v = self.get_velocidade()
		mat = self.get_material()
		HB = mat.get_brinell()
		phi_n = self.get_angulo_pressao_degrees()

		Ns = self.get_dentes_solar()
		Np = self.get_dentes_planeta()
		Na = self.get_dentes_anular()
		ds = self.get_diametro_solar()
		dp = self.get_diametro_planeta()
		da = self.get_diametro_anular()

		Ft = self.get_Ft()

		sfs, shs = np.zeros(4), np.zeros(4)
		for i in range(4):
			if i == 0:
				mensagem = "Calculo dos fatores da planeta engrenando com a solar:"
				Npar, Nimpar = Ns, Np
				dpar, dimpar = ds, dp
			elif i == 1:
				mensagem = "Calculo dos fatores da planeta engrenando com a anular:"
				Npar, Nimpar = Na, Np
				dpar, dimpar = da, dp
			elif i == 2:
				mensagem = "Calculo dos fatores da solar engrenando com a planeta:"
				Npar, Nimpar = Np, Ns
				dpar, dimpar = dp, ds
			elif i == 3:
				mensagem = "Calculo dos fatores da anular engrenando com a planeta:"
				Npar, Nimpar = Np, Na
				dpar, dimpar = dp, da
			if verbose:
				print(mensagem)

			dpinhao = dpar if dpar < dimpar else dimpar
			dw1 = dpinhao

			#print("Flexao")
			K0 			= Tensoes.FatoresDenteFundamentalFlexao.K0()
			Kv 			= Tensoes.FatoresDenteFundamentalFlexao.Kv(Qv = 11, v = v)
			Ks 			= Tensoes.FatoresDenteFundamentalFlexao.Ks(N = Nimpar, b = b, m = m)
			KH 			= Tensoes.FatoresDenteFundamentalFlexao.KH(b = b, d = dimpar)
			mt 			= Tensoes.FatoresDenteFundamentalFlexao.mt(m, 0) # é necessario que esteja em mm, o zero é o angulo gamma
			KB 			= Tensoes.FatoresDenteFundamentalFlexao.KB()
			YJ 			= Tensoes.FatoresDenteFundamentalFlexao.YJ(Npar, Nimpar)
			#print("[K0, Kv, Ks, KH, mt, KB, YJ] = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]" % (K0, Kv, Ks, KH, mt, KB, YJ))
			sigma_flexao = Ft*K0*Kv*Ks*KH*KB/(b*mt*YJ)
			
			sigmaFP 	= Tensoes.FatoresDentePermitidoFlexao.sigmaFP(HB)
			YN 			= Tensoes.FatoresDentePermitidoFlexao.YN()
			Ytheta 		= Tensoes.FatoresDentePermitidoFlexao.Ytheta()
			Yz 			= Tensoes.FatoresDentePermitidoFlexao.Yz()
			SF 			= Tensoes.FatoresDentePermitidoFlexao.SF()
			#print("[sigmaFP, YN, Ytheta, Yz, SF] = [%.2f, %.2f, %.2f, %.2f, %.2f]" % (sigmaFP, YN, Ytheta, Yz, SF))
			sigma_max_flexao = (sigmaFP*YN)/(Ytheta*Yz)  #########

			if verbose:
				print("sigma_flexao = [%.3f / %.3f] MPa" % (sigma_flexao, sigma_max_flexao))

			#print("Contato")
			K0 			= Tensoes.FatoresDenteFundamentalContato.K0()
			Kv 			= Tensoes.FatoresDenteFundamentalContato.Kv(Qv = 11, v = v)
			Ks 			= Tensoes.FatoresDenteFundamentalContato.Ks(N = Nimpar, b = b, m = m)
			KH 			= Tensoes.FatoresDenteFundamentalContato.KH(b = b, d = dimpar)
			ZI 			= Tensoes.FatoresDenteFundamentalContato.ZI(phi_n, 2)
			ZR 			= Tensoes.FatoresDenteFundamentalContato.ZR()
			ZE 			= Tensoes.FatoresDenteFundamentalContato.ZE(mat, mat)
			#print("[K0, Kv, Ks, KH, ZE, ZR, ZI] = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]" % (K0, Kv, Ks, KH, ZE, ZR, ZI))
			sigma_contato = ZE*np.sqrt(Ft*K0*Kv*Ks*KH*ZR/(dw1*ZI*b))

			sigmaHP 	= Tensoes.FatoresDentePermitidoContato.sigmaHP(HB)
			ZN 			= Tensoes.FatoresDentePermitidoContato.ZN()
			ZW 			= Tensoes.FatoresDentePermitidoContato.ZW()
			SH 			= Tensoes.FatoresDentePermitidoContato.SH()
			Ytheta 		= Tensoes.FatoresDentePermitidoContato.Ytheta()
			Yz 			= Tensoes.FatoresDentePermitidoContato.Yz()
			#print("[sigmaHP, ZN, Ytheta, Yz, ZW, SH] = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f]" % (sigmaHP, ZN, Ytheta, Yz, ZW, SF))
			sigma_max_contato = sigmaHP*ZN*ZW/(Ytheta*Yz)
			if verbose:
				print("sigma_contato = [%.3f / %.3f] MPa" % (sigma_contato, sigma_max_contato))

			sfs[i] = (sigma_max_flexao/sigma_flexao)
			shs[i] = (sigma_max_contato/sigma_contato)**2 
	
		return min(sfs), min(shs)
	def msg1(self):
		n_b = self.get_rotacao_braco()
		n_s = self.get_rotacao_solar()
		n_p = self.get_rotacao_planeta()
		n_a = self.get_rotacao_anular()

		w_b = n_b * 2*(np.pi)/60
		w_s = n_s * 2*(np.pi)/60
		w_p = n_p * 2*(np.pi)/60
		w_a = n_a * 2*(np.pi)/60

		da = self.get_diametro_anular()
		dp = self.get_diametro_planeta()
		ds = self.get_diametro_solar()
		db = dp+ds
		
		mensagem  = ("Rotacao   Braço: n_b = %.2f rpm\n" % n_b )
		mensagem += ("Rotacao   Solar: n_s = %.2f rpm\n" % n_s )
		mensagem += ("Rotacao Planeta: n_p = %.2f rpm\n" % n_p )
		mensagem += ("Rotacao  Anular: n_a = %.2f rpm\n" % n_a )
		mensagem += "\n"
		mensagem += ("Vel. cent. plan: v_c = %.2f m/s\n" % (w_b*db/2000) ) # Para passar o diametro para metro e o dois é porque precisamos do raio, nao do diametro
		mensagem += ("Vel. tang. plan: v_p = %.2f m/s\n" % (w_p*dp/2000) )
		mensagem += ("Vel. tang. sol.: v_s = %.2f m/s\n" % (w_s*ds/2000) )
		mensagem += ("Vel. tang. anu.: v_a = %.2f m/s\n" % (w_a*da/2000) )
		mensagem += "\n"
		mensagem += ("Potencia no braço: P_b = %.2f kW\n" % self.get_potencia_braco() )
		mensagem += ("  Torque no braço: T_b = %.2f N*m\n" % self.get_torque_braco() )
		mensagem += ("Potencia na solar: P_s = %.2f kW\n" % self.get_potencia_solar() )
		mensagem += ("  Torque na solar: T_s = %.2f N*m\n" % self.get_torque_solar() )
		return mensagem
	def msg2(self):
		mensagem  = ("Forca do braco na planeta: Fp = %.2f N\n" % self.get_forca_na_planeta() )
		mensagem += (" Forca tangencial das eng: Ft = %.2f N\n" % self.get_Ft())
		mensagem += ("     Forca radial das eng: Fn = %.2f N\n" % self.get_Fn())
		
		return mensagem




if __name__ == "__main__":
	pass

