
from Engrenagens_lib import EngExtDentesRetos, EngExtDentesHelicoidais
from ConjuntoEngrenagem import TEPsimples_modificado
from Eixos_lib import Eixo
#from Mates_lib import EngNoEixo, Engrenamento
from Tensoes_dente import EsforcoNoDente, FatoresDentePermitido
from ProjetoDeEixos import Distribuicao, Forca, Momento, calcula_diagramas_esforcos, get_deslocamentos, plota_diagramas, plota_deslocamentos

import numpy as np
import sys



def settings(numero, eixo):
	distribuicoes	= []
	forcas			= []
	momentos		= []
	if numero == 1:
		# Eixo A
		L1, L2, L3, L4, L5 = 310, 100, 30, 100, 70
		P0		= 500 * gravidade # 500 kg de massa, vezes a gravidade, dá em Newton
		C0		= 300 # mm
		q1		= [0, -P0/C0, 0]
		q2		= [0, -eixo.get_peso_por_comprimento(), 0] # Peso em N por unidade de mm
		F1		= [0, -(tep1.get_numero_planetas()*tep1.get_massa_planeta() + tep1.get_massa_braco())*gravidade, 0]
		M1 		= [tep1.get_torque_braco()*1e+3, 0, 0] # Porque o torque retorna em N*m, e nao em N*mm
		M2		= [-tep1.get_torque_braco()*1e+3, 0, 0] # Porque o torque retorna em N*m, e nao em N*mm
		M3 		= [0, 0, L5*F1[1]]
		
		L		= L1 + L2 + L3 + L4 + L5
		a, b 	= L1, L1+L2 # A posicao dos mancais
		

		distribuicoes.append(Distribuicao(q1, [0, C0]))
		distribuicoes.append(Distribuicao(q2, [0, L]))
		forcas.append(Forca(F1, L))
		momentos.append(Momento(M1, 0))
		momentos.append(Momento(M2, L))
		momentos.append(Momento(M3, L))
	elif numero == 2:
		L6, L7, L8, L9, L10 = 50, 100, 30, 100, 70

		q1		= [0, -eixo.get_peso_por_comprimento(), 0]
		F1		= [0, -tep1.get_massa_solar()*gravidade, 0]
		F2		= [0, -(tep2.get_numero_planetas()*tep2.get_massa_planeta() + tep2.get_massa_braco())*gravidade, 0]
		M1		= [tep1.get_torque_solar()*1e+3, 0, 0] # Porque o torque retorna em N*m, e nao em N*mm
		M2		= [-tep2.get_torque_braco()*1e+3, 0, 0] # Porque o torque retorna em N*m, e nao em N*mm
		M3		= [0, 0, L10*F2[1]]

		a, b 	= L6, L6+L7
		L 		= L6 + L7 + L8 + L9 + L10

		distribuicoes.append(Distribuicao(q1, [0, L]))
		forcas.append(Forca(F1, 0))
		forcas.append(Forca(F2, L))
		momentos.append(Momento(M1, 0))
		momentos.append(Momento(M2, L))
		momentos.append(Momento(M3, L))
	elif numero == 3:
		L11, L12, L13 = 100, 100, 160
		P1		= 150 * gravidade # 150 kg de massa, vezes a gravidade, dá em Newton
		C1		= 150 # mm
		q1		= [0, -P1/C1, 0]
		q2		= [0, -eixo.get_peso_por_comprimento(), 0] # Peso em N por unidade de mm
		F1		= [0, -tep2.get_massa_solar()*gravidade, 0]
		M1 		= [tep2.get_torque_solar()*1e+3, 0, 0] # Porque o torque retorna em N*m, e nao em N*mm
		M2		= [-tep2.get_torque_solar()*1e+3, 0, 0] # Porque o torque retorna em N*m, e nao em N*mm
		
		L		= L11 + L12 + L13
		a, b 	= L11, L11+L12 # A posicao dos mancais

		distribuicoes.append(Distribuicao(q1, [L-C1, L]))
		distribuicoes.append(Distribuicao(q2, [0, L]))
		forcas.append(Forca(F1, L))
		momentos.append(Momento(M1, 0))
		momentos.append(Momento(M2, L))
	else:
		raise ValueError('Nao e possivel o argumento"' + str(numero) + '"')
	return [distribuicoes, forcas, momentos], [a, b, L]

if __name__ == "__main__":
	try:
		args = sys.argv
		teste = args[1]
		teste = int(teste)
		if teste < 0 or teste > 3:
			argumento_invalido = True
		else:
			argumento_invalido = False
	except:
		argumento_invalido = True

	if argumento_invalido:
		raise ValueError("Hey, coloque algo do tipo: python3 arquivo 3")
	
	log 		= open("resultados.txt", "w")
	tep1		= TEPsimples_modificado("DentesRetos")
	tep2		= TEPsimples_modificado("DentesRetos")
	

	eficiencia	= 0.975
	gravidade	= 9.85
	
	

	tep1.set_modulo(9)
	tep1.set_angulo_pressao(20) # angulo de pressao de 20 graus
	tep1.set_razao_transmissao(5, Ns = 28) # uma ampliacao de 5, solar com 20 dentes
	tep1.set_numero_planetas_maximo()
	tep1.set_largura(102)
	tep1.set_material("steel 1045")
	tep1.set_rotacao_braco(50) # rpm
	tep1.set_potencia_braco(420) # kW
	#tep1.set_eficiencia(eficiencia)
	#tep1.set_massa_braco(0)
	

	tep2.set_modulo(4)
	tep2.set_angulo_pressao(20) # angulo de pressao de 20 graus
	tep2.set_razao_transmissao(6, Ns = 32)
	tep2.set_numero_planetas_maximo()
	tep2.set_largura(54)
	tep2.set_material("steel 1045")
	tep2.set_rotacao_braco(tep1.get_rotacao_solar())
	tep2.set_potencia_braco(tep1.get_potencia_solar()) # kW
	#tep2.set_eficiencia(eficiencia)
	#tep2.set_massa_braco(0)

	if teste != 0:
		mancais = np.zeros(6)
		defl_max 	= 0.12
		inc_max 	= 0.0087
		fator_seg 	= 2.82
		EixoA		= Eixo()
		EixoB		= Eixo()
		EixoC		= Eixo()
		EixoA.set_diametro(100) # Apenas uma estimativa inicial
		EixoA.set_material("steel 1045")

		EixoB.set_diametro(100)
		EixoB.set_material("steel 1045")

		EixoC.set_diametro(80)
		EixoC.set_material("steel 1045")
		for i in range(1, 3+1):

			if i == 1:
				eixo = EixoA
			elif i == 2:
				eixo = EixoB
			elif i == 3:
				eixo = EixoC
			
			#log.write('eixo = ' + str(eixo))
			#log.write('type(eixo) = ' + str(type(eixo)))
			materi 	= eixo.get_material()
			E		= materi.get_elasticidade()
			G		= materi.get_elasticidade_torcional()
			Sy		= materi.get_limite_escoamento()
			Sfrat	= materi.get_limite_resistencia_tracao()

			Se 		= (Sfrat/2)


			esforcos, [a, b, L] = settings(i, eixo)
			[N, Vy, Vz], [T, My, Mz], [Fa, Fb] = calcula_diagramas_esforcos(esforcos, [a, b], [0, L])
			

			d_tor, d_est, d_fad, d_def = [0], [0], [0], [0]

			N 		= 1200
			X 		= np.linspace(0, L, N)
			Torqus 	= T(X)
			Momens 	= Mz(X)

			# Projeto para torcao
			if 1:
				#
				d_tor	= ( ((32*180*20)/(np.pi**2)) * (np.abs(Torqus)/G) )**(1/3)

			# Projeto para estatico
			if 1:
				module	= (np.abs(Momens**2 + Torqus**2))**(1/2)
				d_est	= (32*fator_seg*module/(np.pi*Sy))**(1/3)

			# Projeto para fadiga
			if 1:
				Kf  	= 1.0 # Fator para flexao
				Kfs 	= 1.0 # Fator para torcao
				Ma, Ta 	= 0, 0 # Fletor e Torsor alternado
				A 		= np.sqrt(4*(Kf*Ma)**2 		+ 3*(Kfs*Ta)**2)
				B 		= np.sqrt(4*(Kf*Momens)**2 	+ 3*(Kfs*Torqus)**2)
				d_fad 	= (16*fator_seg*np.sqrt(np.abs((A/Se)**2 + (B/Sy)**2))/np.pi)**(1/3)

			# Projeto para deflexao
			if 1:
				dmin = 10
				dmax = 1000
				iterador, iterador_max = 1, 15
				while iterador < iterador_max:
					d = (dmin+dmax)/2

					eixo.set_diametro(d)
					esforcos, [a, b, L] = settings(i, eixo)
					[N, Vy, Vz], [T, My, Mz], [Fa, Fb] = calcula_diagramas_esforcos(esforcos, [a, b], [0, L])
					
					Iz = eixo.get_momento_inercia_z()
					Iy = eixo.get_momento_inercia_y()
					J = eixo.get_momento_polar()
					y, theta 	= get_deslocamentos([My, Mz], E, [Iy, Iz], [0, L], [a, b])

					theta_max 	= max(np.abs(theta))
					y_max 		= max(np.abs(y))
					if (y_max < defl_max) and (theta_max < inc_max):
						dmax = d
					else:
						dmin = d
					iterador += 1
					print("i = %d/%d" % (iterador, iterador_max))
				d_def = [dmax]			

			
			diametro = float(max([max(d_def), max(d_tor), max(d_est), max(d_fad)]))
			#log.write('A, B, C = ' + str(A) + " " + str(B) + " " + str(C))
			if 1:

				if i == 1:
					log.write("Eixo A\n")
				elif i == 2:
					log.write("Eixo B\n")
				elif i == 3:
					log.write("Eixo C\n")
				log.write('max def = %.2f mm\n' % max(d_def))
				log.write('max tor = %.2f mm\n' % max(d_tor))
				log.write('max est = %.2f mm\n' % max(d_est))
				log.write('max fad = %.2f mm\n' % max(d_fad))
				log.write('Minimo diametro: d = %.2f mm\n' % diametro)
				log.write('\n\n')
				mancais[2*i-2] = np.linalg.norm(Fa)
				mancais[2*i-1] = np.linalg.norm(Fb)

			if 1:
				eixo.set_diametro(diametro)
				esforcos, [a, b, L] = settings(i, eixo)
				[N, Vy, Vz], [T, My, Mz], [Fa, Fb] = calcula_diagramas_esforcos(esforcos, [a, b], [0, L])
				
				Iz			= eixo.get_momento_inercia_z()
				Iy			= eixo.get_momento_inercia_y()
				J			= eixo.get_momento_polar()
				y, theta 	= get_deslocamentos([My, Mz], E, [Iy, Iz], [0, L], [a, b])
				plota_diagramas([N, Vy, Vz], [T, My, Mz], [0, L], nome = str(i))
				plota_deslocamentos([0, L], y, theta, nome = str(i))
		for j in range(6):
			log.write("Mancal %d: %.2f N\n" % (j+1, mancais[j]))
		log.write('\n\n\n\n')
	if 1:
		log.write(str(tep1) + '\n')
		log.write(str(tep1.msg1()) + '\n')
		log.write(str(tep1.msg2()) + '\n')
		log.write("[SF, SH] = [%.2f, %.2f]\n" % tep1.get_fatores())

	if 1:
		log.write("\n\n\n")

		log.write(str(tep2) + '\n')
		log.write(str(tep2.msg1()) + '\n')
		log.write(str(tep2.msg2()) + '\n')
		log.write("[SF, SH] = [%.2f, %.2f]\n" % tep2.get_fatores())
	log.close()

	


'''
if __name__ == "__main__":

	m1, m2 = 10, 8
	r1, r2 = 2, 3
	np1, np2 = 30, 20
	P = 420
	b1, b2 = 126, 106# 4*np.pi*m1, 4.2*np.pi*m2

	A, B, C = Eixo(), Eixo(), Eixo()
	E1 = EngExtDentesRetos()
	E2 = EngExtDentesRetos()
	E3 = EngExtDentesRetos()
	E4 = EngExtDentesRetos()

	K1 = EngNoEixo(E1, A, 0.7)
	K2 = EngNoEixo(E2, B, 0.7)
	K3 = EngNoEixo(E3, B, 0.3)
	K4 = EngNoEixo(E4, C, 0.3)

	# M1 = Engrenamento(Eg = E1, Ep = E2)
	# M2 = Engrenamento(Eg = E3, Ep = E4)

	E1.set_modulo(m1)
	E2.set_modulo(m1)
	E3.set_modulo(m2)
	E4.set_modulo(m2)


	E1.set_dentes(r1*np1)
	E2.set_dentes(np1)
	E3.set_dentes(r2*np2)
	E4.set_dentes(np2)
	E1.set_angulo_pressao(20)
	E2.set_angulo_pressao(20)
	E3.set_angulo_pressao(20)
	E4.set_angulo_pressao(20)
	E1.set_largura(b1)
	E2.set_largura(b1)
	E3.set_largura(b2)
	E4.set_largura(b2)
	E1.set_material("steel 1045")
	E2.set_material("steel 1045")
	E3.set_material("steel 1045")
	E4.set_material("steel 1045")

	log.write(E1)
	#log.write("\n\n")
	log.write(E2)
	#log.write("\n\n")
	log.write(E3)
	#log.write("\n\n")
	log.write(E4)
	#log.write("\n\n")


	A.set_rotacao(250) # Em rpm
	B.set_rotacao(250*r1)
	C.set_rotacao(250*r1*r2)
	A.set_potencia(P) # Em kW
	B.set_potencia(P) # Em kW
	C.set_potencia(P) # Em kW
	
	# M1.set_razao_engrenagem(5) # Uma transmissão de 1 volta para 5 no primeiro conjunto
	# M2.set_razao_engrenagem(6) # Uma transmissão de 1 volta para 6 no segundo conjunto

	# Angulo de pressao do primeiro conjunto de 20 graus, e será uma razão de 5
	Verbose = False
	tensoes1 = EsforcoNoDente.get_tensoes_dente(E1, E2, verbose = Verbose)
	tensoes2 = EsforcoNoDente.get_tensoes_dente(E2, E1, verbose = Verbose)
	tensoes3 = EsforcoNoDente.get_tensoes_dente(E3, E4, verbose = Verbose)
	tensoes4 = EsforcoNoDente.get_tensoes_dente(E4, E3, verbose = Verbose)
	FatorS = FatoresDentePermitido.S()
	#log.write("tensoes1 = " + str(tensoes1))
	#log.write("flex = %.2f cont = %.2f" % (tensoes1[4], tensoes1[5]))

	#log.write("tensoes2 = " + str(tensoes2))
	#log.write("flex = %.2f cont = %.2f" % (tensoes2[4], tensoes2[5]))

	#log.write("tensoes3 = " + str(tensoes3))
	#log.write("flex = %.2f cont = %.2f" % (tensoes3[4], tensoes3[5]))

	#log.write("tensoes4 = " + str(tensoes4))
	#log.write("flex = %.2f cont = %.2f" % (tensoes3[4], tensoes3[5]))

	# E2.set_dentes(EngExtDentesRetos.menor_pinhao(r = 5, phi_n = 20))
	# Angulo de pressao do segundo conjunto de 20 graus, e terá uma razão de 6
	# E4.set_dentes(EngExtDentesRetos.menor_pinhao(r = 6, phi_n = 20))
'''