!##################################################
!	Programa principal
!##################################################

module mod_tubos
	use geraRede
	use mod_rndgen
	use mod_tools
	implicit none
	
	integer, allocatable, private :: listAdjN2(:), listAdjTotal(:)
	integer, allocatable, private :: grauN2(:), grauTotal(:)
	integer, allocatable, private :: auxN2(:), auxTotal(:)
	integer, allocatable, private :: matrizTotal(:,:)
	integer, allocatable, private :: stubs_pl(:)	
	type(rndgen), private :: gen3
	
	integer, private :: N2
	integer, private :: i1, i2, i3
	integer, private :: N_pontas
	integer, private :: conta_rejec	
	contains
	
subroutine PL_Tubos(this, N, kMin, kMax, gama, semente1, plus_xN, f_Tubos, PLA, alpha)
!#######################################################################
! Variaveis de entrada
!#######################################################################
	class(grafo_PL_UCM) :: this
	integer, intent(in) :: N
	integer, intent(in) :: kMin
	real(dp), intent(in) :: kMax
	real(dp), intent(in) :: gama
	integer, intent(in) :: semente1
	integer :: semente 
	real(dp), intent(in) :: plus_xN, f_Tubos
	logical, intent(in) :: PLA
	real(dp), intent(in) :: alpha
!#######################################################################
! Auxiliares
!#######################################################################
	integer, allocatable :: pontas(:), tubos_puros(:)
!#######################################################################
!
!#######################################################################
!	real(dp), allocatable :: aux_tubos(:), aux_tubos_prov(:)	
	integer :: i1, i2, i3, last_stub
	integer :: tam_tubo
	integer :: no2_1, no2_2
	integer :: pl_node1, pl_node2
	integer :: edgeTotal
	integer :: nao_nulos
	logical :: tamanho1
	
!#######################################################################	
!	* plus_xN > 1 e mostra qual a proporcao de N sera adicionada
!	  em nos de grau 2 e sera adicionada a power law original;
!		
!#######################################################################
	
!#######################################################################	
!	

!#######################################################################	
!		Rede PL eh criada

		semente  = semente1
		call this%iniciaGrafo(N)
		call this%inicia(kMin, kMax, gama, semente)
		
		semente = semente1
		call this%liga(semente, .True.)


		this%degMin = minval(this%deg)		
		this%degMax = maxval(this%deg)
				
!#######################################################################
		
		call sub_classifica_clusters(this,.False., 000, 'sem_arquivo.dat')



		write(*,*) 'Apos gerar a rede, a componente gigante e ', comp_gigante

if(plus_xN > 0.0_dp)then
		!###############################################################
		!	Guardo o numero de nos de grau 2 e o numero de pontas
			N2 = int(1.0_dp * plus_xN * this%nodes)
			N_pontas = int(1.0_dp * f_Tubos * N2)
			
			if(mod(N_pontas, 2) > 0) N_pontas = N_pontas + 1
		!
		!###############################################################

		!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		!	Ateh aqui parece ok
		!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		if(allocated(matrizTotal)) deallocate(matrizTotal)
			allocate(matrizTotal(this%edge, 2))
		
		do i1 = 1, this%edge
			matrizTotal(i1, 1) = this%matriz(i1, 1)
			matrizTotal(i1, 2) = this%matriz(i1, 2)			
		enddo
		edgeTotal = this%edge

		!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		!	Ateh aqui parece ok tambem
		!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		
		deallocate(this%matriz)
		
		allocate(this%matriz(this%edge + (2 * N2 + N_pontas), 2))

		!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		!	Aqui foi acrescentado um numero de links N_pontas acima
		!	do necessario.
		!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		
		do i1 = 1, this%edge
			this%matriz(i1, 1) = matrizTotal(i1, 1)
			this%matriz(i1, 2) = matrizTotal(i1, 2)			
		enddo
		
		deallocate(matrizTotal)
			
		
!#######################################################################						
!	Listas dos deg = 2 que ligarao ao nucleo.
	
	if(allocated(listAdjN2)) deallocate(listAdjN2)
		allocate(listAdjN2(2 * N2))

		if(allocated(grauN2)) deallocate(grauN2)
			allocate(grauN2((this%nodes +1): (this%nodes + N2)))
			grauN2 = 0

!#######################################################################

		if(allocated(grauTotal)) deallocate(grauTotal)
			allocate(grauTotal(this%nodes))
			
		do i1 = 1, this%nodes
			grauTotal(i1) = this%deg(i1)
		enddo
			
		if(allocated(this%deg)) deallocate(this%deg)
			allocate(this%deg(this%nodes+N2))
			this%deg = 0
			
		do i1 = 1, this%nodes
			 this%deg(i1) = grauTotal(i1)
		enddo			
		
		deallocate(grauTotal)
							
		if(allocated(auxN2)) deallocate(auxN2)
			allocate(auxN2((this%nodes + 1):(this%nodes + N2)))					
		
		auxN2(this%nodes + 1) = 1
		
		do i1 = this%nodes + 2, this%nodes + N2
			auxN2(i1) = auxN2(i1 - 1) + 2
		enddo 
		
!#######################################################################			
				
	
!		if(allocated(pontas)) deallocate(pontas)
!			allocate(pontas(this%nodes+1, this%nodes + (2 * N_tubos)))

!		if(allocated(tubos_puros)) deallocate(tubos_puros)
!			allocate(tubos_puros(this%nodes+ (2 * N_tubos) + 1, this%nodes + N2))					
						
!#######################################################################		
	
!		semente = semente1
		
!		call liga_tubos(semente, .False.)

!#######################################################################		
				
		semente = semente1
		
		call liga_tubos(this, semente, .True., alpha)

!#######################################################################				

		deallocate(this%aux)
		allocate(this%aux(this%nodes + N2))
		
		this%aux(1) = 1
		do i1 = 2, size(this%aux)
			this%aux(i1) = this%aux(i1 - 1) + this%deg(i1 - 1)
		enddo
		
		this%sumDeg = sum(this%deg)
		
		deallocate(this%listAdj)
		allocate(this%listAdj(this%sumDeg))

		if(allocated(grauTotal)) deallocate(grauTotal)
			allocate(grauTotal(this%nodes + N2))
			grauTotal = 0
			
		do i1 = 1, this%edge
			
			this%listAdj(this%aux(this%matriz(i1, 1)) + grauTotal(this%matriz(i1,1))) = this%matriz(i1, 2)
			grauTotal(this%matriz(i1, 1)) = grauTotal(this%matriz(i1, 1)) + 1
			
			this%listAdj(this%aux(this%matriz(i1, 2)) + grauTotal(this%matriz(i1,2))) = this%matriz(i1, 1)
			grauTotal(this%matriz(i1, 2)) = grauTotal(this%matriz(i1, 2)) + 1
		
		enddo
	
		deallocate(grauTotal)

!#######################################################################
!	Agora conto se algum no nao foi ligado
!#######################################################################
		
		nao_nulos = 0
		
		do i1 = 1, this%nodes + N2
			if(this%deg(i1) > 0) nao_nulos = nao_nulos + 1
		enddo
		
		write(*,*) 'Sobraram ', this%nodes + N2 - nao_nulos, ' sitios sem qualquer ligacao.'

		this%nodes = nao_nulos
		
		call sub_classifica_clusters(this,.False., 000, 'sem_arquivo.dat')

endif

end subroutine

!#######################################################################
		
subroutine liga_tubos(this,semente, conecta, alpha)
type(grafo_PL_UCM) :: this
integer :: semente
logical, intent(in) :: conecta
real(dp), intent(in) :: alpha
real(dp) :: prob
integer :: viz1, viz2
integer :: N_tubos_prov, size_tubos_puros
integer :: nao_nulos
integer :: stubs, stub1, stub2, no1, no2
integer, allocatable :: grauTotal(:)
integer :: grauMaxTotal
integer, allocatable :: stubs_N2(:)


!#######################################################################
	
!#######################################################################			
!	Muito importante!

call gen3%init(semente)


!#######################################################################
!	O nucleo ja existe, basta conectar as pontas a eles.
!#######################################################################

if(alpha == 0.0_dp) PLA = .False.


write(*,*) "o valor de alpha e : ", alpha
if(PLA)then
	write(*,*) "escolheu PLA"
	write(*,*) "grau maximo eh ", this%degMax
	do i1 = this%nodes + 1, this%nodes + N_pontas

retete:		do
				i2 = gen3%int(1, this%nodes)

				if(lista_de_clusters(i2) /= i_comp_gigante) cycle retete

				prob = gen3%rnd()

!				if((1.0_dp * (this%deg(i2)/this%degMax))**(alpha) < prob) cycle retete

				if(((1.0_dp * this%deg(i2))/(1.0_dp * this%degMax))**(alpha) < prob) cycle retete

				exit retete

			enddo retete
			
			this%edge = this%edge + 1
			
			this%matriz(this%edge, 1) = i1
			this%matriz(this%edge, 2) = i2
									
			listAdjN2(auxN2(i1) + this%deg(i1)) = i2
			
			this%deg(i1) = this%deg(i1) + 1
			this%deg(i2) = this%deg(i2) + 1
			
			this%degMax = max(this%degMax,this%deg(i2))
			
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!	Tinha algo errado aqui. A ordem estava assim
!
!				this%deg(i1) = this%deg(i1) + 1
!				listAdjN2(auxN2(i1) + this%deg(i1)) = i2
!	
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@			
	
	enddo
else

	do i1 = this%nodes + 1, this%nodes + N_pontas

retete1:	do
				i2 = gen3%int(1, this%nodes)
				if(lista_de_clusters(i2) /= i_comp_gigante) cycle retete1
				exit retete1
			enddo retete1
			
			this%edge = this%edge + 1
			
			this%matriz(this%edge, 1) = i1
			this%matriz(this%edge, 2) = i2

			listAdjN2(auxN2(i1) + this%deg(i1)) = i2
			this%deg(i1) = this%deg(i1) + 1
			this%deg(i2) = this%deg(i2) + 1						

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!	Tinha algo errado no PLA 1, sim. A ordem estava diferente.
!
!				this%deg(i1) = this%deg(i1) + 1
!				listAdjN2(auxN2(i1) + this%deg(i1)) = i2
!	
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@			


	enddo

endif

this%degMax = maxval(this%deg)

!#######################################################################
!	Pontas conectadas. Agora temos que ligar os outros deg = 2
!	Para isto, vamos preparar uma lista de pontas que saem dos
!	nos de grau 2.
!#######################################################################


if(allocated(stubs_N2)) deallocate(stubs_N2)
	allocate(stubs_N2((2 * N2 - N_pontas)))
	
	i2 = 1
	do i1 = this%nodes + 1, this%nodes + N_pontas
		stubs_N2(i2) = i1
		i2 = i2 + 1
	enddo
	
	do i1 = this%nodes + N_pontas + 1, this%nodes + N2

		do i3 = 1, 2	
			stubs_N2(i2) = i1
			i2 = i2 + 1
		enddo
	enddo
	
!#######################################################################
!	A lista das pontas esta preparada, basta agora liga-las
!#######################################################################

stubs = size(stubs_N2)

willow: do while(stubs > 0)

	stub1 = gen3%int(1, stubs)	
	no1 = stubs_N2(stub1)
	
	conta_rejec = 0
	
willow2:	do

				if(conta_rejec > 100 * stubs) exit willow
				
				stub2 = gen3%int(1, stubs)	
				
				no2 = stubs_N2(stub2)
			
				if(no1 == no2)then
					conta_rejec = conta_rejec + 1
					cycle willow2
				endif
				
				do i1 = auxN2(no1), auxN2(no1) + this%deg(no1) - 1
					if(listAdjN2(i1) == no2)then
						conta_rejec = conta_rejec + 1
						cycle willow2
					endif	
				enddo
					
				exit willow2
		
			enddo willow2
	
	this%edge = this%edge + 1
	
	this%matriz(this%edge, 1) = no1
	this%matriz(this%edge, 2) = no2
	
	listAdjN2(auxN2(no1) + this%deg(no1)) = no2

	this%deg(no1) = this%deg(no1) + 1

	listAdjN2(auxN2(no2) + this%deg(no2)) = no1

	this%deg(no2) = this%deg(no2) + 1	
		
	if(stub1 > stub2)then
		
		no1 = stubs_N2(stubs)
		stubs_N2(stub1) = no1
		stubs = stubs - 1
		
		no2 = stubs_N2(stubs)
		
		stubs_N2(stub2) = no2
		stubs = stubs - 1		
	else
		
		no2 = stubs_N2(stubs)
		
		stubs_N2(stub2) = no2
		stubs = stubs - 1
		
		no1 = stubs_N2(stubs)
		
		stubs_N2(stub1) = no1
		stubs = stubs - 1			
	endif
		
enddo willow

	deallocate(auxN2)
	deallocate(listAdjN2)

write(*,*) " "
write(*,*) 'Sobraram ', stubs, ' stubs.'
write(*,*) " "
!#######################################################################
!	Se saiu do willow, nao tem muito mais o que fazer.
!	Agora eh atualizar as listas
!	Para saber se sobrou alguem, contamos os itens nao nulos
!	da lista de grau.
!	Algo que pode ocorrer eh a formacao de loops.
!	Isso vai acontecer, portanto, eh importante recalcularmos
!	a componente gigante.
!	Temos que atualizar a lista de adjacencias e a auxiliar
!#######################################################################	
this%degMin = minval(this%deg)
this%degMax = maxval(this%deg)


end subroutine


end module


program main
	use mod_tubos
	use mod_rndgen
	use geraRede
	use dynamics !, only: contactProcess
	use getters_metricos_dinamicos
	use mod_tools
	
	implicit none
	
	!##########################################
	!	Parametros internos ao probrama
	!##########################################	
!	integer, parameter :: dp = kind(0.0d0)

	!##########################################
	!	Medidor de tempo computacional
	!##########################################	

!	real(dp) :: start, finish, start2, finish2, tempSis

	!##########################################
	!	Substrato ou rede e seus
	!	parametros
	!##########################################	


	type(grafo_PL_UCM) :: rede
		
	integer :: i8, i9, i10, i11, i21, i22, i100, i101
	integer :: i7, j7
	!##########################################
	!	Semente para o gerador de numeros
	!	pseudo aleatorios.
	!##########################################	
	integer :: seed, label1
	character(len=255) :: arquivo1, diretorio, diretorioPai, subDiretorio, comp_gigante_char,&
	id_amostra_Char, nInfMax_n_redundante_char, nInfMax_redundante1_char, nInfMax_redundante2_char
	character(len=300) :: arquivo, dist_tubos_file, dist_grau_file, &
	dist_grau2_grau1_file, dist_grau2_grau2_file, dist_grau2_grauGT2_file,&
	dist_stubs_grau1, dist_stubs_grau2

	real(dp), allocatable :: conta_Vizinhos_Nos_Grau2(:), conta_Vizinhos_Nos_Grau1(:)

	integer :: conta_stubs_Grau2, conta_stubs_Grau1
	real(dp) :: Xi_Max_prov, Xi_last, lambda_C, lambda_last, rho_QS_last, t_LS_last
	integer :: io_Status

	seed =347361833
	imuniza_tubos = .False.

!##############################################
!	Inicializacao da rotina
!##############################################


!##########################################################################################
!		Aloca lista de tamanhos
!##########################################################################################

	call get_Arg_PL_UCM_com_Tubos_RBS_Ensemble_PLA_Random()

	diretorio=''
	
	write(id_amostra_char, '(I0)') id_amostra

!##########################################################################################
!				Inicia grafo
!##########################################################################################

			
			write(id_amostra_Char, '(I0)') id_amostra


			call PL_Tubos(rede, N, kMin2, kMax2, gama2, seed_Ensemble, plus_xN,f_Tubos, PLA, alpha)

			write(*,*) " "
			write(*,*) "Rede inicializada"


!#######################################################################	
!			call sub_classifica_clusters(rede,.False., 000, 'sem_arquivo.dat')
			
			write(*,*) "O valor da componente gigante raiz e ", comp_gigante
			
			label1=50
			arquivo1 = 'tam_tubos.dat'
			call statistic_Tubes_c_detalhe(rede, label1, arquivo1)
			call calcula_k_nn(rede,.True., 222, 'knn_vs_k.dat')
			call clustering(rede,.True., 223, 'clustering_vs_k.dat')
!			call calcula_distancias_no_nucleo(rede)			
!#######################################################################			
 
 			write(*,*) " "
			write(*,*) "Numero de nos= ", rede%nodes, "grau por no= ", rede%degMean
			write(*,*) " "
			

			write(*,*) " "
			write(*,*) "Rede gerada com sucesso"
			write(*,*) " "


!			open(139, file='nos_imunizados.dat', status='unknown')	
			
				write(*,*) ' '
				write(*,*) 'O tamanho da componente gigante e: ', comp_gigante
				write(*,*) ' '
			

			!#######################################################################
!				deallocate(rede%aux)
!				deallocate(rede%listAdj)
!				deallocate(rede%deg)
			!#######################################################################
			
			call allocaMedias_An_QS2(rede)
			
			write(comp_gigante_char, '(I0)') comp_gigante
			write(nInfMax_n_redundante_char, '(I0)') nInfMax_n_redundante
			write(nInfMax_redundante2_char, '(I0)') nInfMax_redundante2
						
			open(unit=88, file = trim(adjustl(diretorio))//'Xi_Max_vs_lambda'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(comp_gigante_char))//'.dat', status='unknown')
			open(unit=89, file = trim(adjustl(diretorio))//'rho_vs_lambda_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(comp_gigante_char))//'.dat', status = 'unknown')
			open(unit=90, file = trim(adjustl(diretorio))//'Xi_vs_lambda_amostra_'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(comp_gigante_char))//'.dat', status = 'unknown')			
			open(unit=91, file = trim(adjustl(diretorio))//'n_vs_QS_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(comp_gigante_char))//'.dat', status = 'unknown')
!			open(unit=93, file = trim(adjustl(diretorio))//'conexoes_Rede.csv', status = 'unknown')
			open(unit=94, file = trim(adjustl(diretorio))//'tempoVidaMedia_vs_lambda_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(comp_gigante_char))//'.dat', status = 'unknown')
			open(unit=95, file = trim(adjustl(diretorio))//'dev_rho_MedioQS_vs_lambda_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(comp_gigante_char))//'.dat', status = 'unknown')
			open(unit=101, file = trim(adjustl(diretorio))//'rho_Max_vs_lambda'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(comp_gigante_char))//'.dat', status='unknown')
			open(unit=102, file = trim(adjustl(diretorio))//'S_vs_lambda'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(comp_gigante_char))//'.dat', status='unknown')		


!#######################################################################
!			Bloco de arquivos para analise nao redundante
			
			open(unit=103, file = trim(adjustl(diretorio))//'S_n_redundante_vs_lambda'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_n_redundante_char))//'.dat', status='unknown')
			open(unit=104, file = trim(adjustl(diretorio))//'Xi_n_redundante_vs_lambda_amostra_'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_n_redundante_char))//'.dat', status = 'unknown')
			open(unit=105, file = trim(adjustl(diretorio))//'Xi_n_redundante_Max_vs_lambda'//trim(adjustl(id_amostra_Char))//'_CG_n_redun_'//trim(adjustl(nInfMax_n_redundante_char))//'.dat', status='unknown')
			open(unit=106, file = trim(adjustl(diretorio))//'rho_n_redundante_vs_lambda_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_n_redundante_char))//'.dat', status = 'unknown')
			open(unit=107, file = trim(adjustl(diretorio))//'n_n_redundantes_vs_QS_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_n_redundante_char))//'.dat', status = 'unknown')
			open(unit=108, file = trim(adjustl(diretorio))//'tempoVidaMedia_n_redundantes_vs_lambda_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_n_redundante_char))//'.dat', status = 'unknown')
			open(unit=109, file = trim(adjustl(diretorio))//'dev_rho_n_redundantes_MedioQS_vs_lambda_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_n_redundante_char))//'.dat', status = 'unknown')
			open(unit=110, file = trim(adjustl(diretorio))//'rho_Max_n_redundantes_vs_lambda'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_n_redundante_char))//'.dat', status='unknown')

!			open(unit=111, file = trim(adjustl(diretorio))//'S_redundante1_vs_lambda'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante1_char))//'.dat', status='unknown')
!			open(unit=112, file = trim(adjustl(diretorio))//'Xi_redundante1_vs_lambda_amostra_'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante1_char))//'.dat', status = 'unknown')
!			open(unit=113, file = trim(adjustl(diretorio))//'Xi_redundante1_Max_vs_lambda'//trim(adjustl(id_amostra_Char))//'_CG_n_redun_'//trim(adjustl(nInfMax_redundante1_char))//'.dat', status='unknown')
!			open(unit=114, file = trim(adjustl(diretorio))//'rho_redundante1_vs_lambda_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante1_char))//'.dat', status = 'unknown')
!			open(unit=115, file = trim(adjustl(diretorio))//'n_redundantes1_vs_QS_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante1_char))//'.dat', status = 'unknown')
!			open(unit=116, file = trim(adjustl(diretorio))//'tempoVidaMedia_redundantes1_vs_lambda_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante1_char))//'.dat', status = 'unknown')
!			open(unit=117, file = trim(adjustl(diretorio))//'dev_rho_redundantes1_MedioQS_vs_lambda_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante1_char))//'.dat', status = 'unknown')
!			open(unit=118, file = trim(adjustl(diretorio))//'rho_Max_redundantes1_vs_lambda'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante1_char))//'.dat', status='unknown')

			open(unit=119, file = trim(adjustl(diretorio))//'S_redundante2_vs_lambda'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante2_char))//'.dat', status='unknown')
			open(unit=120, file = trim(adjustl(diretorio))//'Xi_redundante2_vs_lambda_amostra_'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante2_char))//'.dat', status = 'unknown')
			open(unit=121, file = trim(adjustl(diretorio))//'Xi_redundante2_Max_vs_lambda'//trim(adjustl(id_amostra_Char))//'_CG_n_redun_'//trim(adjustl(nInfMax_redundante2_char))//'.dat', status='unknown')
			open(unit=122, file = trim(adjustl(diretorio))//'rho_redundante2_vs_lambda_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante2_char))//'.dat', status = 'unknown')
			open(unit=123, file = trim(adjustl(diretorio))//'n_redundantes2_vs_QS_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante2_char))//'.dat', status = 'unknown')
			open(unit=124, file = trim(adjustl(diretorio))//'tempoVidaMedia_redundantes2_vs_lambda_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante2_char))//'.dat', status = 'unknown')
			open(unit=125, file = trim(adjustl(diretorio))//'dev_rho_redundantes2_MedioQS_vs_lambda_amostra'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante2_char))//'.dat', status = 'unknown')
			open(unit=126, file = trim(adjustl(diretorio))//'rho_Max_redundantes2_vs_lambda'//trim(adjustl(id_amostra_Char))//'_CG_'//trim(adjustl(nInfMax_redundante2_char))//'.dat', status='unknown')
			open(unit=888, file='lista_conexoes_com_tubos.csv', status='unknown')

!			write(888,*) 'Source', ',' , 'Target'
			do i7 = 1, rede%edge
				if(lista_de_clusters(rede%matriz(i7,1)) /= i_comp_gigante) cycle
				write(888,*) rede%matriz(i7,1), ',', rede%matriz(i7,2)
			enddo
			close(888)

!#######################################################################
!	Estatistica pros nos de grau 2

	if(pok(2) > 0.0_dp)then
		if(allocated(conta_Vizinhos_Nos_Grau2)) deallocate(conta_Vizinhos_Nos_Grau2)
			allocate(conta_Vizinhos_Nos_Grau2(rede%degMax))
			
			conta_Vizinhos_Nos_Grau2 = 0.0_dp

		do i100 = 1, rede%nodes
			
			if(lista_de_clusters(i100) /= i_comp_gigante) cycle
			
			if(rede%deg(i100) /= 2) cycle
			
			do i101 = rede%aux(i100), rede%aux(i100) + 1
				conta_Vizinhos_Nos_Grau2(rede%deg(rede%listAdj(i101))) = conta_Vizinhos_Nos_Grau2(rede%deg(rede%listAdj(i101))) + 1.0_dp
			enddo
		
		enddo
		
		conta_stubs_Grau2 = 2 * nos_grau2
		
		conta_Vizinhos_Nos_Grau2 = 1.0_dp * conta_Vizinhos_Nos_Grau2/conta_stubs_Grau2

		open(663, file='probabilidade_k_dado_grau2.dat', status='unknown')
		
		do i100 = 1, size(conta_Vizinhos_Nos_Grau2)
			if(conta_Vizinhos_Nos_Grau2(i100) > 0.0_dp)then
				write(663, *) i100, conta_Vizinhos_Nos_Grau2(i100) 
			endif
		enddo
		
		close(663)			
	endif
!#######################################################################
				
!				call alocaQS_IOGA(rede)
!				call allocaMedias_An_QS2(rede)

!#######################################################################


!#######################################################################
			
		do i7 = 1, npts

				if(lambda == 0.0_dp) lambda = dlambda
				
!				call condicaoInicial_CG_IOGA(rede)
				
				call condicaoInicial_CG2(rede, imuniza_tubos)							
				! imunizacao foi escolhido .False., portanto, tudo bem.
	
			!###########################################################
			!	Calcula tamanho da rede nao imunizada
			!###########################################################
							

			

			!###########################################################
			!	Escreve nos imunizados
			!###########################################################

!			do i21 = 1, rede%nodes
!				if(sigma(i21) == 2)then
!					write(139,*) i21, 2
!				endif
!			enddo
!			close(139)		

!			if ( npts == 1) then										!Aqui gasta pouco tempo mesmo
!				do i8 = 1, 913
!					esquenta = gen%rnd()
!				enddo
!			endif

			call sisProcessRBS1(rede,.True., seed_Ensemble)
!		    call sisProcessRBS_IOGA(rede, .True., seed_Ensemble)
		
!#############################################################################################
!			Aqui escrevo os resultados da simulacao
!#############################################################################################

			!#####################################################################
			!	Esse e um chute grosseiro para lambda_C, baseado no tamanho
			!       de Xi, porem, devido as flutuacoes, este procedimento nao
			!	e acurado.
			!#####################################################################
			read(88,*, iostat = io_Status) lambda_C, Xi_Max_prov

			if(io_Status < 0)then
				Xi_Max_prov = 0.0_dp; lambda_C = lambda
			endif
			
			rewind(88)
!#############################################################################################

!			write(*,*) 'Os candidatos atuais a lambda_C e susceptibilidade maxima sao ', lambda_C, ' e ', Xi_Max_prov

			if(Xi > Xi_Max_prov)then
				 lambda_C = lambda
		       		 Xi_Max_prov = Xi			
				 write(88, *) lambda_C, Xi_Max_prov
				 write(101, *) lambda_C, rho_medioQS
			endif
			rewind(88)
			rewind(101)
			!#####################################################################
			!	Aqui termino o bloco onde vou implementar essa tentativa
			!	de achar o pico na malicia. rs
			!	Quando ficar bom, vou por numa subrotina.
			!#####################################################################
			write(102, *) lambda, S_Schanon
			write(89,*) lambda, rho_medioQS
			write(95,*) lambda, dev_rho_medioQS
			write(90, *) lambda, Xi

!			do i8 = 1, nInfMax
!				if(Pn_QS(i8) > 0.0_dp)then
!					write(91, *) i8, Pn_QS(i8) 
!				endif
!			enddo


!#######################################################################
!	Estatistica parte n redundante
!#######################################################################


			read(105,*, iostat = io_Status) lambda_C, Xi_Max_prov

			if(io_Status < 0)then
				Xi_Max_prov = 0.0_dp; lambda_C = lambda
			endif
			
			rewind(105)
!#############################################################################################

!			write(*,*) 'Os candidatos atuais a lambda_C e susceptibilidade maxima sao ', lambda_C, ' e ', Xi_Max_prov

			if(Xi_n_redundante > Xi_Max_prov)then
				 lambda_C = lambda
				 Xi_Max_prov = Xi_n_redundante			
				 write(105, *) lambda_C, Xi_Max_prov
				 write(110, *) lambda_C, rho_medioQS_n_redundante
			endif
			rewind(105)
			rewind(110)
			!#####################################################################
			!	Aqui termino o bloco onde vou implementar essa tentativa
			!	de achar o pico na malicia. rs
			!	Quando ficar bom, vou por numa subrotina.
			!#####################################################################

			write(103, *) lambda, S_Schanon_n_redundante
			write(104, *) lambda, Xi_n_redundante
			write(106,*) lambda, rho_medioQS_n_redundante
			write(109,*) lambda, dev_rho_medioQS_n_redundante

!			do i8 = 1, nInfMax_n_redundante
!				if(Pn_QS_n_redundante(i8) > 0.0_dp)then
!					write(107, *) i8, Pn_QS_n_redundante(i8) 
!				endif
!			enddo


!#######################################################################
!	Estatistica parte redundante 2
!#######################################################################

			read(121,*, iostat = io_Status) lambda_C, Xi_Max_prov

			if(io_Status < 0)then
				Xi_Max_prov = 0.0_dp; lambda_C = lambda
			endif
			
			rewind(121)
!#############################################################################################

!			write(*,*) 'Os candidatos atuais a lambda_C e susceptibilidade maxima sao ', lambda_C, ' e ', Xi_Max_prov

			if(Xi_redundante2 > Xi_Max_prov)then
				 lambda_C = lambda
				 Xi_Max_prov = Xi_redundante2			
				 write(121, *) lambda_C, Xi_Max_prov
				 write(126, *) lambda_C, rho_medioQS_redundante2
			endif
			rewind(121)
			rewind(126)
			
			!#####################################################################
			!	Aqui termino o bloco onde vou implementar essa tentativa
			!	de achar o pico na malicia. rs
			!	Quando ficar bom, vou por numa subrotina.
			!#####################################################################

			write(119, *) lambda, S_Schanon_redundante2
			write(120, *) lambda, Xi_redundante2
			write(122,*) lambda, rho_medioQS_redundante2
			write(125,*) lambda, dev_rho_medioQS_redundante2

!			do i8 = 1, nInfMax_n_redundante
!				if(Pn_QS_n_redundante(i8) > 0.0_dp)then
!					write(107, *) i8, Pn_QS_n_redundante(i8) 
!				endif
!			enddo



!#############################################################################################
!			Aqui concluo
!#############################################################################################

			lambda = lambda + dlambda
	enddo


	close(88)
	close(89)
	close(90)
	close(91)
	close(94)
	close(95)
	close(101)
	
	close(102)
	close(103)
	close(104)
	close(105)
	close(106)
	close(107)
	close(108)
	close(109)
	close(110)
	
!	close(111)
!	close(112)
!	close(113)
!	close(114)
!	close(115)
!	close(116)
!	close(117)
!	close(118)

	close(119)	
	close(120)
	close(121)
	close(122)
	close(123)
	close(124)
	close(125)
	close(126)
	
end program
