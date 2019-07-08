module getters_metricos_dinamicos
	use dynamics
	use geraRede
	implicit none
!####################################################################################
!		Pede parametros RBS ao usuario
!####################################################################################

contains

		subroutine promptUserRBS()
				write(*,*) " Digite a porcentagem de nos infectados inicialmente"
				write(*,*) ""
				read(*,*) nInf0
						
				write(*,*) " Digite o tempo de relaxacao da epidemia"
				read(*,*) t_relax
				write(*,*) " "

				write(*,*) " Digite o numero de intervalos de tempo para calcular as medias da epidemia"
				write(*,*) " "
				read(*,*) t_Media
				write(*,*) " "
				
				tMax = t_relax + t_Media			
		end subroutine


!#####################################################################################
!	Get Arg RBS
!#####################################################################################

		subroutine get_Arg_RRN_RBS()
			integer :: N0
			character(len=20) :: leitura, buffer
			
			N0 = iargc()
			
			if(N0 /= 10) stop " Numero de nos N, ki, nIn0 real <= 1, lambda, lambda_i, lambda_f, dlambda, npts, t_relax integer, t_Media integer"			

				call getarg(1, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) N 			
				

				call getarg(2, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) ki 

				call getarg(3, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) nInf0 

				call getarg(4, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) lambda 
				
				call getarg(5, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) lambda_i 
				
				call getarg(6, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) lambda_f
				 
				call getarg(7, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) dlambda

				call getarg(8, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) npts 
				
				call getarg(9, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) t_relax 
				
				
				call getarg(10, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) t_Media 
				
				
				tMax = t_relax + t_Media		
		end subroutine


!#####################################################################################

		subroutine get_Arg_RRN_RBS_Mistura()
			implicit none
			integer :: N0
			character(len=20) :: leitura, buffer
			
			N0 = iargc()
			
			if(N0 /= 12) stop " Numero de nos N integer, k1 integer, k2 integer, p1 real <=1, nIn0 real <= 1, lambda, lambda_i, lambda_f, dlambda, npts, t_relax integer, t_Media integer"			
				
				
				call getarg(1, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) p1  			
				write(*,*) p1
				
				call getarg(2, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) nInf0  			
				write(*,*) nInf0
				
				call getarg(3, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) lambda 				
				write(*,*) lambda
				
				call getarg(4, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) lambda_i 				
				write(*,*) lambda_i
				
				call getarg(5, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) lambda_f
				write(*,*) lambda_f 

				call getarg(6, leitura)
				buffer = trim(adjustl(leitura))
		 		read(buffer,*) dlambda   			
				write(*,*) dlambda
				
				call getarg(7, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) N
				write(*,*) N
				
				call getarg(8, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) k1
				write(*,*) k1
				
				call getarg(9, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) k2
				write(*,*) k2
				
				call getarg(10, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) npts
				write(*,*) npts  

				
				call getarg(11, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) t_relax
				write(*,*) t_relax  
								
				call getarg(12, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) t_Media
				write(*,*) t_Media  
				
				
				tMax = t_relax + t_Media		
		end subroutine


!####################################################################################
!
!####################################################################################

	subroutine get_Arg_Bimod_RBS_Ensemble()
		implicit none
		integer :: N0, recuperaDin_int
		character(len=20) :: leitura, buffer
		
		N0 = iargc()
		
		if(N0 /= 15) stop " Numero de nos N integer, k1 integer, k2 integer, p1 real <=1, nIn0 real <= 1, lambda, lambda_i, lambda_f, dlambda, npts, t_relax integer, t_Media integer, id_amostra integer, seed integer, recuperaDin_int integer"
			
			
			call getarg(1, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) p1  			
			write(*,*) p1
			
			call getarg(2, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) nInf0  			
			write(*,*) nInf0
			
			call getarg(3, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda 				
			write(*,*) lambda
			
			call getarg(4, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda_i 				
			write(*,*) lambda_i
			
			call getarg(5, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda_f
			write(*,*) lambda_f 

			call getarg(6, leitura)
			buffer = trim(adjustl(leitura))
	 		read(buffer,*) dlambda   			
			write(*,*) dlambda
			
			call getarg(7, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) N
			write(*,*) N
			
			call getarg(8, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) k1
			write(*,*) k1
			
			call getarg(9, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) k2
			write(*,*) k2
			
			call getarg(10, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) npts
			write(*,*) npts  

			
			call getarg(11, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) t_relax
			write(*,*) t_relax  
							
			call getarg(12, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) t_Media
			write(*,*) t_Media  
	
			call getarg(13, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) id_amostra
			write(*,*) id_amostra

			call getarg(14, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) seed_Ensemble
			write(*,*) seed_Ensemble
			
			call getarg(15, leitura)
                        buffer = trim(adjustl(leitura))
                        read(buffer,*) recuperaDin_int
                        write(*,*) recuperaDin_int

                        if(recuperaDin_int == 1)then
                                recuperaDinamica = .True.
                        elseif(recuperaDin_int == 0)then
                                recuperaDinamica = .False.
                        else
                                stop "Sao aceitos apenas os valores 1 ou 0"
                        endif

			tMax = t_relax + t_Media		
	end subroutine

!#######################################################################
!	Obtem parametros para rede PL_UCM com Tubos
!	e condicao RBS em ensemble
!#######################################################################

	subroutine get_Arg_PL_UCM_com_Tubos_RBS_Ensemble()
		implicit none
		integer :: N0, recuperaDin_int, ressorteia_int
		real(dp) :: omega
		character(len=20) :: leitura, buffer
					
		N0 = iargc()
		
		if(N0 /= 17) stop " Numero de nos N integer, prob_tube real, &
		kMin2 integer, kMax2 integer, gamma2 real > 0, nIn0 real <= 1,&
		 lambda, lambda_i, lambda_f, dlambda, npts, t_relax integer, &
		  t_Media integer, id_amostra integer, seed integer, recuperaDin_int integer, integer sorteia"			
			
			

			call getarg(1, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) N
			write(*,*) N

			call getarg(2, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) probTube  			
			write(*,*) probTube
						
			call getarg(3, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) kMin2
			write(*,*) kMin2
			
			call getarg(4, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) kMax2
			write(*,*) kMax2
			
			call getarg(5, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) gama2
			write(*,*) gama2
			
			call getarg(6, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) nInf0  			
			write(*,*) nInf0
			
			call getarg(7, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda 				
			write(*,*) lambda
			
			call getarg(8, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda_i 				
			write(*,*) lambda_i
			
			call getarg(9, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda_f
			write(*,*) lambda_f 

			call getarg(10, leitura)
			buffer = trim(adjustl(leitura))
	 		read(buffer,*) dlambda   			
			write(*,*) dlambda
			
			call getarg(11, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) npts
			write(*,*) npts  

			call getarg(12, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) t_relax
			write(*,*) t_relax  
							
			call getarg(13, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) t_Media
			write(*,*) t_Media  
	
			call getarg(14, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) id_amostra
			write(*,*) id_amostra

			call getarg(15, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) seed_Ensemble
			write(*,*) seed_Ensemble
			
			call getarg(16, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) recuperaDin_int
			write(*,*) recuperaDin_int

			call getarg(17, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) ressorteia_int
			write(*,*) ressorteia_int

			if(recuperaDin_int == 1)then
				recuperaDinamica = .True.
			elseif(recuperaDin_int == 0)then
				recuperaDinamica = .False.
			else
				stop "Sao aceitos apenas os valores 1 ou 0"
			endif
				
			if(ressorteia_int == 1)then
				ressorteia = .True.
			elseif(ressorteia_int == 0)then
				ressorteia = .False.
			else
				stop "Sao aceitos apenas os valores 1 ou 0"
			endif				
			
			kMax2 = (1.0_dp * N) ** (1.0_dp/(max(2.0_dp, gama2 - 1.0_dp)))	
				
			tMax = t_relax + t_Media	
				
	end subroutine


!#######################################################################
!	Obtem parametros para rede PL_UCM com Tubos e
!	condicao RBS em ensemble. A rede PL eh criada primeiro
!	e entao os tubos sao adicionados.
!	Escolhe-se se a ligacao dos tubos eh PLA ou random
!#######################################################################

	subroutine get_Arg_PL_UCM_com_Tubos_RBS_Ensemble_PLA_Random()
		implicit none
		integer :: N0, recuperaDin_int, ressorteia_int
		real(dp) :: omega
		character(len=20) :: leitura, buffer
					
		N0 = iargc()
		
		if(N0 /= 21) stop " Numero de nos N integer, prob_tube real, &
		kMin2 integer, kMax2 integer, gamma2 real > 0, nIn0 real <= 1,&
		 lambda, lambda_i, lambda_f, dlambda, npts, t_relax integer, &
		  t_Media integer, id_amostra integer, seed integer, recuperaDin_int integer, integer sorteia, PLA, plus_xN, f_Tubos, alpha"			
			
			

			call getarg(1, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) N
			write(*,*) 'N ', N

			call getarg(2, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) probTube  			
			write(*,*) 'probtube ', probTube
						
			call getarg(3, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) kMin2
			write(*,*) 'kMin ', kMin2
			
			call getarg(4, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) kMax2
			write(*,*) 'kMax ', kMax2
			
			call getarg(5, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) gama2
			write(*,*) 'gama ', gama2
			
			call getarg(6, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) nInf0  			
			write(*,*) 'nInf0 ', nInf0
			
			call getarg(7, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda 				
			write(*,*) 'lambda ', lambda
			
			call getarg(8, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda_i 				
			write(*,*) 'lambda_i ', lambda_i
			
			call getarg(9, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda_f
			write(*,*) 'lambda_f ', lambda_f 

			call getarg(10, leitura)
			buffer = trim(adjustl(leitura))
	 		read(buffer,*) dlambda   			
			write(*,*) 'dlambda ', dlambda
			
			call getarg(11, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) npts
			write(*,*) 'npts ', npts  

			call getarg(12, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) t_relax
			write(*,*) 't_relax ', t_relax  
							
			call getarg(13, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) t_Media
			write(*,*) 't_Media ', t_Media  
	
			call getarg(14, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) id_amostra
			write(*,*) 'id_amostra ', id_amostra

			call getarg(15, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) seed_Ensemble
			write(*,*) 'Semente ', seed_Ensemble
			
			call getarg(16, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) recuperaDin_int
			write(*,*) recuperaDin_int

			call getarg(17, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) ressorteia_int
			write(*,*) 'ressorteia ', ressorteia_int

			call getarg(18, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) PLA_int
			write(*,*) "PLA ", PLA_int
			
			call getarg(19, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) plus_xN
			write(*,*) "plus_xN ", plus_xN				

			call getarg(20, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) f_Tubos
			write(*,*) "fracao de pontas ", f_Tubos

			call getarg(21, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) alpha
			write(*,*) "expoente PLA ", alpha
			
			if(recuperaDin_int == 1)then
				recuperaDinamica = .True.
			elseif(recuperaDin_int == 0)then
				recuperaDinamica = .False.
			else
				stop "Sao aceitos apenas os valores 1 ou 0"
			endif
				
			if(ressorteia_int == 1)then
				ressorteia = .True.
			elseif(ressorteia_int == 0)then
				ressorteia = .False.
			else
				stop "Sao aceitos apenas os valores 1 ou 0"
			endif				

			if(PLA_int == 1)then
				PLA = .True.
			elseif(PLA_int == 0)then
				PLA = .False.
			else
				stop "Sao aceitos apenas os valores 1 ou 0"
			endif

			kMax2 = (1.0_dp * N) ** (1.0_dp/(max(2.0_dp, gama2 - 1.0_dp)))	
				
			tMax = t_relax + t_Media	
				
	end subroutine



!#######################################################################

	subroutine get_Arg_PL_UCM_com_Folhas_RBS_Ensemble()
		implicit none
		integer :: N0, recuperaDin_int, ressorteia_int
		real(dp) :: omega
		character(len=20) :: leitura, buffer
					
		N0 = iargc()
		
		if(N0 /= 17) stop " Numero de nos N integer, prob_folhas real, &
		kMin2 integer, kMax2 integer, gamma2 real > 0, nIn0 real <= 1,&
		 lambda, lambda_i, lambda_f, dlambda, npts, t_relax integer, &
		  t_Media integer, id_amostra integer, seed integer, recuperaDin_int integer, integer sorteia"			
			
			

			call getarg(1, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) N
			write(*,*) N

			call getarg(2, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) probFolhas  			
			write(*,*) probFolhas
						
			call getarg(3, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) kMin2
			write(*,*) kMin2
			
			call getarg(4, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) kMax2
			write(*,*) kMax2
			
			call getarg(5, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) gama2
			write(*,*) gama2
			
			call getarg(6, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) nInf0  			
			write(*,*) nInf0
			
			call getarg(7, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda 				
			write(*,*) lambda
			
			call getarg(8, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda_i 				
			write(*,*) lambda_i
			
			call getarg(9, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda_f
			write(*,*) lambda_f 

			call getarg(10, leitura)
			buffer = trim(adjustl(leitura))
	 		read(buffer,*) dlambda   			
			write(*,*) dlambda
			
			call getarg(11, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) npts
			write(*,*) npts  

			call getarg(12, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) t_relax
			write(*,*) t_relax  
							
			call getarg(13, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) t_Media
			write(*,*) t_Media  
	
			call getarg(14, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) id_amostra
			write(*,*) id_amostra

			call getarg(15, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) seed_Ensemble
			write(*,*) seed_Ensemble
			
			call getarg(16, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) recuperaDin_int
			write(*,*) recuperaDin_int

			call getarg(17, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) ressorteia_int
			write(*,*) ressorteia_int

			if(recuperaDin_int == 1)then
				recuperaDinamica = .True.
			elseif(recuperaDin_int == 0)then
				recuperaDinamica = .False.
			else
				stop "Sao aceitos apenas os valores 1 ou 0"
			endif
				
			if(ressorteia_int == 1)then
				ressorteia = .True.
			elseif(ressorteia_int == 0)then
				ressorteia = .False.
			else
				stop "Sao aceitos apenas os valores 1 ou 0"
			endif				
			
			kMax2 = (1.0_dp * N) ** (1.0_dp/(max(2.0_dp, gama2 - 1.0_dp)))	
				
			tMax = t_relax + t_Media	
				
	end subroutine
	
	
!#######################################################################

	subroutine get_Arg_PL_UCM_com_Folhas_RBS_Ensemble_PLA_Random()
		implicit none
		integer :: N0, recuperaDin_int, ressorteia_int
		real(dp) :: omega
		character(len=20) :: leitura, buffer
					
		N0 = iargc()
		
		if(N0 /= 20) stop " Numero de nos N integer, prob_folhas real, &
		kMin2 integer, kMax2 integer, gamma2 real > 0, nIn0 real <= 1,&
		 lambda, lambda_i, lambda_f, dlambda, npts, t_relax integer, &
		  t_Media integer, id_amostra integer, seed integer, recuperaDin_int integer, integer sorteia, PLA_Random, plus_xN, alpha"			
			
			

			call getarg(1, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) N
			write(*,*) N

			call getarg(2, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) probFolhas  			
			write(*,*) probFolhas
						
			call getarg(3, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) kMin2
			write(*,*) kMin2
			
			call getarg(4, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) kMax2
			write(*,*) kMax2
			
			call getarg(5, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) gama2
			write(*,*) gama2
			
			call getarg(6, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) nInf0  			
			write(*,*) nInf0
			
			call getarg(7, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda 				
			write(*,*) lambda
			
			call getarg(8, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda_i 				
			write(*,*) lambda_i
			
			call getarg(9, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda_f
			write(*,*) lambda_f 

			call getarg(10, leitura)
			buffer = trim(adjustl(leitura))
	 		read(buffer,*) dlambda   			
			write(*,*) dlambda
			
			call getarg(11, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) npts
			write(*,*) npts  

			call getarg(12, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) t_relax
			write(*,*) t_relax  
							
			call getarg(13, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) t_Media
			write(*,*) t_Media  
	
			call getarg(14, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) id_amostra
			write(*,*) id_amostra

			call getarg(15, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) seed_Ensemble
			write(*,*) seed_Ensemble
			
			call getarg(16, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) recuperaDin_int
			write(*,*) recuperaDin_int

			call getarg(17, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) ressorteia_int
			write(*,*) ressorteia_int
			
			call getarg(18, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) PLA_int
			write(*,*) "PLA ", PLA_int			

			call getarg(19, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) plus_xN
			write(*,*) "plus_xN ", plus_xN	

			call getarg(20, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) alpha
			write(*,*) "expoente PLA ", alpha

			if(recuperaDin_int == 1)then
				recuperaDinamica = .True.
			elseif(recuperaDin_int == 0)then
				recuperaDinamica = .False.
			else
				stop "Sao aceitos apenas os valores 1 ou 0"
			endif
				
			if(ressorteia_int == 1)then
				ressorteia = .True.
			elseif(ressorteia_int == 0)then
				ressorteia = .False.
			else
				stop "Sao aceitos apenas os valores 1 ou 0"
			endif
			
			
			if(PLA_int == 1)then
				PLA = .True.
			elseif(PLA_int == 0)then
				PLA = .False.
			else
				stop "Sao aceitos apenas os valores 1 ou 0"
			endif							
			
			kMax2 = (1.0_dp * N) ** (1.0_dp/(max(2.0_dp, gama2 - 1.0_dp)))	
				
			tMax = t_relax + t_Media	
				
	end subroutine	


!#######################################################################
!	Obtem parametros para rede PL_UCM com Tubos
!	e condicao RBS em ensemble
!#######################################################################

	subroutine get_Arg_PL_UCM_RBS_Ensemble()
		implicit none
		integer :: N0, recuperaDin_int
		character(len=20) :: leitura, buffer
		
		N0 = iargc()
		
		if(N0 /= 15) stop " Numero de nos N integer, &
		kMin2 integer, kMax2 integer, gamma2 real > 0, nIn0 real <= 1,&
		 lambda, lambda_i, lambda_f, dlambda, npts, t_relax integer, &
		  t_Media integer, id_amostra integer, seed integer, recuperaDin_int integer"			
			
			call getarg(1, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) N
			write(*,*) N
						
			call getarg(2, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) kMin2
			write(*,*) kMin2
			
			call getarg(3, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) kMax2
			write(*,*) kMax2
			
			call getarg(4, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) gama2
			write(*,*) gama2
			
			call getarg(5, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) nInf0  			
			write(*,*) nInf0
			
			call getarg(6, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda 				
			write(*,*) lambda
			
			call getarg(7, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda_i 				
			write(*,*) lambda_i
			
			call getarg(8, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) lambda_f
			write(*,*) lambda_f 

			call getarg(9, leitura)
			buffer = trim(adjustl(leitura))
	 		read(buffer,*) dlambda   			
			write(*,*) dlambda
			
			call getarg(10, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) npts
			write(*,*) npts  

			call getarg(11, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) t_relax
			write(*,*) t_relax  
							
			call getarg(12, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) t_Media
			write(*,*) t_Media  
	
			call getarg(13, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) id_amostra
			write(*,*) id_amostra

			call getarg(14, leitura)
			buffer = trim(adjustl(leitura))
			read(buffer,*) seed_Ensemble
			write(*,*) seed_Ensemble
			
			call getarg(15, leitura)
                        buffer = trim(adjustl(leitura))
                        read(buffer,*) recuperaDin_int
                        write(*,*) recuperaDin_int

                        if(recuperaDin_int == 1)then
                                recuperaDinamica = .True.
                        elseif(recuperaDin_int == 0)then
                                recuperaDinamica = .False.
			else
				stop "Sao aceitos apenas os valores 1 ou 0"
            endif


			kMax2 = (1.0_dp * N) ** (1.0_dp/(max(2.0_dp, gama2 - 1.0_dp)))	

			tMax = t_relax + t_Media		
	end subroutine



!####################################################################################
!		Pede parametros ao usuario
!####################################################################################


		subroutine promptUser()
				write(*,*) " Digite a porcentagem de nos infectados inicialmente"
				write(*,*) ""
				read(*,*) nInf0
						
				write(*,*) " Digite a taxa de infeccao da epidemia"
				read(*,*) lambda
				write(*,*) " "

				write(*,*) " Digite o numero maximo de intervalos de tempo para a epidemia"
				write(*,*) " "
				read(*,*) tMax
				write(*,*) " "

				write(*,*) " Qual o numero desejado de amostras?"
				write(*,*) " "
				read(*,*) nSamples
		end subroutine

!##############################################################################################################
!		Getter pra metrica apenas
!##############################################################################################################
		
		subroutine get_Arg_RRN_Mistura()
			implicit none
			integer :: N0
			character(len=20) :: leitura, buffer
			
			N0 = iargc()
			
			if(N0 /= 5) stop "p1 real <=1, k1 integer, k2 integer, id_amostra integer, seed_Ensemble integer"			
				
				call getarg(1, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) p1  			
				write(*,*) p1
				
				call getarg(2, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) k1
				write(*,*) k1
				
				call getarg(3, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) k2
				write(*,*) k2
				
				call getarg(4, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) id_amostra
				write(*,*) id_amostra

				call getarg(5, leitura)
				buffer = trim(adjustl(leitura))
				read(buffer,*) seed_Ensemble
				write(*,*) seed_Ensemble
		end subroutine

end module
