        program read_flux
        implicit none
        INTEGER XNUM_K
        !nombre de niveau
        PARAMETER (XNUM_K=202)
        !flux descendant et montant onde courte et onde longue
        REAL,DIMENSION(XNUM_K)::FLUX_SW_DW,FLUX_SW_UP
        REAL,DIMENSION(XNUM_K)::FLUX_LW_DW,FLUX_LW_UP
        REAL,DIMENSION(XNUM_K)::PZZ !height level flux point
        INTEGER::II
        !Variables thermodynamiques
        REAL, DIMENSION(XNUM_K-1)   :: XPTHT    ! THeta at t
        REAL, DIMENSION(XNUM_K-1)   :: XRHODREF    
        REAL, DIMENSION(XNUM_K-1,2) :: XPRT     ! moist variables at t(kg/kg)
        REAL, DIMENSION(XNUM_K-1) :: XPPABST  ! pressure at t(Pa)
        REAL, DIMENSION(XNUM_K-1,2) :: XPSVT    ! scalar variable )
        REAL, DIMENSION(XNUM_K-1):: XPZHAT ! height level mass point
        !contante pour calcul de Cph
        real  XBOLTZ,XAVOGADRO,XMD,XMV,XRD,XRV,XCPD,XCPV,XCI,XCL
        REAL,DIMENSION(XNUM_K-1)::XCPH
        !g:acceleration gravitationnelle
        real XG
        !lecture fichier
        character*72 filename
        character*72 filename1
        character*72 filename2
        Integer NUMFIC
        !AJOUTER VOS VARIABLE ICI
        REAL :: W_sw                           !l'albédo
        PARAMETER (T=289.)                     !la température de surfrace en Kelvin 
        REAL :: EPS_lw                         !émissivité 
        PARAMETER (bolt=1.380658E-23)          !la constante de boltzmann
        REAL :: flux_net_sw                    !flux net ondes courtes  
        REAL :: flux_net_lw                    !flux net ondes longues 
        REAL :: flux_total                     !flux total


!*************************************************************************
        !CHOIX DES PROFILS ATMOSPHERIQUES ET DES FLUX
 999    write(*,*)'NOM DU FICHIER?'
        write(*,*)'ciel clair -------------------> 1'
        write(*,*)'stratocumulus N=100 cm-3 -----> 2'
        write(*,*)'cirrus prop op ice -----------> 3'
        write(*,*)'stratocumulus N=200 cm-3 -----> 4'
        write(*,*)'cirrus prop op water ---------> 5'
        read(*,*)NUMFIC
        if (NUMFIC==1) then
         filename='FICHIER/ATMO/pro_clearsky.dat'
         filename1='FICHIER/FLUX/flux_cas1.out'
         write(*,*)'vous avez choisi: ',trim(filename)
        elseif (NUMFIC==2) then 
         filename='FICHIER/ATMO/pro_strato100.dat'
         filename1='FICHIER/FLUX/flux_cas2.out'
         write(*,*)'vous avez choisi: ',trim(filename)
        elseif (NUMFIC==3) then 
         filename='FICHIER/ATMO/pro_cirrusice.dat'
         filename1='FICHIER/FLUX/flux_cas3.out'
         write(*,*)'vous avez choisi: ',trim(filename)
        elseif (NUMFIC==4) then 
         filename='FICHIER/ATMO/pro_strato200.dat'
         filename1='FICHIER/FLUX/flux_cas4.out'
         write(*,*)'vous avez choisi: ',trim(filename)
        elseif (NUMFIC==5) then 
         filename='FICHIER/ATMO/pro_cirruswater.dat'
         filename1='FICHIER/FLUX/flux_cas5.out'
         write(*,*)'vous avez choisi: ',trim(filename)
        else 
         write(*,*)'ERREUR CHOIX IMPOSSIBLE'
         goto 999
        endif

!*************************************************************************
        !LECTURE DES FLUX
        open(44,file=filename1)
        do II=1,XNUM_K
        read(44,'(5(E21.15,1x))')PZZ(II),FLUX_SW_DW(II),FLUX_SW_UP(II),&
                       FLUX_LW_DW(II),FLUX_LW_UP(II)
        enddo
!*************************************************************************
        !LECTURE DES PROFILS ATMOSPHERIQUES
        open(30,file=filename)
        DO II=1,XNUM_K-1
        read(30,'(7(2x,1e14.6))')XPZHAT(II),XPTHT(II),XPRT(II,1),XPRT(II,2),&
        XPPABST(II),XPSVT(II,2),XRHODREF(II)
        ENDDO
        close(30)
!*************************************************************************
        !CALCUL DE XCPH
        XBOLTZ      = 1.380658E-23
        XAVOGADRO   = 6.0221367E+23
        XMD    = 28.9644E-3
        XMV    = 18.0153E-3
        XRD    =XBOLTZ*XAVOGADRO   /XMD 
        XRV    =XBOLTZ*XAVOGADRO  /XMV 
        XCPD   = 7.* XRD /2.
        XCPV   = 4.* XRV
        XCL    = 4.218E+3
        XCI    = 2.106E+3
        XG=9.81
        if (NUMFIC==1) then
                        XCPH(:) = XCPD + XCPV  *   XPRT(:,1)
        elseif ((NUMFIC==2).or.(NUMFIC==4)) then 
        XCPH(:) = XCPD + XCPV  *   XPRT(:,1)                             &
                     + XCL   * XPRT(:,2)
        else  
        XCPH(:) = XCPD + XCPV  *   XPRT(:,1)                             &
                     + XCI   * XPRT(:,2)
        ENDIF
!*************************************************************************
        !CODER ICI
         W_sw=FLUX_SW_UP(1)/FLUX_SW_DW(1)                                          !calcul de l'albédo 
         print*,W_sw
         EPS_lw=(FLUX_LW_UP(1)-FLUX_LW_DW(1))/((bolt*(T^4))-FLUX_LW_DW(1))         !calcul de l'émissivité  
         print*,EPS
         flux_net_sw=(FLUX_SW_UP(1)+FLUX_SW_UP(2))-(FLUX_SW_DW(1)+FLUX_SW_DW(2))   !calcul du flux net ondes courtes 
         print*,flux_net_sw
         flux_net_lw=(FLUX_LW_UP(1)+FLUX_LW_UP(2))-(FLUX_LW_DW(1)+FLUX_LW_DW(2))   !calcul du flux net ondes longues 
         print*,flux_net_lw
         flux_total=flux_net_sw+flux_net_lw                                        !calcul du flux total 
         print*,flux_total 

        END
        


        

