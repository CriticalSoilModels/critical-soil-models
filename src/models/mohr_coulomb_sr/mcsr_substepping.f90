module mod_SRMC_Substepping
    use mod_csm_kinds, only: wp
    use mod_SRMC_funcs, only:  MatVec,  DotProduct_2
    use mod_state_params, only: check4crossing, Update_GK, get_dilation, Check_Unloading
    use mod_strain_invariants, only: calc_eps_inv
    use mod_strain_invar_deriv, only: calc_deps_q_by_deps
    use mod_stress_invariants, only : calc_sig_inv                           
    use mod_yield_function, only: calc_dF_by_dsig, calc_yield_function
    use mod_plastic_potential, only: calc_dg_plas_by_dsig
    use mod_state_params_deriv, only: calc_ddil_by_deps_p, calc_ddil_by_dI

    implicit none
    
contains
    !___________________________________________________________________________________________
   !     ##    ## ######## ##      ## ########  #######  ##    ##
   !     ###   ## ##       ##  ##  ##    ##    ##     ## ###   ##
   !     ####  ## ##       ##  ##  ##    ##    ##     ## ####  ##
   !     ## ## ## ######   ##  ##  ##    ##    ##     ## ## ## ##
   !     ##  #### ##       ##  ##  ##    ##    ##     ## ##  ####
   !     ##   ### ##       ##  ##  ##    ##    ##     ## ##   ###
   !     ##    ## ########  ###  ###     ##     #######  ##    ##
   !########     ###    ########  ##     ##  ######   #######  ##    ##
   !##     ##   ## ##   ##     ## ##     ## ##    ## ##     ## ###   ##
   !##     ##  ##   ##  ##     ## ##     ## ##       ##     ## ####  ##
   !########  ##     ## ########  #########  ######  ##     ## ## ## ##
   !##   ##   ######### ##        ##     ##       ## ##     ## ##  ####
   !##    ##  ##     ## ##        ##     ## ##    ## ##     ## ##   ###
   !##     ## ##     ## ##        ##     ##  ######   #######  ##    ##

   subroutine Newton_Raphson(G, K, eta_y, dilation, G_0, nu, M, M_tc, No, D_min, h, dilationar, &
    G_s, eps_q, I_coeff, I_act, I_0, k_K, k_G, k_D,&
    Gu, Ku, eta_yu, dil_u, dEps, MAXITER, FTOL, &
    F0, Sig_0, alpha)
    !**************************************************************************************
    !  This subroutine determine the elastic proportion with consistent state and elastic *
    !  parameters                                                                         *
    !**************************************************************************************
    implicit none
    !input
    integer, intent(in):: MAXITER
    real(wp), intent(in):: G_0, nu, M, M_tc, No, D_min, h, dilationar, G_s, &
       eps_q, k_G, k_K, k_D, I_act, I_0, FTOL
    real(wp), dimension(6), intent(in):: dEps
    !output
    real(wp), intent(inout):: G, K, eta_y, dilation, I_coeff, &
       Gu, Ku, eta_yu, dil_u, F0, Sig_0(6)
    real(wp), intent(out):: alpha
    !local variables
    logical:: ApplyStrainRateUpdates
    integer:: n, i
    real(wp):: FT, dG, dK, deta, ddil, dI, &
       dI_alpha, I_alpha, &
       p_alpha, q_alpha, dummyvar, F_prime
    real(wp):: dEps_alpha(6), dSig_alpha(6), Sig_alpha(6), &
       n_vec(6), L, dSigdAlpha(6)
    real(wp):: D1, D2, DE(6,6), ddildG(6,6), ddildK(6,6), aux(6,6)

    !Initialize parameters
    FT=1000
    alpha=0.5d0
    n=0
    !Determine Changes in state parameters
    dG=Gu-G
    dK=Ku-K
    deta=eta_yu-eta_y
    ddil=dil_u-dilation
    dI=I_act-I_coeff
    do while (FT >= - FTOL.and.(n<=MAXITER))
       !Store variables
       Gu=G
       Ku=K
       eta_yu=eta_y
       dil_u=dilation
       F0=FT

       n=n+1
       !___________________________________________________________________
       !Compute trial strains and stresses
       dEps_alpha=alpha*dEps
       dI_alpha=alpha*dI
       I_alpha=I_coeff+dI_alpha
       !___________________________________________________________________
       !Evaluate rate crossing
       call check4crossing(I_coeff, I_alpha, dI_alpha,I_0, ApplyStrainRateUpdates)
       if (ApplyStrainRateUpdates) then !Update parameters
          call Update_GK(G_0, nu, I_alpha, I_0, k_G, k_K, Gu, Ku)
          call get_dilation(h, D_min, I_alpha, I_0, eps_q, k_D, ApplyStrainRateUpdates, dil_u)
          eta_yu=M-dil_u*(1.0-No)
       endif
       !___________________________________________________________________
       !Update trial stress
       !Ensemble elastic matrix
       D1  = Ku+(4*Gu/3)
       D2  = Ku-(2*Gu/3)
       DE  = 0.0
       DE(1:3,1:3) = D2
       DE(1,1) = D1
       DE(2,2) = D1
       DE(3,3) = D1
       DE(4,4) = Gu
       DE(5,5) = Gu
       DE(6,6) = Gu
       call MatVec(DE,6,dEps_alpha,6,dSig_alpha)
       Sig_alpha = Sig_0 + dSig_alpha
       !___________________________________________________________
       !Evaluate yield function
       call calc_sig_inv(Sig_alpha, p_alpha, q_alpha, dummyvar)
       call calc_yield_function(q_alpha, p_alpha, eta_yu, FT)

       !Evaluate n=dF/dSig and L=dF/dXs
       call calc_dF_by_dsig(M_tc, eta_yu, Sig_alpha, n_vec)!dF/dSig
       L=-p_alpha*(1.0-No) !L=dF/dXs
       !____________________________________________________________
       !Evaluate F'(alpha)
       !Elastic matrix derivatives
       ddildK=0.0
       ddildK(1:3,1:3) = 1.0
       ddildG=0.0
       ddildG(1:3,1:3)=-2./3.
       ddildG(1,1)=4./3.
       ddildG(2,2)=4./3.
       ddildG(3,3)=4./3.
       ddildG(4,4)=1.0
       ddildG(5,5)=1.0
       ddildG(5,5)=1.0
       aux=alpha* (dK*ddildK+dG*ddildG)
       aux=DE+aux
       call MatVec(aux, 6, dEps, 6, dSigdAlpha)
       F_prime=0.0d0
       do i=1,6 !dot product
          F_prime=F_prime+n_vec(i)*dSigdAlpha(i)
       enddo
       F_prime=F_prime+L*ddil!F'(alpha)
       !______________________________________________________________
       !Update alpha
       alpha=alpha-FT/F_prime
       if ((alpha>1.0d0).or.(alpha<0.0d0)) then
          alpha=0.0d0
       endif
    end do
    !Update variables to return
    G=Gu
    K=Ku
    eta_y=eta_yu
    dilation=dil_u
    I_coeff=I_alpha
    Sig_0=Sig_alpha
 end subroutine Newton_Raphson
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************

!______________________________________________________________________________________
!######## ##     ## ##       ######## ########      ######
!##       ##     ## ##       ##       ##     ##    ##    ##
!##       ##     ## ##       ##       ##     ##    ##
!######   ##     ## ##       ######   ########      ######
!##       ##     ## ##       ##       ##   ##            ##
!##       ##     ## ##       ##       ##    ##     ##    ##
!########  #######  ######## ######## ##     ##     ######
!   ###    ##        ######    #######  ########  #### ######## ##     ## ##     ##
!  ## ##   ##       ##    ##  ##     ## ##     ##  ##     ##    ##     ## ###   ###
! ##   ##  ##       ##        ##     ## ##     ##  ##     ##    ##     ## #### ####
!##     ## ##       ##   #### ##     ## ########   ##     ##    ######### ## ### ##
!######### ##       ##    ##  ##     ## ##   ##    ##     ##    ##     ## ##     ##
!##     ## ##       ##    ##  ##     ## ##    ##   ##     ##    ##     ## ##     ##
!##     ## ########  ######    #######  ##     ## ####    ##    ##     ## ##     ##

 subroutine Euler_Algorithm(G_0, nu, M_tc, M, No,  D_min, h, dilationart, Gs,&
    G, K, eta_y, dilation, &
    erate, I_0, I, dI, k_G, k_K, k_D, dtime, DT, &
    switch_original, &
    Sig, EpsP, dEps, ddil, dEpsp, dSig)
    !************************************************************************
    ! Euler's algorithm (finite difference) integration						*
    ! Returns:																*
    ! Increment of stresses and change of state parameters                  *
    !************************************************************************
    implicit none
    !input
    logical, intent(in):: switch_original
    real(wp), intent(in):: G_0, nu, M_tc, M, No, D_min, h, dilationart, Gs,&
       erate(6), I_0, k_G, k_K, k_D, dtime, DT,&
       dEps(6)
    real(wp), intent(inout):: I, dI
    !output
    real(wp), intent(inout):: G, K, eta_y, dilation, EpsP(6) , Sig(6)
    real(wp), intent(out):: ddil, dSig(6), dEpsp(6)
    !local variables
    logical:: ApplyStrainRateUpdate
    real(wp):: I_f, DE(6,6), D1, D2, dSig_el(6), num, lambda, &
       n_vec(6), p, q, dummyvar, L, b, m_vec(6), a(6), &
       epsq_p, epsv_p, Hard, den, dummyvec(6), Hvp, erate_q

    !________________________________________________________________________
    !Evaluate for strain rate updating
    I_f=I+dI
    call check4crossing(I,  I_f, dI, I_0, ApplyStrainRateUpdate)

    if (ApplyStrainRateUpdate) then
       call Update_GK(G_0, nu, I, I_0, k_G, k_K, G, K)
    endif
    !________________________________________________________________________

    !________________________________________________________________________
    !Compute invariants and derivatives
    call calc_sig_inv(Sig, p, q, dummyvar)
    call calc_eps_inv(EpsP, epsv_p, epsq_p)
    call calc_dF_by_dsig(M_tc, eta_y, Sig, n_vec) !n=dFdSig
    call calc_dg_plas_by_dsig(dilation, Sig, m_vec) !m=dP/dSig
    L = -p*(1.0-No) !L=dF/dXs
    call calc_ddil_by_deps_p(D_min, h, I_0, k_D, epsq_p, epsv_p, &
       EpsP, I, ApplyStrainRateUpdate, a) !a=ddil/dEpsq^p
    call calc_ddil_by_dI(D_min, h, I_0, k_D, epsq_p, I, b) !b=dXs/dI in place of dXs/dErate
    !________________________________________________________________________

    !________________________________________________________________________
    !Compute elastic predictor
    !Ensemble elastic matrix
    D1  = K+(4*G/3)
    D2  = K-(2*G/3)
    DE  = 0.0
    DE(1:3,1:3) = D2
    DE(1,1) = D1
    DE(2,2) = D1
    DE(3,3) = D1
    DE(4,4) = G
    DE(5,5) = G
    DE(6,6) = G

    !Calculate the predictor
    call MatVec(DE, 6, dEps, 6, dSig_el)
    !________________________________________________________________________

    !________________________________________________________________________
    !compute numerator
    call DotProduct_2(dSig_el, n_vec, 6, num)! dSig_el . n_vec

    if ((ApplyStrainRateUpdate).and.(.not.switch_original)) then !Dashpot visc.
       num=num+L*b*dI
    endif
    !_______________________________________________________________________

    !_______________________________________________________________________
    !Compute denominator

    !Compute hardening modulus H=-L.a.m
    call DotProduct_2(a, m_vec, 6, dummyvar)
    Hard=-L*dummyvar !hardening modulus

    !compute Den=n.D.m +H
    call MatVec(DE, 6, m_vec, 6, dummyvec)
    call DotProduct_2(n_vec, dummyvec, 6, den)
    den=den+Hard

    if ((ApplyStrainRateUpdate).and.(switch_original)) then !original
       !compute viscous hardening modulus Hvp=-L.bv.m/(dt*DT)
       call calc_eps_inv(erate,dummyvar,erate_q)
       dummyvec=calc_deps_q_by_deps(erate_q,erate)
       dummyvec=b*dilationart*sqrt(Gs/abs(p))*dummyvec
       call DotProduct_2(dummyvec, m_vec, 6, dummyvar)
       Hvp=-L*dummyvar/(dtime*DT)
       den=den+Hvp
    endif
    !_____________________________________________________________________

    !_____________________________________________________________________
    !compute viscoplastic mult. lambda, dEpsp, dSig, and ddil

    lambda=num/den! viscoelastic multiplier
    dEpsp=lambda*m_vec! Flow rule
    dummyvec=dEps-dEpsp
    call MatVec(DE, 6, dummyvec, 6, dSig) !stress increment
    ddil=0.0
    call DotProduct_2(a, dEpsp, 6, ddil)!plastic hard/softening
    ddil=ddil+b*dI !rate hard/softening
    !_____________________________________________________________________

    !_____________________________________________________________________
    !Update EpsP, Sig, D, and eta_y
    EpsP=EpsP+dEpsp
    Sig=Sig+dSig
    dilation=dilation+ddil
    eta_y=M-dilation*(1.0-No)
    !____________________________________________________________________
 end subroutine Euler_Algorithm
!*************************************************************************


!________________________________________________________________________

! ######  ######## ########  ########  ######   ######
!##    ##    ##    ##     ## ##       ##    ## ##    ##
!##          ##    ##     ## ##       ##       ##
! ######     ##    ########  ######    ######   ######
!      ##    ##    ##   ##   ##             ##       ##
!##    ##    ##    ##    ##  ##       ##    ## ##    ##
! ######     ##    ##     ## ########  ######   ######
!########  ########  #### ######## ########
!##     ## ##     ##  ##  ##          ##
!##     ## ##     ##  ##  ##          ##
!##     ## ########   ##  ######      ##
!##     ## ##   ##    ##  ##          ##
!##     ## ##    ##   ##  ##          ##
!########  ##     ## #### ##          ##

 subroutine Stress_Drift_Correction(G_0, nu, M_tc, M, No, D_min, h, dilationart, Gs, &
    G, K, eta_y, dilation, &
    erate, I_0, I_f, I, dI, k_G, k_K, k_D, dtime, DT, &
    switch_original, MAXITER, F0, FTOL, &
    Sig, EpsP, dEps)
    !************************************************************************
    ! Corrects any drift caused during the Mod. Euler's procedure       	*
    ! Returns corrected stresses and state parameters                       *
    !************************************************************************
    implicit none
    !input
    logical, intent(in):: switch_original
    integer, intent(in):: MAXITER
    real(wp), intent(in):: G_0, nu, M_tc, M, No, D_min, h, dilationart, Gs,&
       erate(6), I_0, k_G, k_K, k_D, dtime, DT,&
       FTOL, dEps(6)
    real(wp), intent(inout):: I, I_f, dI
    !output
    real(wp), intent(inout):: G, K, eta_y, dilation, EpsP(6) , Sig(6), F0
    !local variables
    logical:: ApplyStrainRateUpdate
    integer:: n
    real(wp)::epsq_p, epsv_p, n_vec(6), m_vec(6), L, a(6), b,&
       DE(6,6), D1, D2, Den,  Hard, Hvp, dlambda, &
       dummyvec(6), dummyvar, FC, &
       dil_u, eta_yu, Sigu(6), ddil, dSig(6), dEpsp(6),Epspu(6), &
       p, q, erate_q
    !________________________________________________________________________
    !Evaluate for strain rate updating
    call check4crossing(I, I_f, dI, I_0, ApplyStrainRateUpdate)

    if (ApplyStrainRateUpdate) then
       call Update_GK(G_0, nu, I_f, I_0, k_G, k_K, G, K)
    endif
    !________________________________________________________________________

    !________________________________________________________________________
    !Initialize values
    n=0 !counter
    FC=F0
    !________________________________________________________________________

    do while ((abs(FC)>FTOL).and.(n<MAXITER))
       !___________________________________________________________________
       !update values
       F0=FC
       n=n+1
       !___________________________________________________________________

       !________________________________________________________________________
       !Compute invariants and derivatives
       call calc_sig_inv(Sig, p, q, dummyvar)
       call calc_eps_inv(EpsP, epsv_p, epsq_p)
       call calc_dF_by_dsig(M_tc, eta_y, Sig, n_vec) !n=dFdSig
       call calc_dg_plas_by_dsig(dilation, Sig, m_vec) !m=dP/dSig
       L=-p*(1.0-No) !L=dF/dXs
       call calc_ddil_by_deps_p(D_min, h, I_0, k_D, epsq_p, epsv_p, &
          EpsP, I_f, ApplyStrainRateUpdate, a) !a=ddil/dEpsq^p
       call calc_ddil_by_dI(D_min, h, I_0, k_D, epsq_p, I_f, b) !b=dXs/dI in place of dXs/dErate
       !________________________________________________________________________

       !___________________________________________________________________
       !compute denominator

       !Ensemble elastic matrix
       D1  = K+(4*G/3)
       D2  = K-(2*G/3)
       DE  = 0.0
       DE(1:3,1:3) = D2
       DE(1,1) = D1
       DE(2,2) = D1
       DE(3,3) = D1
       DE(4,4) = G
       DE(5,5) = G
       DE(6,6) = G

       !Compute hardening modulus H=-L.a.m
       call DotProduct_2(a, m_vec, 6, dummyvar)
       Hard=-L*dummyvar !hardening modulus

       !compute Den=n.D.m +H
       call MatVec(DE, 6, m_vec, 6, dummyvec)
       call DotProduct_2(n_vec, dummyvec, 6, den)
       den=den+Hard

       if ((ApplyStrainRateUpdate).and.(switch_original)) then !original
          !compute viscous hardening modulus Hvp=-L.bdI/dErate.m/(dt*DT)
          call calc_eps_inv(erate,dummyvar,erate_q)
          dummyvec=calc_deps_q_by_deps(erate_q,erate)
          dummyvec=b*dilationart*sqrt(Gs/abs(p))*dummyvec
          call DotProduct_2(dummyvec, m_vec, 6, dummyvar)
          Hvp=-L*dummyvar/(dtime*DT)
          den=den+Hvp
       endif
       !___________________________________________________________________

       !___________________________________________________________________
       !Compute DeltaLambda=F0/den and updated Sig, D, and eta_y
       dlambda=F0/den      !error of plastic multiplier
       dEpsp=dlambda*m_vec !error in plastic strain tensor
       Epspu=EpsP-dEpsp
       call MatVec(DE, 6, dEpsp, 6, dSig) !stress drift
       ddil=0.0
       call DotProduct_2(a, dEpsp, 6, ddil)!plastic hard/softening
       !if (switch_original) ddil=ddil+b*dI

       !update
       Sigu=Sig-dSig
       dil_u=dilation-ddil
       eta_yu=M-dil_u*(1.0-No)
       !__________________________________________________________________

       !__________________________________________________________________
       !Evaluate yield function
       call calc_sig_inv(Sigu, p, q, dummyvar)
       call calc_yield_function(q, p, eta_yu, FC)
       !__________________________________________________________________
       !Evaluate change direction
       if (abs(FC)>abs(F0)) then
          call DotProduct_2(n_vec, n_vec, 6, den)
          dlambda=F0/den
          dEpsp=dlambda*m_vec !error in plastic strain tensor
          Epspu=EpsP-dEpsp
          call MatVec(DE, 6, dEpsp, 6, dSig) !stress drift
          Sigu=Sig-dSig
          dil_u=dilation
          eta_yu=eta_y
       end if
       !__________________________________________________________________

       !__________________________________________________________________
       !Update the state parameters
       EpsP=Epspu
       Sig=Sigu
       dilation=dil_u
       eta_y=eta_yu
       !__________________________________________________________________
    end do


 end subroutine Stress_Drift_Correction
!******************************************************************************

end module mod_SRMC_Substepping