!**********************************************************************************************************
!                                                                                                         !      
!                                             UEL for PF-CZM                                              !
!                                        BFGS quasi-newton solver                                         !                                       
!                                                                                                         !
!**********************************************************************************************************
*  Copyright (C) 2019 South China University of Technology, China. All rights reserved.
*  
*  This subroutine implemens the phase field regularized coheisve zone model
*
*  Status: only 2D plane stress and CPS4 elements are considered.
*  
*  Author: J. Y. Wu (jywu@scut.edu.cn) and Y. Huang
*  Date: 31 Oct. 2019
*
*
*  If you want to use this subroutine (research ONLY), please cite our papers:
*  1. Wu, J. Y., 2017. A unified phase-field theory for the mechanics of damage and quasi-brittle failure. 
*     Journal of the Mechanics and Physics of Solids, 103: 72-99.
*  2. Wu, J. Y., 2018. A geometrically regularized gradient damage model with energetic equivalence. 
*     Computer Methods in Applied Mechanics and Engineering, 328: 612-637.
*  3. Wu, J. Y. and Nguyen, V.P., 2018. A length scale insensitive phase-field damage model for brittle 
*     fracture. Journal of the Mechanics and Physics of Solids, 119: 20-42.
*  4. Wu, J. Y., Huang, Y. and Nguyen, V. P., 2019. On the BFGS monolithic algorithm for the unified 
*     phase-field damage theory. Computer Methods in Applied Mechanics and Engineering, 112704.
*  5. Wu, J. Y. and Huang, Y., 2019. Comprehensive ABAQUS implementation of phase-field damage models  
*     for fracture in solids. Theoretical and Applied Fracutre Mechancis, in press.
*      
!**********************************************************************************************************
!
      module NumKind
!
!**********************************************************************************************************
        implicit none
        integer (kind(1)), parameter :: ikind = kind(1), 
     &                                  rkind = kind(0.D0), 
     &                                  lkind = kind(.true.)
!       
      end module Numkind
      

!**********************************************************************************************************
!
      module ModelParam
!
!**********************************************************************************************************
        use NumKind
        implicit none

        ! Constants

        ! Flag of initilization
        logical (lkind) :: bInitialized = .false.

        ! Tolerance
        real    (rkind), parameter :: TOL = 1.0d-12 
        ! number of guass points
        integer (ikind), parameter :: ngp = 3
        
        ! 
        real(rkind) :: thk, EA, mu, l
        real(rkind) :: C(3, 3), D(6, 6)        
        real(rkind) :: gp(2*ngp), gw(ngp)
        real(rkind) :: QQ(18,18) 
        !
        contains

        !===================================
          subroutine Initialize(props, nprops, istype)
        !===================================

            integer (ikind), intent (in) :: nprops, istype
            real    (rkind), intent (in) :: props(nprops)
            integer(ikind) :: indexq(18), i

            !********************************************
            real(rkind) :: G0, K11, K12

            ! material properties
            EA   =  props(1) ! props(1) -- Young's modulus
            mu   =  props(2) ! props(2) -- Poisson's ratio
            l    =  props(3) ! props(5) -- length scale            
            thk  =  props(4) ! props(6) -- thickness
            if (thk < TOL) thk = 1.0
            
            ! elastic stiffness matrix
            G0      = EA / (2.d0*(1.d0 + mu))
            K11     = EA / (1.D0 - mu * mu)
            K12     = mu * K11
            C(:,1) = (/ K11,  K12, 0.D0/)
            C(:,2) = (/ K12,  K11, 0.D0/)
            C(:,3) = (/0.D0, 0.D0,   G0/)

            D(:,1) = (/ l**2*K11, l**2*K12,   0.D0,    0.D0,   0.D0,   0.D0/)
            D(:,2) = (/ l**2*K12, l**2*K11,   0.D0,    0.D0,   0.D0,   0.D0/)
            D(:,3) = (/    0.D0,    0.D0, l**2*G0,    0.D0,   0.D0,   0.D0/)
            D(:,4) = (/    0.D0,    0.D0,   0.D0, l**2*K11, l**2*K12,  0.D0/)
            D(:,5) = (/    0.D0,    0.D0,   0.D0, l**2*K12, l**2*K11,  0.D0/)
            D(:,6) = (/    0.D0,    0.D0,   0.D0,    0.D0,   0.D0, l**2*G0/)
            
            ! integration points
            gp = (/  1.d0/6.d0, 1.d0/6.d0, 2.d0/3.d0, 1.d0/6.d0, 2.d0/3.d0, 1.d0/6.d0 /)
            gw = (/  1.d0, 1.d0, 1.d0 /) / 6.d0
            
            ! dof interchange
            indexq = (/ 1,2,13, 3,4,14, 5,6,15, 7,8,16, 9,10,17, 11,12,18 /)
            ! interchange the locations of dofs
            QQ = 0.d0
            do i = 1, 18
              QQ(indexq(i),i) = 1.d0
            end do             

            bInitialized = .true.
            
            return
          end subroutine Initialize
      !========================================================================= 
      end module ModelParam

!**********************************************************************************************************
!
      module FEM
!
!**********************************************************************************************************
        use NumKind
        implicit none

        contains      
          !==================shape function and its derivative with xi and eta======================    
          subroutine shapefuc(n, dn_xieta, ddn_xieta, nbar, dnbar_xieta, xi, eta)
          
            implicit none      
            real(rkind) :: n(6), dn_xieta(2, 6), ddn_xieta(3, 6), nbar(3), dnbar_xieta(2, 3), xi, eta

            n(1) = 1.d0-3.d0*xi-3.d0*eta+2.d0*xi**2.d0+4.d0*xi*eta+2.d0*eta**2.d0
            n(2) = -xi+2.d0*xi**2.d0
            n(3) = -eta+2.d0*eta**2.d0
            n(4) = 4.d0*xi-4.d0*xi**2.d0-4.d0*xi*eta
            n(5) = 4.d0*xi*eta
            n(6) = 4.d0*eta-4.d0*xi*eta-4.d0*eta**2.d0
            
            dn_xieta(1, 1) = -3.d0+4.d0*xi+4.d0*eta
            dn_xieta(1, 2) = -1.d0+4.d0*xi
            dn_xieta(1, 3) =  0.d0
            dn_xieta(1, 4) =  4.d0-8.d0*xi-4.d0*eta
            dn_xieta(1, 5) =  4.d0*eta
            dn_xieta(1, 6) = -4.d0*eta
            dn_xieta(2, 1) = -3.d0+4.d0*xi+4.d0*eta
            dn_xieta(2, 2) =  0.d0
            dn_xieta(2, 3) = -1.d0+4.d0*eta
            dn_xieta(2, 4) = -4.d0*xi
            dn_xieta(2, 5) =  4.d0*xi
            dn_xieta(2, 6) =  4.d0-8.d0*eta-4.d0*xi

            ddn_xieta(1, 1) = 4.d0
            ddn_xieta(1, 2) = 4.d0
            ddn_xieta(1, 3) = 0.d0
            ddn_xieta(1, 4) = -8.d0
            ddn_xieta(1, 5) = 0.d0
            ddn_xieta(1, 6) = 0.d0
            ddn_xieta(2, 1) = 4.d0
            ddn_xieta(2, 2) = 0.d0
            ddn_xieta(2, 3) = 4.d0
            ddn_xieta(2, 4) = 0.d0
            ddn_xieta(2, 5) = 0.d0
            ddn_xieta(2, 6) = -8.d0
            ddn_xieta(3, 1) = 4.d0
            ddn_xieta(3, 2) = 0.d0
            ddn_xieta(3, 3) = 0.d0
            ddn_xieta(3, 4) = -4.d0
            ddn_xieta(3, 5) = 4.d0
            ddn_xieta(3, 6) = -4.d0

            nbar(1) = 1.d0-xi-eta
            nbar(2) = xi
            nbar(3) = eta
            
            dnbar_xieta(1, 1) = -1.d0
            dnbar_xieta(1, 2) =  1.d0
            dnbar_xieta(1, 3) =  0.d0
            dnbar_xieta(2, 1) = -1.d0
            dnbar_xieta(2, 2) =  0.d0
            dnbar_xieta(2, 3) =  1.d0

            return 
          end subroutine shapefuc

          !===============traditional b matrix==============================================      
          subroutine b_matrix(nu,bu,du,det_jacb,coords,xi,eta)
          
            implicit none
            real(rkind) :: nu(2,6), bu(3,12), du(6,12)
            real(rkind) :: jacb(2,2), inv_jacb(2,2), coords(2, 6), coords1(2, 3)
            real(rkind) :: det_jacb, xi, eta, a1, a2, b1, b2
            
            !local varibles
            real(rkind) :: n(6), dn_xieta(2,6), ddn_xieta(3,6), nbar(3), dnbar_xieta(2, 3)
            real(rkind) :: dn_x(6), dn_y(6), dn_xx(6), dn_yy(6), dn_xy(6)
            integer(ikind) :: i, j
             
            ! shape functions 
            call shapefuc(n, dn_xieta, ddn_xieta, nbar, dnbar_xieta, xi, eta)
            
            coords1(1, :) = (/ coords(1, 1), coords(1, 2), coords(1, 3) /)
            coords1(2, :) = (/ coords(2, 1), coords(2, 2), coords(2, 3) /)

            ! jacob matrix
            jacb = matmul(dnbar_xieta, transpose(coords1))      
         !   jacb = matmul(dn_xieta, transpose(coords))          
            det_jacb = jacb(1,1)*jacb(2,2) - jacb(1,2)*jacb(2,1)
            inv_jacb(1, 1) = jacb(2, 2)
            inv_jacb(1, 2) =-jacb(1, 2)
            inv_jacb(2, 1) =-jacb(2, 1)
            inv_jacb(2, 2) = jacb(1, 1)
            inv_jacb = 1.d0/det_jacb*inv_jacb    

            a1 = inv_jacb(1, 1)
            a2 = inv_jacb(1, 2)
            b1 = inv_jacb(2, 1)   
            b2 = inv_jacb(2, 2)
            
            !initialize varibles
            do i = 1,6
              dn_x(i) = inv_jacb(1,1)*dn_xieta(1,i)
     &                + inv_jacb(1,2)*dn_xieta(2,i)
              dn_y(i) = inv_jacb(2,1)*dn_xieta(1,i)
     &                + inv_jacb(2,2)*dn_xieta(2,i)
            end do

            do i = 1,6
              dn_xx(i) = a1**2*ddn_xieta(1,i) + a2**2*ddn_xieta(2,i) + 2*a1*a2*ddn_xieta(3,i)
              dn_yy(i) = b1**2*ddn_xieta(1,i) + b2**2*ddn_xieta(2,i) + 2*b1*b2*ddn_xieta(3,i)
              dn_xy(i) = a1*b1*ddn_xieta(1,i) + a2*b2*ddn_xieta(2,i) + (a2*b1+a1*b2)*ddn_xieta(3,i)
            end do
            
            ! Nu matrix for displacement
            do j = 1, 6
              nu(1, j) = n(j)
              nu(2, j) = n(j)
            end do

            ! Bu matrix for displacement
            bu = 0.d0
            do j = 1, 6
              bu(1, 2*(j-1) + 1) = dn_x(j)
              bu(2, 2*(j-1) + 2) = dn_y(j)
              bu(3, 2*(j-1) + 1) = dn_y(j)
              bu(3, 2*(j-1) + 2) = dn_x(j)
            end do
            
            ! Du matrix for damage
            du = 0.d0
            do j = 1,6
              du(1, 2*(j-1) + 1) = dn_xx(j)
              du(2, 2*(j-1) + 2) = dn_xy(j)
              du(3, 2*(j-1) + 1) = dn_xy(j)
              du(3, 2*(j-1) + 2) = dn_xx(j)
              du(4, 2*(j-1) + 1) = dn_xy(j)
              du(5, 2*(j-1) + 2) = dn_yy(j)
              du(6, 2*(j-1) + 1) = dn_yy(j)
              du(6, 2*(j-1) + 2) = dn_xy(j)
            end do
          
            return
          end subroutine b_matrix
      
        !********************************************************************
        ! define the dyadic function
          function dyadic(vector1,vector2, vlen)
        !********************************************************************
            integer (ikind) :: vlen, i, j
            real    (rkind) :: vector1(vlen),vector2(vlen)
            real    (rkind) :: dyadic(vlen,vlen)
          
            do i = 1, vlen
              do j = 1, vlen
                dyadic(i,j) = vector1(i) * vector2(j)
              end do
            end do

            return
          end function dyadic

      end module FEM
      

!**********************************************************************************************************
!
      subroutine pfczm(rhs,amatrix,coords,u,svars)
!
!**********************************************************************************************************
        use NumKind
        use ModelParam
        use FEM
        implicit none

        real(rkind):: rhs(18), amatrix(18,18), coords(2,6)
        real(rkind):: svars(4), u(18)
       
        ! local varibles
        real(rkind):: nu(2,6), bu(3,12), du(6,12)
        real(rkind):: uu(12), ru(12), kuu(12,12)
        real(rkind):: rr(18), kk(18,18)
        real(rkind):: det_jacb, dvol
        integer(ikind):: i, j, k
        
        ! initialize varibles
        do i = 1, 6
          uu(2*i - 1) = u(3*i - 2)
          uu(2*i    ) = u(3*i - 1)
        end do

        kuu = 0.d0
        do i = 1, ngp
          do j = 1, ngp      
            call b_matrix(nu,bu,du,det_jacb,coords,gp(i),gp(ngp+i))
            dvol= gw(i)*det_jacb*thk
            kuu = kuu + dvol*(matmul(matmul(transpose(bu),C),bu)+matmul(matmul(transpose(du),D),du))
          end do
        end do   

        write (*, *) 'kuu = ', kuu

        ru = -matmul(kuu,uu) ! applies to hybrid formulation
        rr = 0.d0
        kk = 0.d0

        rr(1:12 ) = ru
        rr(13:15) = 0.d0
          
        kk = 0.d0
        kk(1:12 , 1:12 ) = kuu
        kk(13:15, 13:15) = 0.d0

        do i=16,18
          rr(i) = 0.d0
          kk(i,i) = 1.d0 
        end do
          
        rhs     = matmul(transpose(QQ),rr)
        amatrix = matmul(matmul(transpose(QQ),kk),QQ)
  
        return 
      end subroutine pfczm
     
!**********************************************************************************************************
      subroutine UEL(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars,
     &               props, nprops, coords, mcrd, nnode, 
     &               u, du, v, a, jtype, time, dtime, kstep, 
     &               kinc, jelem, params, ndload, jdltyp, adlmag,
     &               predef, npredf, lflags, mlvarx, ddlmag, mdload,
     &               pnewdt, jprops,njprop,period)
!**********************************************************************************************************

        use NumKind
        use ModelParam
        implicit none

!**********************************************************************************************************
      ! interface of uel, DO NOT change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! variables passed in
        integer (ikind), intent (in    ) :: ndofel, mlvarx, nrhs, 
     &    nsvars, nprops, mcrd, nnode, jtype, kstep, kinc, jelem, 
     &    ndload,npredf, mdload, njprop
     
        integer (ikind), intent (in    ) :: jdltyp(mdload,*), 
     &    lflags(*), jprops(njprop)
     
        real    (rkind), intent (in    ) :: props(nprops), 
     &    coords(mcrd,nnode), u(ndofel), du(mlvarx,*), v(ndofel), 
     &    a(ndofel), time(2), params(3), adlmag(mdload,*),
     &    ddlmag(mdload,*), predef(2,npredf,nnode), dtime, period
     
        ! variables can be updated
        real    (rkind), intent (in out) :: pnewdt
  
        ! variables to be updated (the update of energy(8) is optional)
        real    (rkind), intent (in out) :: rhs(mlvarx,nrhs), 
     &    amatrx(ndofel,ndofel), svars(nsvars), energy(8)
      ! interface of uel, DO NOT change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !********************************************************************************************************
      !                   
      ! user coding to define rhs, amatrx, svars, energy and pnewdt (optional for the last two)
      !
        
        write (*, *) 'TIME = ', time(2)

        ! initialize parameters, etc.                     
        if (.not. bInitialized) call Initialize(props, nprops, jtype)
        
        ! right now only Q4 element is implemented
        call pfczm(rhs(:,1),amatrx,coords,u,svars)
        
        return

      end subroutine uel
