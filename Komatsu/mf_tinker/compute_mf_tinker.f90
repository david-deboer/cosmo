PROGRAM Compute_Mf_Tinker
  USE linearpk
  USE sigma
  ! A sample program for computing dn/dlnRh and dn/dlnMh using 
  ! Tinker et al.'s formula.
  ! Units: 
  ! - Rh in h^-1 Mpc, 
  ! - Mh in omega_matter h^-1 M_solar, 
  ! - dndlnRh in h^3 Mpc^-3, and
  ! - dndlnMh in h^3 Mpc^-3
  ! August 25, 2008: E.Komatsu
  IMPLICIT none
  double precision :: deltac=1.6865d0,mf,dndlnRh,dndlnMh
  double precision :: lnnu,dlnnudlnRh,Mh
  real :: chebev
  real :: lnsigma2,lnRh,Rh,dlnsigma2dlnRh
  character(len=128) :: filename
  integer :: n
! read in and tabulate P(k)
  filename='wmap5baosn_max_likelihood_matterpower.dat'
  n=896 ! # of lines in the file
  CALL open_linearpk(filename,n)
! fit sigma^2(R) to Chebyshev polynomials
  CALL compute_sigma2
  CALL close_linearpk
! now output dn/dlnRh [in h^3 Mpc^-3] and dn/dlnMh [h^3 Mpc^-3] 
! as a function of R [h^-1 Mpc] and M [omega_matter h^-1 M_solar]...'
  open(1,file='Rh_dndlnRh.txt')
  open(2,file='Mhom0_dndlnMh.txt')
  do lnRh=lnR1,lnR2,0.01
     lnsigma2 = CHEBEV(lnR1,lnR2,c,ndim,lnRh)          ! ln(sigma^2)
     dlnsigma2dlnRh = CHEBEV(lnR1,lnR2,cder,ndim,lnRh) ! dln(sigma^2)/dlnRh
     lnnu = 2d0*dlog(deltac)-dble(lnsigma2) ! ln(nu)
     dlnnudlnRh = -dble(dlnsigma2dlnRh)     ! dln(nu)/dlnRh
     Rh=exp(lnRh) ! h^-1 Mpc
     dndlnRh = (3d0/4d0/3.1415926535d0)*dlnnudlnRh*mf(lnnu)/dble(Rh)**3d0
     write(1,'(2E13.5)') Rh,dndlnRh
     dndlnMh = dndlnRh/3d0 ! in units of h^3 Mpc^-3
     Mh = (4d0*3.1415926535d0/3d0)*2.775d11*Rh**3d0 ! in units of omega_matter h^-1 M_solar
     write(2,'(2E13.5)') Mh,dndlnMh
  enddo
  close(1)
  close(2)
END PROGRAM Compute_Mf_Tinker
