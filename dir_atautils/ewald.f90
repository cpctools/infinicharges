!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    InfiniCharges: A program to generate partial charges for periodic systems
!    Copyright (C) 2015  Andrea Gabrieli and Marco Sant                        
! 
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 3 of the License, or
!    (at your option) any later version.
! 
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
! 
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301, USA.
! 
!    See also http://www.fsf.org/licensing/licenses/gpl.html
! 
!    InfiniCharges, including its sources and pointers to the authors
!    can be found at  http://www.physchem.uniss.it/cpc/
!
!    Contact Address:
!    agabrieli@uniss.it
!    msant@uniss.it
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    This module computes the electrostatic potential according to the Ewald
!    formalism.
!
!    It has been developed starting from the corresponding one in
!    Music -- Multipurpose Simulation Code:
!    Object-oriented molecular simulation software
!    Copyright (C) 2002  Amit Gupta, Marty Sanborn, Louis Clark, Shaji 
!    Chempath, Lev Sarkisov, Randy Snurr
!
!    Refer to Frenkel-Smit Chapter 12 ("FS" in comments),
!    Dlpoly 2.20 manual, and Allen-Tildesley for theory and equations of Ewald.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Module ewald
  Use parameters, Only: Pi, RDbl, zero, one, twopi, e2kcal
  Implicit None

  !Product of real space cutoff and sqrt(alpha):
  !afct=dsqrt(alpha)*rcut
  !It is used to compute alpha starting from the real space cutoff.
  !In dlpoly 2.20 manual it should be 3.2 or larger. 3.35 ensures
  !that the resulting alpha will give a well converged
  !real space summation for the given cutof, at the expense of
  !longer reciprocal computation, which is not a problem thanks to "locut".
  Real(Kind=RDbl), Parameter :: afct = 3.35d0

  Type EwaldParams
    Real(kind=RDbl)        :: kappa     !dlpoly_alpha=afct/rcut=sqrt(alpha_FS)
    Real(Kind=RDbl)        :: sqrcut    !Square of real space spherical cutoff
    Integer                :: nk, kmax  !No. of accepted kvecs No. of recip. space images in each dir
    Integer                :: nkmax     !max no. kvecs
    Integer                :: kx_hi, ky_hi, kz_hi !max kx or ky or kz used
    Real(kind=RDbl)        :: locut     !Low cutoff for reciprocal space sum
    !Locut has a strong effect on convergence of results, if changed
    !it should be well tested against a reference having a very small locut, e.g., smaller than 10^(-20).
    Real(Kind=RDbl), Dimension(:,:), Pointer :: s
    Real(Kind=RDbl), Dimension(:,:), Pointer :: si
    Integer, Dimension(:), Pointer :: natoms_per_kind
    Integer, Dimension(:), Pointer :: kind_idx
    Integer, Dimension(:,:), Pointer :: kused
  End Type EwaldParams

  !Ewald parameters variable
  Type(EwaldParams), Target, Save :: eparams

  !Ewald kvectors array for double-sum implementation
  Real(Kind=RDbl), Save :: IGAMMA, IBETA, IALPHA
  
  Real(Kind=RDbl), Dimension(:,:), Allocatable :: us, usi
  Integer, Dimension(:,:), Allocatable :: kus


Contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize the Ewald variables
  ! The real space cutoff is automatically set if rcut <= 0.0.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine ewald_init(kmax, locut, rcut, ell, natoms_per_kind, kind_idx)
    Integer :: i, error
    Real(Kind=RDbl),dimension(3), Intent(In) :: ell
    Real(Kind=RDbl), Intent(In) :: locut, rcut
    Integer, Intent(In) :: kmax
    Integer, Dimension(:), Intent(In) :: natoms_per_kind
    Integer, Dimension(:), Intent(In) :: kind_idx

    !number of k "boxes" in reciprocal space along each direction
    eparams%kmax = kmax
    !max number of storable eikr elements
    eparams%nkmax = (kmax*2+1)**3
    
    !set a default REAL SPACE cutoff to half smallest cell vector
    eparams%sqrcut = (0.5d0*minval(ell))**2
    !real space cutoff, if user specified
    if (rcut.gt.0.0d0) eparams%sqrcut = rcut*rcut
    if (eparams%sqrcut.lt.16.0d0) print *,"WARNING! Are you sure about this cutoff? It seems small."
    !test if kmax is large enough, according to dlpoly 2.20 manual
    if (kmax < (afct*maxval(ell)/dsqrt(eparams%sqrcut))) then
      print *,"ERROR! KMAX= ", kmax, " is too small. At least one dimension is very large."
      print *,"ERROR! KMAX should be bigger than: ", (afct*maxval(ell)/dsqrt(eparams%sqrcut))
      stop
    end if
    eparams%kappa=afct/dsqrt(eparams%sqrcut)

    !cutoff for smallest elements considered in the reciprocal space summation
    eparams%locut = locut

    !allocate space for the natoms_per_kind vector array
    !natoms_per_kind is the number of atoms in each atom-type group
    Allocate(eparams%natoms_per_kind(size(natoms_per_kind)), stat=error)
    If (error /= 0) Then
      Write(0,*) ' Cannot allocate space for eparams%natoms_per_kind'
      Stop              
    End If
    eparams%natoms_per_kind = natoms_per_kind 

    !allocate space for the kind_idx vector array
    !kind_idx: is a vector ordered according to atom-type groups (as in natoms_per_kind),
    !given the index of the atom in the group it returns the index of the atom
    !in the cube file
    Allocate(eparams%kind_idx(sum(natoms_per_kind)), stat=error)
    If (error /= 0) Then
      Write(0,*) ' Cannot allocate space for eparams%kind_idx'
      Stop              
    End If
    eparams%kind_idx = kind_idx

    !allocate space for the kused vector array
    Allocate(eparams%kused(3,eparams%nkmax), stat=error)
    If (error /= 0) Then
      Write(0,*) ' Cannot allocate space for eparams%kused'
      Stop              
    End If
    eparams%kused=-777

    !initialize variables for smallest and largest kx, ky, kz
    eparams%kx_hi=-kmax 
    eparams%ky_hi=-kmax
    eparams%kz_hi=-kmax

    !this vector holds the real part of exp(ik.r) for the atomic structure factor S(k)
    Allocate(eparams%s(eparams%nkmax,Size(natoms_per_kind)),stat=error)
    If (error /= 0) Then
      Write(0,*) ' Cannot allocate space for eparams%s'
      Stop              
    End If
    eparams%s = 0.0_RDbl

    !this vector holds the imaginary part of exp(ik.r) for the atomic structure factor S(k)
    Allocate(eparams%si(eparams%nkmax,Size(natoms_per_kind)),stat=error)
    If (error /= 0) Then
      Write(0,*) ' Cannot allocate space for eparams%si'
      Stop              
    End If
    eparams%si = 0.0_RDbl

    !initialize inverse of of the simulation cell edge lengths
    IALPHA = one/ell(1)
    IBETA  = one/ell(2)
    IGAMMA = one/ell(3)

  End Subroutine ewald_init


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculates the atomic structure factor S(k), independently for each
  ! atom-type group.
  !
  ! It requires the vectors of all ionic atoms (atvecs) and their respective
  ! charges (charges).
  !
  ! The structure factor is used because it converts the Ewald summation from
  ! a double sum over i and j to a single summation over i.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine ewald_sfactor(atvecs,charges)
    Real(Kind=RDbl), Intent(In), Dimension(:,:) :: atvecs
    Real(Kind=RDbl), Intent(In), Dimension(:) :: charges
    Real(Kind=RDbl) :: kalpha, kbeta, kgamma
    Real(Kind=RDbl) :: ikap, ikapPiSq, iVolPi
    Real(Kind=RDbl) :: ka2pi, kb2pi, kg2pi
    Real(Kind=RDbl) :: ksqx, ksqy, ksqz, ksq
    Real(Kind=RDbl) :: factor, sym_factor
    Integer :: kmax, x, y, z, error

    !initialize the number of accepted k vectors
    !the final value will vary according to the atoms coordinates
    !and the kappa and cutoff values, it is updated inside "ewald_sfSum"
    eparams%nk = 0

    !set local kmax
    kmax = eparams%kmax

    !inverse of kappa (units of length)
    ikap = one/eparams%kappa

    !(Pi/kappa)^2 (units of length^2)
    ikapPiSq = (pi*ikap)**2

    !1/Volume * 1/Pi (units of 1/length^3)
    iVolPi = IALPHA*IBETA*IGAMMA/pi

    !loop over the reciprocal space images
    Do x = 0,kmax
      If ( x == 0) then
          sym_factor = 1.0_RDbl
      Else 
          sym_factor = 2.0_RDbl
      End If

      kalpha = Real(x,RDbl)*IALPHA     !units of 1/length
      ka2pi = kalpha*twopi             !units of 1/length
      ksqx = kalpha**2                 !units of 1/length^2

      Do y = -kmax,kmax
        kbeta = Real(y,RDbl)*IBETA
        kb2pi = kbeta*twopi
        ksqy = kbeta**2

        Do z = -kmax,kmax
          !we skip the calculation of the k vectors for the central image
          If (x == 0 .And. y == 0 .And. z == 0) Cycle
          kgamma = Real(z,RDbl)*IGAMMA
          kg2pi = kgamma*twopi
          ksqz = kgamma**2
          ksq = ksqx + ksqy + ksqz    ! units of 1/length^2
          factor = iVolPi*Exp(-ikapPiSq*ksq)/ksq
          !look for too small elements, if they are above "locut" add them to the kvecs
          If (factor > eparams%locut) Then
            Call ewald_sfSum(atvecs,charges,(/ka2pi,kb2pi,kg2pi/), &
                 factor, sym_factor,x,y,z)
          End If
        End Do
      End Do
    End Do

    !trim the relevant arrays
    allocate(us(eparams%nk,size(eparams%natoms_per_kind)), &
            usi(eparams%nk,size(eparams%natoms_per_kind)), &
            kus(3,eparams%nk)) 
    us =eparams%s (1:eparams%nk,:)
    usi=eparams%si(1:eparams%nk,:)
    kus=eparams%kused(:,1:eparams%nk)

    !Deallocate the large arrays
    Deallocate(eparams%kused,stat=error)
    If (error/=0) Then
      Write(*,*) " Could not deallocate eparams%kused"
      Stop
    End If
    Deallocate(eparams%s,stat=error)
    If (error/=0) Then
      Write(*,*) " Could not deallocate eparams%s"
      Stop
    End If
    Deallocate(eparams%si,stat=error)
    If (error/=0) Then
      Write(*,*) " Could not deallocate eparams%si"
      Stop
    End If


  End Subroutine ewald_sfactor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculates the structure factor and store k vectors
  ! for the given reciprocal space k-vectors. 
  ! kvals= 2*pi*n_i/L = k in FS , n_i=-3,-2,-1,0,1,2...., for i=x,y,z
  ! factor = 1/pi/V * exp(-k^2/4/alpha) /( (n_x^2 + n_y^2 +n_z^2)/L^2)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine ewald_sfSum(atvecs, charges, kvals, factor, sym_factor, kx, ky, kz)
    Real(kind=RDbl), Intent(In), Dimension(:,:) :: atvecs !system atom coordinates
    Real(kind=RDbl), Intent(In), Dimension(:) :: charges !system atom charges
    Integer, Intent(In) :: kx, ky, kz
    Real(Kind=RDbl), Dimension(:), Intent(In) :: kvals
    Real(Kind=RDbl) :: smag, rdotk
    Real(Kind=RDbl), Dimension(:), Allocatable :: sum, sumi
    Real(Kind=Rdbl), Intent(InOut) :: factor
    Real(Kind=Rdbl), Intent(In) :: sym_factor
    Integer :: i, j, id, nk, first, last
    Logical :: good


    !calculate the real and imaginary parts of the structure factor
    !for each atom-type group (i), if any of the ith S(-k) pass the 
    !filter k is accepted for all groups (good = .True.) 
    good = .False.
    Allocate(sum(size(eparams%natoms_per_kind)), sumi(size(eparams%natoms_per_kind)))
    sum = zero
    sumi = zero
    !loop over atom-type groups
    Do i = 1, size(eparams%natoms_per_kind)
      if (i == 1) then 
          first = 1
          last = 0
      else
          first = last + 1
      end if
      last = last + eparams%natoms_per_kind(i)
      !print *, "first, last", first, last
      
      !loop over atoms of group i
      Do id = first, last
       j=eparams%kind_idx(id) !index in atvecs of atom id of atom-type group i
       rdotk = dot_product(kvals,atvecs(:,j))  !dimensionless 
       sum(i) = sum(i) + Cos(rdotk)*charges(j) !real part of S(-k) units of |e|
       ! following FS, exp(-ik.r)
       sumi(i) = sumi(i) - Sin(rdotk)*charges(j) !imaginary part of S(-k) Units of |e|
      End Do


      !verify that the value is within our desired bounds. 
      !The 1.0e-3_RDbl value comes from the original MUSIC implementation.
      !In the future it should be further investigated.
      smag = factor**2 * (sum(i)**2 + sumi(i)**2) !units of |e|^2/length^2
      If (smag <= eparams%locut*1.0e-3_RDbl) Cycle
      good = .True.
    End Do
    If (.not. good) Return !no k passed the filter so it is discarded

    !k has survived, increment the number of k vectors calculated by 1
    eparams%nk = eparams%nk + 1
    nk = eparams%nk
    If (nk > eparams%nkmax) Then
      Write(0, *) "Maximum number of k vectors exceed. NKMAX=", eparams%nkmax
      Stop
    End If
    
    !computation of real and imaginary part of S(-k)
    eparams%s(nk,:)  = sum(:)*factor*sym_factor  !units of |e|/length
    eparams%si(nk,:) = sumi(:)*factor*sym_factor !units of |e|/length

    !store the smallest and largest kx,ky,kz accepted. They
    !will be used in the precompute_rdotk subroutine to avoid
    !the computation of more k than necessary
    eparams%kx_hi=max(eparams%kx_hi,abs(kx))
    eparams%ky_hi=max(eparams%ky_hi,abs(ky))
    eparams%kz_hi=max(eparams%kz_hi,abs(kz))
    !map kx, ky, kz to this nk, again to speedup
    !precompute_rdotk
    eparams%kused(:,nk) = (/kx, ky, kz/)
    
    Deallocate(sum, sumi)

  End Subroutine ewald_sfSum


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculates crdotk and srdotk for given grid point for each k used during
  ! the computation of S(-k) in "ewald_sfactor". Finally, compute the reciprocal
  ! space potential by dot product.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine ewald_get_reciprocal(r,tpiell,rec_pot,q)
    Real(Kind=RDbl), Dimension(:), Intent(In) :: tpiell, r
    Real(Kind=RDbl), Intent(In) :: q
    Real(rdbl), Dimension(:), Intent(Out) :: rec_pot
    Integer :: i, kx, ky, kz
    Real(rdbl) :: kxrx, kyry, kzrz
    Real(Kind=RDbl), Dimension(eparams%nk) :: crdotk, srdotk
    Double Complex, Dimension(0:eparams%kx_hi) :: eikx
    Double Complex, Dimension(-eparams%ky_hi:eparams%ky_hi) :: eiky
    Double Complex, Dimension(-eparams%kz_hi:eparams%kz_hi) :: eikz
    Double Complex :: eikr

    !Construct exp(ik . r) for k vectors up to the maximum number
    !used for the computation of S(-k) in "ewald_sfactor"

    !e^0 = 1, for kx=ky=kz=0
    eikx(0) = (1.0, 0.0)
    eiky(0) = (1.0, 0.0)
    eikz(0) = (1.0, 0.0)

    !eikx, eiky, eikz, for kx=ky=kz=1
    kxrx=tpiell(1)*(r(1))
    eikx(1) = Cmplx ( Cos(kxrx) , Sin(kxrx ) )! e(i kx. rx)

    kyry=tpiell(2)*(r(2))
    eiky(1) = Cmplx ( Cos(kyry) , Sin(kyry ) )! e(i ky. ry)

    kzrz=tpiell(3)*(r(3))
    eikz(1) = Cmplx ( Cos(kzrz) , Sin(kzrz ) )! e(i kz. rz)

    !set exp(-ik.r), there is no need to do this for eikx
    !since we loop only over positive kx
    eiky(-1) = Conjg (eiky(1))
    eikz(-1) = Conjg (eikz(1))

    !Calculate remaining kx, ky and kz by recurrence
    !Message from THE MAN: "This looks cool to me (2-28-91)"
    !Is the above quote From Theodorou to a hapless student?

    Do kx = 2,eparams%kx_hi
      eikx( kx) = eikx(kx-1) * eikx(1)
    End Do

    Do ky = 2,eparams%ky_hi
      eiky( ky) = eiky(ky-1) * eiky(1)
      eiky(-ky) = Conjg (eiky(ky))
    End Do

    Do kz = 2,eparams%kz_hi
      eikz( kz) = eikz(kz-1) * eikz(1)
      eikz(-kz) = Conjg (eikz(kz))
    End Do

    !get crdotk and srdotk
    Do i = 1, eparams%nk
      eikr = eikx(kus(1,i))*eiky(kus(2,i))*eikz(kus(3,i))
      crdotk(i) = Real(eikr)
      srdotk(i) = Aimag(eikr)
    End Do

    !compute the reciprocal space potential
    rec_pot = 0.0d0
    do i=1,size(rec_pot)
     rec_pot(i) = rec_pot(i) + dot_product(crdotk(:), us(:,i))  - &
                               dot_product(srdotk(:),usi(:,i))
    end do
    rec_pot = rec_pot * (q * e2kcal)

  End Subroutine ewald_get_reciprocal

 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculates the real space part of the potential at atvec1 with q=charge1 due
  ! to the presence of charges2 at atvecs2. 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine ewald_get_real(ell, ellh2, atvec1, atvecs2, &
      charge1, charges2, pot, at)
    Real(Kind=RDbl), Dimension(:,:), Intent(In) :: atvecs2
    Real(Kind=RDbl), Intent(In) :: charge1
    Real(Kind=RDbl), Dimension(:), Intent(In) :: atvec1, charges2
    Real(Kind=RDbl), Intent(Out) :: pot
    Real(Kind=RDbl) :: rpot

    Real(Kind=RDbl), Allocatable, Dimension(:) :: di, ch
    Integer :: natoms2
    Integer :: i, j, cnt
    Integer, Intent(In) :: at
    Real(Kind=RDbl) :: dist2
    Real(Kind=RDbl), Dimension(3) :: dat2, d

    !variables for minimum image
    Real(Kind=RDbl), Dimension(3), Intent(In) :: ell, ellh2

    !zero the potential
    pot = zero
    rpot = zero


    !get the number of atoms in atom-type group
    natoms2 = Size(atvecs2,2)

    !get the separation distance between the two sets of coords,
    !we use the minimum image convention.
    allocate(di(natoms2), ch(natoms2))
    cnt=0
    !loop over atoms in given atom-type group
    Do i=1,  natoms2
        d = atvec1 - atvecs2(:,i)
        dat2=d**2
        do j = 1, 3
            if (dat2(j) > ellh2(j)) then
                if (d(j) > 0.0_RDbl) then
                    d(j) = d(j) - ell(j)
                else 
                    d(j) = d(j) + ell(j)
                end if
                dat2(j) = d(j)**2
                if (dat2(j) > ellh2(j)) then
                    write(*,*) "ERROR: pbc problem along direction ", j
                    stop
                end if
            end if
        end do
       
        dist2=sum(dat2)
        if (dist2 < eparams%sqrcut) then
         cnt=cnt+1
         !store values to be used in erfc only for the accepted points
         di(cnt) = sqrt(dist2)*eparams%kappa
         ch(cnt) = charges2(i)
        end if

    End Do


    !compute the real space sum
    rpot = get_real_ewald(di(1:cnt), ch(1:cnt), charge1, eparams%kappa)
    
    !the potential returned is in kcal/mol
    pot = rpot*e2kcal

  End Subroutine ewald_get_real

  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine computes the real space part of Ewald by using
  ! the approximation to the complementary error function taken from:
  ! Abramowitz and Stegun, Handbook of Mathematical Functions,
  ! National Bureau of Standards, formula 7.1.28        
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Function get_real_ewald(x1,ch,charge1,kappa)
    Real(Kind=RDbl), Dimension(:), Intent(In) :: x1, ch
    Real(Kind=RDbl), Intent(In) :: charge1, kappa
    Real(Kind=RDbl)  :: get_real_ewald
    Real(Kind=RDbl), Parameter :: b1 = 0.0705230784, b2 = 0.0422820123, &
                     b3 = 0.0092705272, b4 = 0.0001520143, b5 = 0.0002765672, b6 =  0.0000430638

    get_real_ewald = sum(charge1*kappa*ch/((1.0d0+ &
                            x1*b1+ &
                            x1*x1*b2+ &
                            x1*x1*x1*b3+ &
                            x1*x1*x1*x1*b4+ &
                            x1*x1*x1*x1*x1*b5+ &
                            x1*x1*x1*x1*x1*x1*b6)**16*x1))
    Return
  End Function get_real_ewald


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Print Ewald parameters used throughout the module
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine print_ewald_info(ou)
    Integer, Intent(In) :: ou
    Integer :: fu
    Logical :: isopen
    Inquire(Unit=ou, Opened=isopen)
    If (isopen) Then
        fu = ou
    Else
        fu = 5
        Write(5,*) "File unit ", ou, " is not open, writing on standard output."
    End If
    write(fu,*), "Ewald info..."
    write(fu,*), "kappa     ", eparams%kappa
    write(fu,*), "locut     ", eparams%locut
    write(fu,*), "sqrcut    ", eparams%sqrcut
    write(fu,*), "rcut      ", dsqrt(eparams%sqrcut)
    write(fu,*), "kmax      ", eparams%kmax
    write(fu,*), "nk        ", eparams%nk
    write(fu,*), "nkmax     ", eparams%nkmax
    write(fu,*)  "kx_hi     ", eparams%kx_hi
    write(fu,*)  "ky_hi     ", eparams%ky_hi
    write(fu,*)  "kz_hi     ", eparams%kz_hi

  End Subroutine
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Cleans up any allocated memory
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine ewald_cleanup()
    Integer :: error

    Deallocate(eparams%natoms_per_kind,stat=error)
    If (error/=0) Then
      Write(*,*) " Could not deallocate eparams%natoms_per_kind"
      Stop
    End If
    Deallocate(eparams%kind_idx,stat=error)
    If (error/=0) Then
      Write(*,*) " Could not deallocate eparams%kind_idx"
      Stop
    End If

    Deallocate(us,usi,kus,stat=error)
    If (error/=0) Then
      Write(*,*) " Could not deallocate auxiliary sf arrays"
      Stop
    End If

  End Subroutine ewald_cleanup


End Module ewald
