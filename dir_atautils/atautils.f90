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
!    This module contains a series of subroutines related to the generation of
!    ESP data for the InfiniCharges program.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module atautils
    use cubecruncher, only: cube_type, init_cube, deallocate_cube, read_cube
    use parameters, only: RDbl, twopi, autokcal, bohrtoangs, dtoea
    implicit none
    private rm_mask, a, b, ou, first_ata_called, ierr, current_radius, cube, cucoords
    ! pointer must be hidden from f2py, use private.
    public :: init, filter, get_ata_matrix, get_ata_matrix_post, release, tdf_smooth, &
        get_zatom, autokcal, read_xyz
    !f2py real(RDbl), parameter :: autokcal
    
    !mask for the grid points, identifying the ones to be kept
    integer, allocatable, dimension(:,:,:) :: rm_mask
    !For linear system of equations Aq=B
    real(RDbl), allocatable, dimension(:,:) :: a !model data A
    real(RDbl), allocatable, dimension(:) :: b !reference data B
    ! this is the unit where output from this module will be written
    integer, parameter :: ou=50
    logical :: first_ata_called = .False.
    integer :: ierr
    real(RDbl) :: current_radius !vdw scaling factor gamma actually used

    type(cube_type), pointer :: cube !relevant data from cube
    !double precision atom coordinates from cube%coords
    real(RDbl), allocatable, dimension(:,:) :: cucoords
    
    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read and store all data from cube file.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine init(cube_name_in)
        integer iunit_cube_in
        character(len=80) cube_name_in
        !f2py intent in cube_name_in
        logical op

        !opening the output file for the entire module
        inquire(unit=ou, opened=op)
        if (.not. op) open(unit=ou, file="atautils.log")

        iunit_cube_in = 701
        allocate(cube)
        call init_cube(cube)

        write(ou,fmt='(A)',advance="No") "Reading cube..."
        open(iunit_cube_in,file=trim(cube_name_in))
        !read relevant data
        call read_cube(cube,iunit_cube_in)
        close(iunit_cube_in)
        
        !store coordinates in double precision
        allocate(cucoords(size(cube%coords,1),size(cube%coords,2)),stat=ierr)
        if (ierr /= 0) then
          write(ou,*), "Could not allocate cucoords."
          stop
        end if
        cucoords= cube%coords

        write(ou,*) "Done."

        ! print useful dimensions
        write(ou,*), "bin size (Bohr): "
        write(ou,'(3f12.6)'), cube%dh(1,:)
        write(ou,'(3f12.6)'), cube%dh(2,:)
        write(ou,'(3f12.6)'), cube%dh(3,:)
        write(ou,*), "box size (Bohr): "
        write(ou,'(3f12.6)'), cube%h(1,:)
        write(ou,'(3f12.6)'), cube%h(2,:)
        write(ou,'(3f12.6)'), cube%h(3,:)
        write(ou,*), "box size (A): "
        write(ou,'(3f12.6)'), cube%h(1,:)*bohrtoangs
        write(ou,'(3f12.6)'), cube%h(2,:)*bohrtoangs
        write(ou,'(3f12.6)'), cube%h(3,:)*bohrtoangs
        write(ou,*), ""
    
    end subroutine
   

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Remove all grid points lying inside a sphere 
    ! surrounding each atom and having radius = rvdW*atomic_radius.
    ! The matrix rm_mask contains the index of the grid points not excluded,
    ! and -1 for the excluded ones.
    ! This subroutine can be called more than once ONLY for increasing rvdW (gamma). 
    ! In this case points already excluded are not recomputed.
    ! n_good_grid_pts stores the number of included points.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine filter(rvdW, print_idx, n_good_grid_pts, e2_tot)
        use parameters, only: atm_radi
        logical inside
        integer :: icnt, scnt, gcnt, xx, yy, zz, i, j, error, temp
        real(kind=8), allocatable, dimension(:) :: rcut, rsqcut
        real(kind=8) :: cx, cy, cz, rr
        real(kind=8), dimension(3) :: cubeh2, cubev, rrv, rrv2, cc, ctemp
        real(kind=8), allocatable, dimension(:,:) :: coords
        integer, allocatable, dimension(:) :: Zatom
        real(kind=8), intent(in) :: rvdW
        !cf2py intent in rvdW
        logical, intent(in) :: print_idx
        !cf2py intent in print_idx  
        integer, intent(out) :: n_good_grid_pts
        !cf2py intent out n_good_grid_pts
        real(kind=8), intent(out) :: e2_tot
        !cf2py intent out e2_tot

        write(ou,*) "Working for gamma: ", rvdW

        ! initialize variables for PBC
        do i = 1,3
          cubeh2(i)=(cube%h(i,i)*0.5d0)**2
          cubev(i) =(cube%h(i,i))
        enddo
        
        ! fold atoms to fall between 0 and L
        do i = 1, cube%Natom
          cucoords(:,i) = cucoords(:,i) - cubev * anint(cucoords(:,i)/cubev-0.5d0)
        end do

        ! sort the atoms based on their atomic number
        ! we expect small z atoms to be on the surface of the system
        ! thus may speedup the filtering procedure if scanned first
        allocate(coords(size(cucoords,1),size(cucoords,2)),stat=error)
        allocate(Zatom(size(cube%Zatom)),stat=error)
        Zatom = cube%Zatom
        coords= cucoords

        do i = 2, cube%Natom
          j = i - 1
          temp = Zatom(i)
          ctemp = coords(:,i)
          do
            if (j < 1) exit
            if (Zatom(j)<=temp) exit
            Zatom(j+1)    = Zatom(j)
            coords(:,j+1) = coords(:,j)
            j = j - 1
          end do
          Zatom(j+1) = temp
          coords(:,j+1) = ctemp
        end do
       
        ! initialize mask
        write(ou,*), "Generating mask for points to be removed"
        if ( .not. allocated(rm_mask) ) then
          allocate(rm_mask(cube%npoints(1), cube%npoints(2), cube%npoints(3)), stat=ierr)
          if ( ierr /= 0 ) then 
            write(ou,*), "Could not allocate rm_mask"
            stop
          end if
          rm_mask = -999

          ! set the initial radius
          current_radius = rvdW
          write(ou,*) "Generating the inital mask."
        end if
        n_good_grid_pts = 0
        ! store the current vdW radius, the module works for increasing radii only
        if (current_radius <= rvdW) then
         current_radius = rvdW 
        else
          write(ou,*) "Error . The module works for increasing radii only"
          stop
        end if
        
        !initialize filter
        allocate(rcut(size(atm_radi)),rsqcut(size(atm_radi)))
        
        ! points lying within rcut (in bohr) from an atom center are skipped
        ! one value for each atomic number
        rcut=rvdW*atm_radi
        rsqcut=rcut*rcut
        
        ! counter of processed points
        icnt=0
        ! counter of removed points
        scnt=0
        ! counter of good points
        gcnt=0
        
        if (print_idx) then
          open(111,file="keep_idx."//trim(adjustl("1")))
        end if

        ! loop over grid points
        do xx=1,cube%npoints(1)
         do yy=1,cube%npoints(2)
          do zz=1,cube%npoints(3)
           icnt=icnt+1
           ! check if the xx,yy,zz has been already excluded by a previous run
           ! with smaller radius
           if (rm_mask(xx,yy,zz) >= 0 .or. rm_mask(xx,yy,zz) == -999) then
            cx = (xx-1)*cube%dh(1,1)
            cy = (yy-1)*cube%dh(2,2)
            cz = (zz-1)*cube%dh(3,3)
            cc = (/ cx, cy, cz /)
            inside = .False.
            ! loop over atoms
            do j=1,cube%Natom
             ! this works ONLY for orthorhombic simulation box
             rrv=cc-coords(:,j)
             rrv2=rrv*rrv
             do i = 1, 3
               if (rrv2(i).gt.cubeh2(i)) then
                 if (rrv(i).gt.0.0d0) then
                   rrv(i)=rrv(i)-cubev(i)
                 else
                   rrv(i)=rrv(i)+cubev(i)
                 end if
                 rrv2(i)=rrv(i)*rrv(i)
                 if (rrv2(i).gt.cubeh2(i)) then
                   write(*,*) i, " pbc problem"
                   write(*,*) cubev(i)
                   write(*,*) cc, cc-coords(:,j)
                   write(*,*) coords(i,j)
                   write(*,*) rrv(i)
                   stop
                 end if
               end if
             end do
             
             rr=sum(rrv2)

             ! rsqcut is a vector having the squared cutoff for each atom kind
             if ( rr < rsqcut(Zatom(j)) ) then
              scnt=scnt+1
              inside = .True.
              rm_mask(xx,yy,zz) = -1
              exit
             end if
            end do
           else
            ! the point is excluded
            inside = .True.
            scnt=scnt+1
           end if
           if (.not. inside) then
              ! the point is valid, update the mask
              gcnt = gcnt + 1 
              if (rm_mask(xx,yy,zz) == -999) rm_mask(xx,yy,zz) = gcnt
              if (print_idx) write(111,"(i9)") icnt
           endif
          end do
         end do
        end do
        
        if (print_idx) then
          close(111)
        end if

        n_good_grid_pts = icnt - scnt
        if (n_good_grid_pts /= gcnt) then
          write(ou,*) "Error: wrong number of good grid points.", n_good_grid_pts, gcnt
          stop
        end if
        write(ou,*), "number of elements processed               = ", icnt
        write(ou,*), "number of spatial elements removed      = ", scnt
        write(ou,*), "number of useful elements      = ", n_good_grid_pts

        ! compute the sum of the squared ESP for all grid pts
        ! for future check of cube files used
        e2_tot=sum(cube%grid**2)
        write(ou,*) "Total sum of ESP**2: ", e2_tot
        
    end subroutine 
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate arrays when done
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine release()
        use ewald, only: ewald_cleanup
        call deallocate_cube(cube)
        deallocate(cucoords, stat=ierr)
        if ( ierr /= 0 ) then
          write(ou,*), "Could not release cucoords"
          stop
        end if
        if (allocated(rm_mask).or.allocated(a).or.allocated(b)) then
          deallocate(rm_mask, a, b, stat=ierr)
          if ( ierr /= 0 ) then
              write(ou,*), "Could not release."
              stop
          end if
          call ewald_cleanup()
          write(ou,*), "rm_mask, a, and b deleted and ewald cleaned."
        end if
        first_ata_called = .False.
        write(ou,*), "Cube deleted."
        write(ou,*), ""
    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Return the ATA, ATB, BTB and related averages.
    ! For the first vdW scaling factor (gamma) the computation is performed via Ewald
    ! summation method, afterwards the ESP data already computed are simply filtered to
    ! obtain the values for the larger gamma.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_ata_matrix(kmax, locut, rcut, n_good_grid_pts, natoms_per_kind, nkinds, kind_idx, natoms, &
               print_pot, esp_sign, atb, ata, btb, avg_a)
        use ewald, only: ewald_init, ewald_sfactor, print_ewald_info, ewald_get_reciprocal,  ewald_get_real
        real(kind=RDbl),dimension(3) :: ell, iell, tpiell, ellh2
        real(kind=RDbl), dimension(cube%Natom) :: charges
        real(kind=RDbl),dimension(3) :: atvec1 !grid point coordinates
        real(kind=RDbl),dimension(3, cube%Natom) :: atvecs2 !atoms coordinates
        real(kind=RDbl) :: pot, cx, cy, cz, t1, t2, test_charge
        integer :: i, j, xx, yy, zz, icnt, at, ccnt, first, last, error
        
        !WARNING: to use RDbl the file .f2py_f2cmap is necessary, otherwise the code 
        !will not work.
        integer, intent(in) :: kmax
        !cf2py intent in kmax
        real(kind=RDbl), intent(in) :: locut, rcut
        !cf2py intent in locut, rcut
        integer, intent(in) :: n_good_grid_pts
        !cf2py intent in n_good_grid_pts
        integer, dimension(nkinds), intent(in) :: natoms_per_kind
        !cf2py intent in natoms_per_kind
        integer, intent(in) :: nkinds
        !cf2py intent in nkinds
        integer, dimension(natoms), intent(in) :: kind_idx
        !cf2py intent in kind_idx
        integer, intent(in) :: natoms
        !cf2py intent in natoms
        logical, intent(in) :: print_pot
        !cf2py intent in print_pot
        real(kind=RDbl), intent(in) :: esp_sign
        !cf2py intent in esp_sign
        real(kind=RDbl), dimension(nkinds), intent(out) :: atb
        !cf2py intent out atb
        real(kind=RDbl), dimension(nkinds,nkinds), intent(out) :: ata
        !cf2py intent out ata
        real(kind=RDbl), intent(out) :: btb
        !cf2py intent out btb
        real(kind=RDbl), dimension(nkinds), intent(out) :: avg_a
        !cf2py intent out avg_a
        
        real(kind=RDbl) :: avg_b
        
        real(kind=RDbl), dimension(nkinds) :: rec_pot

        !model ESP should be computed via Ewald only for the first rvdW (gamma)
        if (.not. first_ata_called) then

          if (.not. allocated(a) .and. .not. allocated(b)) then
           allocate(a(n_good_grid_pts,nkinds), b(n_good_grid_pts), stat=error)
           if (error /= 0) then
               write(ou,*), "ERROR: could not allocate a and b matrices."
               stop
           end if
           first_ata_called= .True. 
          else 
           write(ou,*), "ERROR: this function can be called only once for a given cube file."
           stop
          end if
           
          !check number of atoms for consistency
          if (sum(natoms_per_kind) /= cube%Natom) then
              write(ou,*), "Wrong total number of atoms in natoms_per_kind"
              stop
          endif
          if (size(kind_idx) /= cube%Natom) then
              write(ou,*), "Wrong total number of atoms in kind_idx"
              stop
          endif

          !all atoms' charges are set to 1.0 according to the method
          charges = 1.0_RDbl
          !to perform tests it is possible to use any set of charges.
          !include 'charges.dat'

          !test charge value
          test_charge = 1.0_RDbl
          
          !convert atoms coordinates
          do i = 1, cube%Natom
            do j = 1, 3
              atvecs2(j,i) = cucoords(j,i)*bohrtoangs
            end do
          end do
          
          !convert simulation box cell and precompute useful quantities
          ell(1) = cube%h(1,1)*bohrtoangs
          ell(2) = cube%h(2,2)*bohrtoangs
          ell(3) = cube%h(3,3)*bohrtoangs
          iell = 1.0_RDbl/ell
          tpiell = twopi*iell

          !square of half box
          ellh2 = (ell*0.5d0)**2

          !allocate relevant arrays and set Ewald alpha
          call ewald_init(kmax, locut, rcut, ell, natoms_per_kind, kind_idx)

          !compute only once the eikr for each atom-type group
          call ewald_sfactor(atvecs2, charges)
          call print_ewald_info(ou)
          
          !main loop for Ewald computation
          call cpu_time(t1)
          ccnt = 0 !total grid points counter
          icnt = 0 !good grid points counter
          a = 0.0_RDbl !A matrix
          b = 0.0_RDbl !B matrix
          avg_a = 0.0_RDbl !average of each column of A
          avg_b = 0.0_RDbl !average of B
          do xx=1,cube%npoints(1)
           do yy=1,cube%npoints(2)
            do zz=1,cube%npoints(3)
             ccnt= ccnt+1
             if (rm_mask(xx,yy,zz) >= 0) then
              !good grid point
              icnt=icnt+1
              if (mod(icnt,100000) == 0) then
               write(ou,'(a,i9,a,f6.2,a)') "Map progress: ",icnt," points calculated", &
                    (icnt*100.0)/(n_good_grid_pts),"% of total"
              end if

              !actual grid point coordinates
              cx = (xx-1)*cube%dh(1,1)*bohrtoangs
              cy = (yy-1)*cube%dh(2,2)*bohrtoangs
              cz = (zz-1)*cube%dh(3,3)*bohrtoangs
              atvec1 = (/cx,cy,cz/)

              !compute the reciprocal space part of Ewald summation, accumulated into "pot"
              call ewald_get_reciprocal(atvec1,tpiell, rec_pot, test_charge)
              a(icnt,:) = rec_pot(:)

              !loop over atom-type groups
              do at = 1, size(natoms_per_kind)
               !set indexes of first and last atom in this group
               if (at == 1) then
                   first = 1
                   last = 0
               else
                   first = last + 1
               end if
               last = last + natoms_per_kind(at)
              
               !compute the real part of Ewald and give final potential "pot"
               call ewald_get_real(ell, ellh2, atvec1, atvecs2(:, kind_idx(first:last)), &
                                            test_charge, charges(kind_idx(first:last)), pot, &
                                            at)
               !kind_idx: is a vector ordered according to atom-type groups (as in natoms_per_kind),
               !given the index of the atom in the group it returns the index of the atom
               !in the cube file

               !accumulate model data
               a(icnt,at) = a(icnt,at) + pot
              end do
              !accumulate reference data
              !esp_sign depends on the program used to generate the data
              b(icnt) = cube%grid(xx,yy,zz) * (autokcal) * (esp_sign)
             endif
            end do
           end do
          end do
          call cpu_time(t2)

          write(ou,*), "Time to process points: ", t2-t1
          write(ou,*), ""
          
          ! print ESP map if wanted
          if (print_pot) then
            do i = 1, n_good_grid_pts
              do at = 1, size(natoms_per_kind)
               write(7588,'(f12.6, 1X)', advance="no") a(i,at)
              end do
              write(7588,*)
            end do
          end if
        end if

        ! once A and B are available, we can generate ATA, ATB, BTB to be used by the solver
        call  get_ata_matrix_post(n_good_grid_pts, nkinds, print_pot, atb, ata, btb, avg_a)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Actually computes ATA, ATB, BTB and relevant averages, filtering data from
    ! ESP values computed for first vdW scaling factor (gamma).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_ata_matrix_post(n_good_grid_pts, nkinds, print_pot, atb, ata, btb, avg_a)
      integer :: xx, yy, zz, icnt, i, at
      real(RDbl), allocatable, dimension(:,:) :: tmp_a
      real(RDbl), allocatable, dimension(:) :: tmp_b
        
      integer, intent(in) :: n_good_grid_pts
      !cf2py intent in n_good_grid_pts
      integer, intent(in) :: nkinds
      !cf2py intent in nkinds
      logical, intent(in) :: print_pot
      !cf2py intent in print_pot
      real(kind=RDbl), dimension(nkinds), intent(out) :: atb
      !cf2py intent out atb
      real(kind=RDbl), dimension(nkinds,nkinds), intent(out) :: ata
      !cf2py intent out ata
      real(kind=RDbl), intent(out) :: btb
      !cf2py intent out btb
      real(kind=RDbl), dimension(nkinds), intent(out) :: avg_a
      !cf2py intent out avg_a
      
      real(kind=RDbl) :: avg_b

      !allocate temporary matrices for gammas after the first
      if (first_ata_called) then 
        allocate(tmp_a(n_good_grid_pts,nkinds), tmp_b(n_good_grid_pts), stat=ierr)
         if (ierr /= 0) then
             write(ou,*), "Could not allocate tmp_a and tmp_b matrices."
             stop
         end if
      end if

      !loop over all grid points
      icnt = 0
      do xx=1,cube%npoints(1)
       do yy=1,cube%npoints(2)
        do zz=1,cube%npoints(3)

         !for good grid points get ESP data from first gamma
         if (rm_mask(xx,yy,zz) > 0) then
           icnt = icnt + 1
           !write(1111,*) icnt, rm_mask(xx,yy,zz)
           tmp_b(icnt) = b(rm_mask(xx,yy,zz))
           do at = 1, size(tmp_a,2)
            tmp_a(icnt,at) = a(rm_mask(xx,yy,zz),at)
           end do              
         end if
        end do
       end do
      end do

      if (icnt /= n_good_grid_pts) then
        write(ou,*) "ERROR: computed points do not agree with n_good_grid_pts"
        stop
      end if
      if (print_pot) then
        open(unit=8588,file="pot_from_post.dat")
        do i = 1, size(tmp_a,1)
          do at = 1, size(tmp_a,2)
           write(8588,'(f12.6, 1X)', advance="no") tmp_a(i,at)
          end do
          write(8588,*)
        end do
        close(8588)
      end if

      !compute ATA, ATB, BTB and averages
      avg_a = sum(tmp_a, 1)/n_good_grid_pts
      do at = 1, size(tmp_a,2)
        tmp_a(:,at) = tmp_a(:,at) - avg_a(at)
      end do
      avg_b = sum(tmp_b)/n_good_grid_pts
      tmp_b = tmp_b - avg_b
      !write(ou,*), "avg_b", avg_b
      ata = matmul(transpose(tmp_a),tmp_a)
      atb = matmul(transpose(tmp_a),tmp_b)
      btb = dot_product(tmp_b,tmp_b)

      deallocate(tmp_a,tmp_b)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Take dipole signal and fold it into the branch of the first dipole.
    ! Input data in Debye, output in eA.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tdf_smooth(tdf_name_in, is_cp2k, dpx, dpy, dpz, ndip, print_tdf)
      character(len=80), intent(in) :: tdf_name_in
      !f2py intent in tdf_name_in
      logical :: is_cp2k
      !f2py intent in is_cp2k
      real(kind=RDbl), dimension(ndip), intent(out) :: dpx
      !cf2py intent out dpx
      real(kind=RDbl), dimension(ndip), intent(out) :: dpy
      !cf2py intent out dpy
      real(kind=RDbl), dimension(ndip), intent(out) :: dpz
      !cf2py intent out dpz
      integer, intent(in) :: ndip
      !cf2py intent in ndip
      logical, optional :: print_tdf
      !cf2py intent in print_tdf
      integer           :: i, ierr, nstep
      real(kind=RDbl)   :: ftrash, px, py, pz, ipx, ipy, ipz,     &
                           xpre, ypre, zpre, xaft, yaft, zaft,    &
                           rx, ry, rz, x0, y0, z0
      character(len=80) :: dum
      logical           :: end_reached, op


      !open the output file for the entire module
      inquire(unit=ou, opened=op)
      if (.not. op) open(unit=ou, file="atautils.log")
      write(ou,*) "Processing TDF..."

      !get periodic cell data
      open(75,file=tdf_name_in,action='read',status='old')
      end_reached = .FALSE.
      if (is_cp2k) then !reading from CP2K output
          write(ou,*) "Assuming total dipole fluctuations coming from CP2K"
        do i=1,5
         read (75,*)
        end do
        read (75,*) dum,dum,px
        read (75,*) dum,ftrash,py !the '=' makes only one 'dum'
        read (75,*) dum,dum,ftrash,ftrash,pz
        write(ou,'(3(A5,f20.10))') '  px=', px,', py=', py,', pz=', pz
      else !reading from .tdf format
        write(ou,*) "Assuming InfiniCharges TDF format (in Debye)"
        read (75,*) px, ftrash, ftrash
        read (75,*) ftrash, py, ftrash
        read (75,*) ftrash, ftrash, pz
      end if

      !convert Debye to eA
      px=px*dtoea
      py=py*dtoea
      pz=pz*dtoea

      ipx=1.0d0/px !inverse of period along x
      ipy=1.0d0/py
      ipz=1.0d0/pz

      read(75,*) !header for both file types
      !read first total dipole for refolding, then rewind
      if (is_cp2k) then
          read(75,*) dum,xpre,dum,ypre,dum,zpre
          rewind(75)
      else !InfiniCharges .tdf
          read(75,*) xpre, ypre, zpre
          rewind(75)
          read(75,*) !skip header
          read(75,*) !skip header
          read(75,*) !skip header
          read(75,*) !skip header
      end if

      !store first dipole value
      xpre=xpre*dtoea
      ypre=ypre*dtoea
      zpre=zpre*dtoea
      x0=0.0d0
      y0=0.0d0
      z0=0.0d0
      write(ou,'(3f20.10)')  xpre, ypre, zpre
      
      !loop over all dipoles
      nstep=0
      do
        if (is_cp2k) then !for CP2K there is one dipole every 10 lines 
          do i=1,9
            read(75,*,iostat=ierr)
            if (ierr < 0) then 
              !print *, "End of file reached."
              end_reached = .TRUE.
              exit
             endif
          end do
          if (.not. end_reached) then !the dipole is the last line of the 10
            read(75,*,iostat=ierr) dum,xaft,dum,yaft,dum,zaft 
            !print *, ierr, dum,xaft,dum,yaft,dum,zaft
            if (ierr /= 0) then
              print *, "Something went wrong."
              print *, ierr, xaft,yaft,zaft
              exit
            endif
          else
            exit
          end if
        else !InfiniCharges .tdf file (read line by line)
          read(75,*,iostat=ierr) xaft, yaft, zaft
          if (ierr < 0) exit
        end if
        !actual dipole
        xaft=xaft*dtoea
        yaft=yaft*dtoea
        zaft=zaft*dtoea
        nstep=nstep+1

        !delta in dipoles
        rx=xaft-xpre
        ry=yaft-ypre
        rz=zaft-zpre
        !write(*,'(i10,3f20.10)') nstep, rx,ry,rz
        !folding
        rx=rx-px*dnint(rx*ipx)
        ry=ry-py*dnint(ry*ipy)
        rz=rz-pz*dnint(rz*ipz)
        dpx(nstep) = rx+xpre-x0
        dpy(nstep) = ry+ypre-y0
        dpz(nstep) = rz+zpre-z0
        !x0 can be used in writing final results to have origin in 0.0

        !reset pre
        xpre=rx+xpre !preceding value, plus shift due to minimum image
        ypre=ry+ypre
        zpre=rz+zpre
      end do
      close(75)
        
      if (nstep /= ndip) then
          write(ou,*) "ERROR: wrong number of dipoles: ", nstep, " while expected: ", ndip
          stop
      end if

      !print data if wanted
      if ( present(print_tdf) .and. print_tdf )then
        open(78,file='check_tdf.dat')
        do i = 1, nstep
          write(78,'(I10,3F20.10)') i, dpx(i), dpy(i), dpz(i)
        end do
        close(78)
      end if 

      write(ou,*) 'Dipoles counted = ', nstep
      write(ou,*) 'TDF processing completed.'
      write(ou,*) ''
   
    end subroutine tdf_smooth


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Return a vector containing the atomic number of the atoms in the cube file.
    ! Required for the python interface, which cannot work with Fortran pointers.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_zatom(natoms, zatom)
      integer, intent(in) :: natoms
      !cf2py intent in natoms
      integer, dimension(natoms), intent(out) :: zatom
      !cf2py intent out Zatom
      if ( .not. associated(cube) ) then
        write(ou,*), "ERROR: In Zatom, cube file not allocated"
        stop
      end if
      if (natoms /= cube%Natom) then
        write(ou,*), "ERROR: In Zatom, wrong number of atoms"
        stop
      end if
        
      zatom = cube%Zatom

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read an .xyz trajectory from fname.
    ! Requires number of atoms and number of frames.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_xyz(fname, natoms, nstep, types, rx, ry, rz)
      character(len=80), intent(in) :: fname
      !f2py intent in fname
      integer, intent(in) :: natoms
      !cf2py intent in natoms
      integer, intent(in) :: nstep
      !cf2py intent in nstep
      !character(len=5), dimension(natoms), intent(out) :: types
      !!f2py character(len=5), intent(out), dimension(natoms) :: types
      character, dimension(natoms,80), intent(out) :: types
      !f2py character, intent(out), dimension(natoms,80) :: types
      real(kind=rdbl), dimension(natoms,nstep), intent(out) :: rx, ry, rz
      !cf2py intent out rx, ry, rz
      integer :: tmp_nat, i, j
      character(len=80), dimension(natoms) :: dum
      
      open(75,file=fname,action='read',status="old")

      !loop over number of frames
      do i = 1, nstep
          read(75,*)tmp_nat
          if (tmp_nat /= natoms) then
              write(ou,*)"ERROR: number of atoms given ", natoms,&
                  " differs from the number of atoms read ", tmp_nat, &
                  " in ", fname
              stop
          end if
          read(75,*)
          !loop over atoms in given frame
          do j = 1, natoms
              read(75,*)dum(j), rx(j,i), ry(j,i), rz(j,i)
          end do
      end do
      close(75)

      !trick to pass the atom types through the python interface
      do j = 1, natoms
          do i = 1,80
            types(j,i) = dum(j)(i:i)
          end do
      end do

    end subroutine read_xyz


end module
