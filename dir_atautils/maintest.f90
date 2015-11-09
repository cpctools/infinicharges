program maintest
use atautils
      character(len=80) :: filename
      integer :: n_good_grid_pts
      integer :: kmax = 15
      real(kind=8) :: kappa = 6.7d0
      real(kind=8) :: locut=1.0d-10
      real(kind=8) :: rcut = -1.0d0
      integer, parameter :: nkinds = 7
      integer, parameter :: natoms = 276
      integer, dimension(nkinds) :: natoms_per_kind
      integer, dimension(natoms) :: kind_idx
      logical :: print_pot = .True.
      real(kind=8), dimension(nkinds) :: atb
      real(kind=8), dimension(nkinds,nkinds) :: ata
      real(kind=8)  :: btb
      real(kind=8), dimension(nkinds) :: avg_a

      natoms_per_kind = (/24,48,24,48,72,48,12/)
      kind_idx = (/ (i, i=1,natoms)/)

      filename = "../energy-force-v_hartree-1_0.cube"
      call init(filename)
      call filter(1.0d0, .False., n_good_grid_pts)
      call get_ata_matrix(kmax, kappa, locut, rcut, n_good_grid_pts, natoms_per_kind, nkinds, &
          kind_idx, natoms, print_pot, atb, ata, btb, avg_a)
      call release()
end program
