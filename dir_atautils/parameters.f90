module parameters
 implicit none
 integer,         parameter :: RDbl  = selected_real_kind(10, 50)
 real(kind=RDbl), parameter :: pi       =  3.14159265358979323846264338_RDbl ! Pi
 !real(kind=RDbl), parameter :: music_pi    = 3.1415927_RDbl  ! this is the value used in MUSIC

 real(kind=RDbl), parameter :: twopi = 2.0_RDbl*pi
 real(kind=RDbl), parameter :: zero  = 0.0_RDbl
 real(kind=Rdbl), parameter :: one   = 1.0_RDbl
 real (RDbl), parameter :: bohrtoangs = 0.52917720859d0
 real (RDbl), parameter :: dtoea = 0.20822678 !debye to e*Angstrom
 real (RDbl), parameter :: autokcal=627.509468713739d0
                             
 real (RDbl), parameter :: ech = 1.60217648700000e-19_RDbl ! Elementary charge [C]
 real (RDbl), parameter :: eps0 = 8.85418781762039e-12_RDbl ! Electric constant or permittivity of vacuum [F/m] or [C^2/(J*m)]
 real (RDbl), parameter :: avo = 6.02214179000000e+23_RDbl ! Avogadro constant [1/mol]
 real (RDbl), parameter :: mtoa = 1.0e10_RDbl ! Meters to Angstrom
 real (RDbl), parameter :: kcaltoj = 4.184e3_RDbl
 !e2kcal converts from e^2/Ang to kcal/mol. Remember that the charges'
 !values are in elementary charge units "ech".  
 !Usual derivations of ewald (included the one in atautils) drop the 4*pi*eps0 term
 !from the equations, where eps0 is in C^2/(J*m).
 !This yields units of e^2/Ang. The conversion factor, therefore, is
 !332.06370407376090 kcal/mol:
 real (RDbl), parameter :: e2kcal = ech**2*mtoa*avo/(kcaltoj*eps0*4.0_RDbl*pi)
 
 !real(kind=RDbl), parameter :: music_e2kcal = 331.5 ! this is the value used in MUSIC


 !van der Waals radii (Bohr units) taken from UFF database
 !through r = x_I/bohrtoangs/2
 !J. Am. Chem. Soc., 1992, 114 (25), pp 10024–10035, DOI:10.1021/ja00051a040
 double precision atm_radi(103)
 parameter(atm_radi=&
(/     2.72687, &
       2.23177, &
       2.31586, &
       2.59365, &
       3.85788, &
       3.63867, &
       3.4582, &
       3.30702, &
       3.17852, &
       3.06419, &
       2.81853, &
       2.85443, &
       4.25094, &
       4.05819, &
       3.91835, &
       3.81252, &
       3.72937, &
       3.65473, &
       3.60182, &
       3.21159, &
       3.11332, &
       2.99994, &
       2.97065, &
       2.85632, &
       2.79774, &
       2.75144, &
       2.71365, &
       2.67774, &
       3.3023, &
       2.61066, &
       4.14133, &
       4.04401, &
       3.99677, &
       3.97315, &
       3.95803, &
       3.91268, &
       3.88717, &
       3.44025, &
       3.16057, &
       2.95175, &
       2.99049, &
       2.88372, &
       2.8327, &
       2.79963, &
       2.7675, &
       2.73916, &
       2.97443, &
       2.69097, &
       4.21692, &
       4.14984, &
       4.17629, &
       4.22354, &
       4.25188, &
       4.16118, &
       4.26795, &
       3.49883, &
       3.32781, &
       3.35993, &
       3.40718, &
       3.37789, &
       3.35143, &
       3.32592, &
       3.30041, &
       3.1823, &
       3.26072, &
       3.23899, &
       3.22104, &
       3.20403, &
       3.18797, &
       3.17002, &
       3.4393, &
       2.96781, &
       2.99522, &
       2.89978, &
       2.79113, &
       2.94797, &
       2.68341, &
       2.60215, &
       3.11143, &
       2.55585, &
       4.10732, &
       4.06008, &
       4.12905, &
       4.44936, &
       4.4881, &
       4.50227, &
       4.62983, &
       3.47426, &
       3.28623, &
       3.20875, &
       3.23521, &
       3.20781, &
       3.23521, &
       3.23521, &
       3.19458, &
       3.14261, &
       3.1549, &
       3.13033, &
       3.1171, &
       3.10482, &
       3.09348, &
       3.06892, &
       3.05758/))

 contains
 
end module parameters
