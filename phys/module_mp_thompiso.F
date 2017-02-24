!+---+-----------------------------------------------------------------+
!.. This subroutine computes the moisture tendencies of water vapor,
!.. cloud droplets, rain, cloud ice (pristine), snow, and graupel.
!.. Prior to WRFv2.2 this code was based on Reisner et al (1998), but
!.. few of those pieces remain.  A complete description is now found in
!.. Thompson, G., P. R. Field, R. M. Rasmussen, and W. D. Hall, 2008:
!.. Explicit Forecasts of winter precipitation using an improved bulk
!.. microphysics scheme. Part II: Implementation of a new snow
!.. parameterization.  Mon. Wea. Rev., 136, 5095-5115.
!.. Prior to WRFv3.1, this code was single-moment rain prediction as
!.. described in the reference above, but in v3.1 and higher, the
!.. scheme is two-moment rain (predicted rain number concentration).
!..
!.. Most importantly, users may wish to modify the prescribed number of
!.. cloud droplets (Nt_c; see guidelines mentioned below).  Otherwise,
!.. users may alter the rain and graupel size distribution parameters
!.. to use exponential (Marshal-Palmer) or generalized gamma shape.
!.. The snow field assumes a combination of two gamma functions (from
!.. Field et al. 2005) and would require significant modifications
!.. throughout the entire code to alter its shape as well as accretion
!.. rates.  Users may also alter the constants used for density of rain,
!.. graupel, ice, and snow, but the latter is not constant when using
!.. Paul Field's snow distribution and moments methods.  Other values
!.. users can modify include the constants for mass and/or velocity
!.. power law relations and assumed capacitances used in deposition/
!.. sublimation/evaporation/melting.
!.. Remaining values should probably be left alone.
!..
!..Author: Greg Thompson, NCAR-RAL, gthompsn@ucar.edu, 303-497-2805
!..Last modified: 06 Sep 2013
!+---+-----------------------------------------------------------------+
!wrft:model_layer:physics
!+---+-----------------------------------------------------------------+
!
      MODULE module_mp_thompiso
      USE module_wrf_error
      USE module_mp_radar
     USE module_utility, ONLY: WRFU_Clock, WRFU_Alarm
      USE module_domain, ONLY : HISTORY_ALARM, Is_alarm_tstep
      use module_mp_water_isotopologues
      IMPLICIT NONE
      LOGICAL, PARAMETER, PRIVATE:: iiwarm = .false.
      INTEGER, PARAMETER, PRIVATE:: IFDRY = 0
      REAL, PARAMETER, PRIVATE:: T_0 = 273.15
      REAL, PARAMETER, PRIVATE:: PI = 3.1415926536
!..Densities of rain, snow, graupel, and cloud ice.
      REAL, PARAMETER, PRIVATE:: rho_w = 1000.0
      REAL, PARAMETER, PRIVATE:: rho_s = 100.0
      REAL, PARAMETER, PRIVATE:: rho_g = 500.0
      REAL, PARAMETER, PRIVATE:: rho_i = 890.0
!..Prescribed number of cloud droplets.  Set according to known data or
!.. roughly 100 per cc (100.E6 m^-3) for Maritime cases and
!.. 300 per cc (300.E6 m^-3) for Continental.  Gamma shape parameter,
!.. mu_c, calculated based on Nt_c is important in autoconversion
!.. scheme.
!bloss(2014-11): Add namelist for setting Nt_c and possibly other variables.
      REAL, PRIVATE:: Nt_c = 100.E6
!..Generalized gamma distributions for rain, graupel and cloud ice.
!.. N(D) = N_0 * D**mu * exp(-lamda*D);  mu=0 is exponential.
      REAL, PARAMETER, PRIVATE:: mu_r = 0.0
      REAL, PARAMETER, PRIVATE:: mu_g = 0.0
      REAL, PARAMETER, PRIVATE:: mu_i = 0.0
      REAL, PRIVATE:: mu_c
!..Sum of two gamma distrib for snow (Field et al. 2005).
!.. N(D) = M2**4/M3**3 * [Kap0*exp(-M2*Lam0*D/M3)
!..    + Kap1*(M2/M3)**mu_s * D**mu_s * exp(-M2*Lam1*D/M3)]
!.. M2 and M3 are the (bm_s)th and (bm_s+1)th moments respectively
!.. calculated as function of ice water content and temperature.
      REAL, PARAMETER, PRIVATE:: mu_s = 0.6357
      REAL, PARAMETER, PRIVATE:: Kap0 = 490.6
      REAL, PARAMETER, PRIVATE:: Kap1 = 17.46
      REAL, PARAMETER, PRIVATE:: Lam0 = 20.78
      REAL, PARAMETER, PRIVATE:: Lam1 = 3.29
!..Y-intercept parameter for graupel is not constant and depends on
!.. mixing ratio.  Also, when mu_g is non-zero, these become equiv
!.. y-intercept for an exponential distrib and proper values are
!.. computed based on same mixing ratio and total number concentration.
      REAL, PARAMETER, PRIVATE:: gonv_min = 1.E4
      REAL, PARAMETER, PRIVATE:: gonv_max = 3.E6
!..Mass power law relations:  mass = am*D**bm
!.. Snow from Field et al. (2005), others assume spherical form.
      REAL, PARAMETER, PRIVATE:: am_r = PI*rho_w/6.0
      REAL, PARAMETER, PRIVATE:: bm_r = 3.0
      REAL, PARAMETER, PRIVATE:: am_s = 0.069
      REAL, PARAMETER, PRIVATE:: bm_s = 2.0
      REAL, PARAMETER, PRIVATE:: am_g = PI*rho_g/6.0
      REAL, PARAMETER, PRIVATE:: bm_g = 3.0
      REAL, PARAMETER, PRIVATE:: am_i = PI*rho_i/6.0
      REAL, PARAMETER, PRIVATE:: bm_i = 3.0
!..Fallspeed power laws relations:  v = (av*D**bv)*exp(-fv*D)
!.. Rain from Ferrier (1994), ice, snow, and graupel from
!.. Thompson et al (2008). Coefficient fv is zero for graupel/ice.
      REAL, PARAMETER, PRIVATE:: av_r = 4854.0
      REAL, PARAMETER, PRIVATE:: bv_r = 1.0
      REAL, PARAMETER, PRIVATE:: fv_r = 195.0
      REAL, PARAMETER, PRIVATE:: av_s = 40.0
      REAL, PARAMETER, PRIVATE:: bv_s = 0.55
      REAL, PARAMETER, PRIVATE:: fv_s = 100.0
      REAL, PARAMETER, PRIVATE:: av_g = 442.0
      REAL, PARAMETER, PRIVATE:: bv_g = 0.89
      REAL, PARAMETER, PRIVATE:: av_i = 1847.5
      REAL, PARAMETER, PRIVATE:: bv_i = 1.0
!..Capacitance of sphere and plates/aggregates: D**3, D**2
      REAL, PARAMETER, PRIVATE:: C_cube = 0.5
      REAL, PARAMETER, PRIVATE:: C_sqrd = 0.3
!..Collection efficiencies.  Rain/snow/graupel collection of cloud
!.. droplets use variables (Ef_rw, Ef_sw, Ef_gw respectively) and
!.. get computed elsewhere because they are dependent on stokes
!.. number.
      REAL, PARAMETER, PRIVATE:: Ef_si = 0.05
      REAL, PARAMETER, PRIVATE:: Ef_rs = 0.95
      REAL, PARAMETER, PRIVATE:: Ef_rg = 0.75
      REAL, PARAMETER, PRIVATE:: Ef_ri = 0.95
!..Minimum microphys values
!.. R1 value, 1.E-12, cannot be set lower because of numerical
!.. problems with Paul Field's moments and should not be set larger
!.. because of truncation problems in snow/ice growth.
      REAL, PARAMETER, PRIVATE:: R1 = 1.E-18 !1.E-12
      REAL, PARAMETER, PRIVATE:: R2 = 1.E-6
      REAL, PARAMETER, PRIVATE:: eps = 1.E-15
!..Constants in Cooper curve relation for cloud ice number.
      REAL, PARAMETER, PRIVATE:: TNO = 5.0
      REAL, PARAMETER, PRIVATE:: ATO = 0.304
!..Rho_not used in fallspeed relations (rho_not/rho)**.5 adjustment.
      REAL, PARAMETER, PRIVATE:: rho_not = 101325.0/(287.05*298.0)
!..Schmidt number
      REAL, PARAMETER, PRIVATE:: Sc = 0.632
      REAL, PRIVATE:: Sc3
!..Homogeneous freezing temperature
      REAL, PARAMETER, PRIVATE:: HGFR = 235.16
!..Water vapor and air gas constants at constant pressure
      REAL, PARAMETER, PRIVATE:: Rv = 461.5
      REAL, PARAMETER, PRIVATE:: oRv = 1./Rv
      REAL, PARAMETER, PRIVATE:: R = 287.04
      REAL, PARAMETER, PRIVATE:: Cp = 1004.0
!..Enthalpy of sublimation, vaporization, and fusion at 0C.
      REAL, PARAMETER, PRIVATE:: lsub = 2.834E6
      REAL, PARAMETER, PRIVATE:: lvap0 = 2.5E6
      REAL, PARAMETER, PRIVATE:: lfus = lsub - lvap0
      REAL, PARAMETER, PRIVATE:: olfus = 1./lfus
!..Ice initiates with this mass (kg), corresponding diameter calc.
!..Min diameters and mass of cloud, rain, snow, and graupel (m, kg).
      REAL, PARAMETER, PRIVATE:: xm0i = 1.E-12
      REAL, PARAMETER, PRIVATE:: D0c = 1.E-6
      REAL, PARAMETER, PRIVATE:: D0r = 50.E-6
      REAL, PARAMETER, PRIVATE:: D0s = 200.E-6
      REAL, PARAMETER, PRIVATE:: D0g = 250.E-6
      REAL, PRIVATE:: D0i, xm0s, xm0g
!..Lookup table dimensions
      INTEGER, PARAMETER, PRIVATE:: nbins = 100
      INTEGER, PARAMETER, PRIVATE:: nbc = nbins
      INTEGER, PARAMETER, PRIVATE:: nbi = nbins
      INTEGER, PARAMETER, PRIVATE:: nbr = nbins
      INTEGER, PARAMETER, PRIVATE:: nbs = nbins
      INTEGER, PARAMETER, PRIVATE:: nbg = nbins
      INTEGER, PARAMETER, PRIVATE:: ntb_c = 37
      INTEGER, PARAMETER, PRIVATE:: ntb_i = 64
      INTEGER, PARAMETER, PRIVATE:: ntb_r = 37
      INTEGER, PARAMETER, PRIVATE:: ntb_s = 28
      INTEGER, PARAMETER, PRIVATE:: ntb_g = 28
      INTEGER, PARAMETER, PRIVATE:: ntb_g1 = 28
      INTEGER, PARAMETER, PRIVATE:: ntb_r1 = 37
      INTEGER, PARAMETER, PRIVATE:: ntb_i1 = 55
      INTEGER, PARAMETER, PRIVATE:: ntb_t = 9
      INTEGER, PRIVATE:: nic2, nii2, nii3, nir2, nir3, nis2, nig2, nig3
      DOUBLE PRECISION, DIMENSION(nbins+1):: xDx
      DOUBLE PRECISION, DIMENSION(nbc):: Dc, dtc
      DOUBLE PRECISION, DIMENSION(nbi):: Di, dti
      DOUBLE PRECISION, DIMENSION(nbr):: Dr, dtr
      DOUBLE PRECISION, DIMENSION(nbs):: Ds, dts
      DOUBLE PRECISION, DIMENSION(nbg):: Dg, dtg
!..Lookup tables for cloud water content (kg/m**3).
      REAL, DIMENSION(ntb_c), PARAMETER, PRIVATE:: &
      r_c = (/1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
              1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
              1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
              1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
              1.e-2/)
!..Lookup tables for cloud ice content (kg/m**3).
      REAL, DIMENSION(ntb_i), PARAMETER, PRIVATE:: &
      r_i = (/1.e-10,2.e-10,3.e-10,4.e-10, &
              5.e-10,6.e-10,7.e-10,8.e-10,9.e-10, &
              1.e-9,2.e-9,3.e-9,4.e-9,5.e-9,6.e-9,7.e-9,8.e-9,9.e-9, &
              1.e-8,2.e-8,3.e-8,4.e-8,5.e-8,6.e-8,7.e-8,8.e-8,9.e-8, &
              1.e-7,2.e-7,3.e-7,4.e-7,5.e-7,6.e-7,7.e-7,8.e-7,9.e-7, &
              1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
              1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
              1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
              1.e-3/)
!..Lookup tables for rain content (kg/m**3).
      REAL, DIMENSION(ntb_r), PARAMETER, PRIVATE:: &
      r_r = (/1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
              1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
              1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
              1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
              1.e-2/)
!..Lookup tables for graupel content (kg/m**3).
      REAL, DIMENSION(ntb_g), PARAMETER, PRIVATE:: &
      r_g = (/1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
              1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
              1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
              1.e-2/)
!..Lookup tables for snow content (kg/m**3).
      REAL, DIMENSION(ntb_s), PARAMETER, PRIVATE:: &
      r_s = (/1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
              1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
              1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
              1.e-2/)
!..Lookup tables for rain y-intercept parameter (/m**4).
      REAL, DIMENSION(ntb_r1), PARAMETER, PRIVATE:: &
      N0r_exp = (/1.e6,2.e6,3.e6,4.e6,5.e6,6.e6,7.e6,8.e6,9.e6, &
                  1.e7,2.e7,3.e7,4.e7,5.e7,6.e7,7.e7,8.e7,9.e7, &
                  1.e8,2.e8,3.e8,4.e8,5.e8,6.e8,7.e8,8.e8,9.e8, &
                  1.e9,2.e9,3.e9,4.e9,5.e9,6.e9,7.e9,8.e9,9.e9, &
                  1.e10/)
!..Lookup tables for graupel y-intercept parameter (/m**4).
      REAL, DIMENSION(ntb_g1), PARAMETER, PRIVATE:: &
      N0g_exp = (/1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4, &
                  1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5, &
                  1.e6,2.e6,3.e6,4.e6,5.e6,6.e6,7.e6,8.e6,9.e6, &
                  1.e7/)
!..Lookup tables for ice number concentration (/m**3).
      REAL, DIMENSION(ntb_i1), PARAMETER, PRIVATE:: &
      Nt_i = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0, &
               1.e1,2.e1,3.e1,4.e1,5.e1,6.e1,7.e1,8.e1,9.e1, &
               1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2, &
               1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3, &
               1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4, &
               1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5, &
               1.e6/)
!..For snow moments conversions (from Field et al. 2005)
      REAL, DIMENSION(10), PARAMETER, PRIVATE:: &
      sa = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
              0.31255,   0.000204,  0.003199, 0.0,      -0.015952/)
      REAL, DIMENSION(10), PARAMETER, PRIVATE:: &
      sb = (/ 0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
              0.060366,  0.000079,  0.000594, 0.0,      -0.003577/)
!..Temperatures (5 C interval 0 to -40) used in lookup tables.
      REAL, DIMENSION(ntb_t), PARAMETER, PRIVATE:: &
      Tc = (/-0.01, -5., -10., -15., -20., -25., -30., -35., -40./)
!..Lookup tables for various accretion/collection terms.
!.. ntb_x refers to the number of elements for rain, snow, graupel,
!.. and temperature array indices.  Variables beginning with t-p/c/m/n
!.. represent lookup tables.  Save compile-time memory by making
!.. allocatable (2009Jun12, J. Michalakes).
      INTEGER, PARAMETER, PRIVATE:: R8SIZE = 8
      REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:,:)::             &
                tcg_racg, tmr_racg, tcr_gacr, tmg_gacr,                 &
                tnr_racg, tnr_gacr
      REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:,:)::             &
                tcs_racs1, tmr_racs1, tcs_racs2, tmr_racs2,             &
                tcr_sacr1, tms_sacr1, tcr_sacr2, tms_sacr2,             &
                tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2
      REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:)::                 &
                tpi_qcfz, tni_qcfz
      REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:)::               &
                tpi_qrfz, tpg_qrfz, tni_qrfz, tnr_qrfz
      REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:)::                 &
                tps_iaus, tni_iaus, tpi_ide
      REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:):: t_Efrw
      REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:):: t_Efsw
      REAL (KIND=R8SIZE), ALLOCATABLE, DIMENSION(:,:,:):: tnr_rev
!..Variables holding a bunch of exponents and gamma values (cloud water,
!.. cloud ice, rain, snow, then graupel).
      REAL, DIMENSION(3), PRIVATE:: cce, ccg
      REAL, PRIVATE::  ocg1, ocg2
      REAL, DIMENSION(7), PRIVATE:: cie, cig
      REAL, PRIVATE:: oig1, oig2, obmi
      REAL, DIMENSION(13), PRIVATE:: cre, crg
      REAL, PRIVATE:: ore1, org1, org2, org3, obmr
      REAL, DIMENSION(18), PRIVATE:: cse, csg
      REAL, PRIVATE:: oams, obms, ocms
      REAL, DIMENSION(12), PRIVATE:: cge, cgg
      REAL, PRIVATE:: oge1, ogg1, ogg2, ogg3, oamg, obmg, ocmg
!..Declaration of precomputed constants in various rate eqns.
      REAL:: t1_qr_qc, t1_qr_qi, t2_qr_qi, t1_qg_qc, t1_qs_qc, t1_qs_qi
      REAL:: t1_qr_ev, t2_qr_ev
      REAL:: t1_qs_sd, t2_qs_sd, t1_qg_sd, t2_qg_sd
      REAL:: t1_qs_me, t2_qs_me, t1_qg_me, t2_qg_me
      CHARACTER*256:: mp_debug
!+---+
!+---+-----------------------------------------------------------------+
!..END DECLARATIONS
!+---+-----------------------------------------------------------------+
!+---+
!ctrlL
      CONTAINS
      SUBROUTINE thompiso_init
      IMPLICIT NONE
      INTEGER:: i, j, k, m, n
      LOGICAL:: micro_init
      !bloss: allow storage of lookup tables in binary file
      !         Probably not a portable solution, but should work
      !         in development work on a single platform and allow
      !         this to be run as a parcel model.
      logical :: read_lookup_tables_from_file
      integer :: ioerror
      integer :: nbins_test, &
               nbr_test, nbs_test, nbg_test, &
               ntb_r_test, ntb_s_test, ntb_g_test, &
               ntb_g1_test, ntb_r1_test, ntb_t_test
      real :: mu_r_test, am_r_test, bm_r_test, av_r_test, bv_r_test, D0r_test, &
             mu_s_test, Kap0_test, Kap1_test, Lam0_test, Lam1_test, & ! extra snow parameters
             am_s_test, bm_s_test, av_s_test, bv_s_test, fv_s_test, D0s_test, &
             mu_g_test, am_g_test, bm_g_test, av_g_test, bv_g_test, D0g_test
!!$      integer :: nbins_test, nbc_test, &
!!$             nbi_test, nbr_test, nbs_test,  nbg_test, ntb_c_test, &
!!$             ntb_i_test, ntb_r_test, ntb_s_test, ntb_g_test, &
!!$             ntb_g1_test, ntb_r1_test, ntb_i1_test, ntb_t_test
!!$      real :: mu_c_test, mu_i_test, mu_r_test, mu_g_test 

      !bloss: add variables to allow input of Nt_c and possibly other
      !  variables through a namelist.
      logical:: have_namelist_file
      integer :: ierr, ierr_missing_namelist, placeholder

      !bloss: Supply parameters for topography and vertical grid through a
      !         namelist called "mp_thompson"
      NAMELIST /MP_THOMPSON/ Nt_c

      !bloss: Make a namelist that shouldn't ever be in the prm file.
      !  This lets you diagnose whether the namelist was missing or
      !  if there was an error reading the namelist.  My feeling is
      !  that there shouldn't be an error if the (optional) namelist is
      !  missing but that the run should stop if the MP_THOMPSON namelist
      !  was in the file but there was an error reading it.
      NAMELIST /bydswioubnuOZBSda/ placeholder

!..Allocate space for lookup tables (J. Michalakes 2009Jun08).
      micro_init = .FALSE.
      if (.NOT. ALLOCATED(tcg_racg) ) then
         ALLOCATE(tcg_racg(ntb_g1,ntb_g,ntb_r1,ntb_r))
         micro_init = .TRUE.
      endif
      if (.NOT. ALLOCATED(tmr_racg)) ALLOCATE(tmr_racg(ntb_g1,ntb_g,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tcr_gacr)) ALLOCATE(tcr_gacr(ntb_g1,ntb_g,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tmg_gacr)) ALLOCATE(tmg_gacr(ntb_g1,ntb_g,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tnr_racg)) ALLOCATE(tnr_racg(ntb_g1,ntb_g,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tnr_gacr)) ALLOCATE(tnr_gacr(ntb_g1,ntb_g,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tcs_racs1)) ALLOCATE(tcs_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tmr_racs1)) ALLOCATE(tmr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tcs_racs2)) ALLOCATE(tcs_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tmr_racs2)) ALLOCATE(tmr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tcr_sacr1)) ALLOCATE(tcr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tms_sacr1)) ALLOCATE(tms_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tcr_sacr2)) ALLOCATE(tcr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tms_sacr2)) ALLOCATE(tms_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tnr_racs1)) ALLOCATE(tnr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tnr_racs2)) ALLOCATE(tnr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tnr_sacr1)) ALLOCATE(tnr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tnr_sacr2)) ALLOCATE(tnr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
      if (.NOT. ALLOCATED(tpi_qcfz)) ALLOCATE(tpi_qcfz(ntb_c,45))
      if (.NOT. ALLOCATED(tni_qcfz)) ALLOCATE(tni_qcfz(ntb_c,45))
      if (.NOT. ALLOCATED(tpi_qrfz)) ALLOCATE(tpi_qrfz(ntb_r,ntb_r1,45))
      if (.NOT. ALLOCATED(tpg_qrfz)) ALLOCATE(tpg_qrfz(ntb_r,ntb_r1,45))
      if (.NOT. ALLOCATED(tni_qrfz)) ALLOCATE(tni_qrfz(ntb_r,ntb_r1,45))
      if (.NOT. ALLOCATED(tnr_qrfz)) ALLOCATE(tnr_qrfz(ntb_r,ntb_r1,45))
      if (.NOT. ALLOCATED(tps_iaus)) ALLOCATE(tps_iaus(ntb_i,ntb_i1))
      if (.NOT. ALLOCATED(tni_iaus)) ALLOCATE(tni_iaus(ntb_i,ntb_i1))
      if (.NOT. ALLOCATED(tpi_ide)) ALLOCATE(tpi_ide(ntb_i,ntb_i1))
      if (.NOT. ALLOCATED(t_Efrw)) ALLOCATE(t_Efrw(nbr,nbc))
      if (.NOT. ALLOCATED(t_Efsw)) ALLOCATE(t_Efsw(nbs,nbc))
      if (.NOT. ALLOCATED(tnr_rev)) ALLOCATE(tnr_rev(nbr, ntb_r1, ntb_r))
      if (micro_init) then

!bloss(2014-11): Add namelist for setting Nt_c and possibly other variables.
        ! check for file with namelist.  open and read in namelist if it exists.
        INQUIRE(FILE='prm',EXIST=have_namelist_file)
        if(have_namelist_file) then
          OPEN(UNIT=94,FILE='prm',FORM='formatted')
          READ(94,NML=MP_THOMPSON,IOSTAT=ierr)
          if(ierr.eq.0) then
            write(6,*) 'Read MP_THOMPSON namelist from prm file'
          else
            REWIND(94)
            READ(94,NML=bydswioubnuOZBSda,IOSTAT=ierr_missing_namelist)
            if(ierr_missing_namelist.eq.0) then
              write(6,*) 'Please remove the bydswioubnuOZBSda namelist from '
              write(6,*) ' your prm file.  Its absence is important in diagnosing'
              write(6,*) ' errors in reading the MP_THOMPSON namelist.'
              STOP 'REMOVE bydswioubnuOZBSda NAMELIST FROM PRM FILE'
            end if
            if(ierr.eq.ierr_missing_namelist) then
              write(6,*) 'No MP_THOMPSON namelist found in prm file'
              write(6,*) 'Running with default parameter choices'
            else
              write(6,*) 'Error reading MP_THOMPSON namelist from prm file'
              READ(94,NML=MP_THOMPSON)
            end if
          end if
          CLOSE(94)
        end if

        ! write namelist to standard output
        write(6,*) 'The following lists the namelist settings within MP_THOMPSON'
        write(6,NML=MP_THOMPSON)

!..From Martin et al. (1994), assign gamma shape parameter mu for cloud
!.. drops according to general dispersion characteristics (disp=~0.25
!.. for Maritime and 0.45 for Continental).
!.. disp=SQRT((mu+2)/(mu+1) - 1) so mu varies from 15 for Maritime
!.. to 2 for really dirty air.
      mu_c = MIN(15., (1000.E6/Nt_c + 2.))
!..Schmidt number to one-third used numerous times.
      Sc3 = Sc**(1./3.)
!..Compute min ice diam from mass, min snow/graupel mass from diam.
      D0i = (xm0i/am_i)**(1./bm_i)
      xm0s = am_s * D0s**bm_s
      xm0g = am_g * D0g**bm_g
!..These constants various exponents and gamma() assoc with cloud,
!.. rain, snow, and graupel.
      cce(1) = mu_c + 1.
      cce(2) = bm_r + mu_c + 1.
      cce(3) = bm_r + mu_c + 4.
      ccg(1) = WGAMMA(cce(1))
      ccg(2) = WGAMMA(cce(2))
      ccg(3) = WGAMMA(cce(3))
      ocg1 = 1./ccg(1)
      ocg2 = 1./ccg(2)
      cie(1) = mu_i + 1.
      cie(2) = bm_i + mu_i + 1.
      cie(3) = bm_i + mu_i + bv_i + 1.
      cie(4) = mu_i + bv_i + 1.
      cie(5) = mu_i + 2.
      cie(6) = bm_i*0.5 + mu_i + bv_i + 1.
      cie(7) = bm_i*0.5 + mu_i + 1.
      cig(1) = WGAMMA(cie(1))
      cig(2) = WGAMMA(cie(2))
      cig(3) = WGAMMA(cie(3))
      cig(4) = WGAMMA(cie(4))
      cig(5) = WGAMMA(cie(5))
      cig(6) = WGAMMA(cie(6))
      cig(7) = WGAMMA(cie(7))
      oig1 = 1./cig(1)
      oig2 = 1./cig(2)
      obmi = 1./bm_i
      cre(1) = bm_r + 1.
      cre(2) = mu_r + 1.
      cre(3) = bm_r + mu_r + 1.
      cre(4) = bm_r*2. + mu_r + 1.
      cre(5) = mu_r + bv_r + 1.
      cre(6) = bm_r + mu_r + bv_r + 1.
      cre(7) = bm_r*0.5 + mu_r + bv_r + 1.
      cre(8) = bm_r + mu_r + bv_r + 3.
      cre(9) = mu_r + bv_r + 3.
      cre(10) = mu_r + 2.
      cre(11) = 0.5*(bv_r + 5. + 2.*mu_r)
      cre(12) = bm_r*0.5 + mu_r + 1.
      cre(13) = bm_r*2. + mu_r + bv_r + 1.
      do n = 1, 13
         crg(n) = WGAMMA(cre(n))
      enddo
      obmr = 1./bm_r
      ore1 = 1./cre(1)
      org1 = 1./crg(1)
      org2 = 1./crg(2)
      org3 = 1./crg(3)
      cse(1) = bm_s + 1.
      cse(2) = bm_s + 2.
      cse(3) = bm_s*2.
      cse(4) = bm_s + bv_s + 1.
      cse(5) = bm_s*2. + bv_s + 1.
      cse(6) = bm_s*2. + 1.
      cse(7) = bm_s + mu_s + 1.
      cse(8) = bm_s + mu_s + 2.
      cse(9) = bm_s + mu_s + 3.
      cse(10) = bm_s + mu_s + bv_s + 1.
      cse(11) = bm_s*2. + mu_s + bv_s + 1.
      cse(12) = bm_s*2. + mu_s + 1.
      cse(13) = bv_s + 2.
      cse(14) = bm_s + bv_s
      cse(15) = mu_s + 1.
      cse(16) = 1.0 + (1.0 + bv_s)/2.
      cse(17) = cse(16) + mu_s + 1.
      cse(18) = bv_s + mu_s + 3.
      do n = 1, 18
         csg(n) = WGAMMA(cse(n))
      enddo
      oams = 1./am_s
      obms = 1./bm_s
      ocms = oams**obms
      cge(1) = bm_g + 1.
      cge(2) = mu_g + 1.
      cge(3) = bm_g + mu_g + 1.
      cge(4) = bm_g*2. + mu_g + 1.
      cge(5) = bm_g*2. + mu_g + bv_g + 1.
      cge(6) = bm_g + mu_g + bv_g + 1.
      cge(7) = bm_g + mu_g + bv_g + 2.
      cge(8) = bm_g + mu_g + bv_g + 3.
      cge(9) = mu_g + bv_g + 3.
      cge(10) = mu_g + 2.
      cge(11) = 0.5*(bv_g + 5. + 2.*mu_g)
      cge(12) = 0.5*(bv_g + 5.) + mu_g
      do n = 1, 12
         cgg(n) = WGAMMA(cge(n))
      enddo
      oamg = 1./am_g
      obmg = 1./bm_g
      ocmg = oamg**obmg
      oge1 = 1./cge(1)
      ogg1 = 1./cgg(1)
      ogg2 = 1./cgg(2)
      ogg3 = 1./cgg(3)
!+---+-----------------------------------------------------------------+
!..Simplify various rate eqns the best we can now.
!+---+-----------------------------------------------------------------+
!..Rain collecting cloud water and cloud ice
      t1_qr_qc = PI*.25*av_r * crg(9)
      t1_qr_qi = PI*.25*av_r * crg(9)
      t2_qr_qi = PI*.25*am_r*av_r * crg(8)
!..Graupel collecting cloud water
      t1_qg_qc = PI*.25*av_g * cgg(9)
!..Snow collecting cloud water
      t1_qs_qc = PI*.25*av_s
!..Snow collecting cloud ice
      t1_qs_qi = PI*.25*av_s
!..Evaporation of rain; ignore depositional growth of rain.
      t1_qr_ev = 0.78 * crg(10)
      t2_qr_ev = 0.308*Sc3*SQRT(av_r) * crg(11)
!..Sublimation/depositional growth of snow
      t1_qs_sd = 0.86
      t2_qs_sd = 0.28*Sc3*SQRT(av_s)
!..Melting of snow
      t1_qs_me = PI*4.*C_sqrd*olfus * 0.86
      t2_qs_me = PI*4.*C_sqrd*olfus * 0.28*Sc3*SQRT(av_s)
!..Sublimation/depositional growth of graupel
      t1_qg_sd = 0.86 * cgg(10)
      t2_qg_sd = 0.28*Sc3*SQRT(av_g) * cgg(11)
!..Melting of graupel
      t1_qg_me = PI*4.*C_cube*olfus * 0.86 * cgg(10)
      t2_qg_me = PI*4.*C_cube*olfus * 0.28*Sc3*SQRT(av_g) * cgg(11)
!..Constants for helping find lookup table indexes.
      nic2 = NINT(ALOG10(r_c(1)))
      nii2 = NINT(ALOG10(r_i(1)))
      nii3 = NINT(ALOG10(Nt_i(1)))
      nir2 = NINT(ALOG10(r_r(1)))
      nir3 = NINT(ALOG10(N0r_exp(1)))
      nis2 = NINT(ALOG10(r_s(1)))
      nig2 = NINT(ALOG10(r_g(1)))
      nig3 = NINT(ALOG10(N0g_exp(1)))
!..Create bins of cloud water (from min diameter up to 100 microns).
      Dc(1) = D0c*1.0d0
      dtc(1) = D0c*1.0d0
      do n = 2, nbc
         Dc(n) = Dc(n-1) + 1.0D-6
         dtc(n) = (Dc(n) - Dc(n-1))
      enddo
!..Create bins of cloud ice (from min diameter up to 5x min snow size).
      xDx(1) = D0i*1.0d0
      xDx(nbi+1) = 5.0d0*D0s
      do n = 2, nbi
         xDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbi) &
                  *DLOG(xDx(nbi+1)/xDx(1)) +DLOG(xDx(1)))
      enddo
      do n = 1, nbi
         Di(n) = DSQRT(xDx(n)*xDx(n+1))
         dti(n) = xDx(n+1) - xDx(n)
      enddo
!..Create bins of rain (from min diameter up to 5 mm).
      xDx(1) = D0r*1.0d0
      xDx(nbr+1) = 0.005d0
      do n = 2, nbr
         xDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbr) &
                  *DLOG(xDx(nbr+1)/xDx(1)) +DLOG(xDx(1)))
      enddo
      do n = 1, nbr
         Dr(n) = DSQRT(xDx(n)*xDx(n+1))
         dtr(n) = xDx(n+1) - xDx(n)
      enddo
!..Create bins of snow (from min diameter up to 2 cm).
      xDx(1) = D0s*1.0d0
      xDx(nbs+1) = 0.02d0
      do n = 2, nbs
         xDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbs) &
                  *DLOG(xDx(nbs+1)/xDx(1)) +DLOG(xDx(1)))
      enddo
      do n = 1, nbs
         Ds(n) = DSQRT(xDx(n)*xDx(n+1))
         dts(n) = xDx(n+1) - xDx(n)
      enddo
!..Create bins of graupel (from min diameter up to 5 cm).
      xDx(1) = D0g*1.0d0
      xDx(nbg+1) = 0.05d0
      do n = 2, nbg
         xDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbg) &
                  *DLOG(xDx(nbg+1)/xDx(1)) +DLOG(xDx(1)))
      enddo
      do n = 1, nbg
         Dg(n) = DSQRT(xDx(n)*xDx(n+1))
         dtg(n) = xDx(n+1) - xDx(n)
      enddo
!+---+-----------------------------------------------------------------+
!..Create lookup tables for most costly calculations.
!+---+-----------------------------------------------------------------+
      do m = 1, ntb_r
         do k = 1, ntb_r1
            do j = 1, ntb_g
               do i = 1, ntb_g1
                  tcg_racg(i,j,k,m) = 0.0d0
                  tmr_racg(i,j,k,m) = 0.0d0
                  tcr_gacr(i,j,k,m) = 0.0d0
                  tmg_gacr(i,j,k,m) = 0.0d0
                  tnr_racg(i,j,k,m) = 0.0d0
                  tnr_gacr(i,j,k,m) = 0.0d0
               enddo
            enddo
         enddo
      enddo
      do m = 1, ntb_r
         do k = 1, ntb_r1
            do j = 1, ntb_t
               do i = 1, ntb_s
                  tcs_racs1(i,j,k,m) = 0.0d0
                  tmr_racs1(i,j,k,m) = 0.0d0
                  tcs_racs2(i,j,k,m) = 0.0d0
                  tmr_racs2(i,j,k,m) = 0.0d0
                  tcr_sacr1(i,j,k,m) = 0.0d0
                  tms_sacr1(i,j,k,m) = 0.0d0
                  tcr_sacr2(i,j,k,m) = 0.0d0
                  tms_sacr2(i,j,k,m) = 0.0d0
                  tnr_racs1(i,j,k,m) = 0.0d0
                  tnr_racs2(i,j,k,m) = 0.0d0
                  tnr_sacr1(i,j,k,m) = 0.0d0
                  tnr_sacr2(i,j,k,m) = 0.0d0
               enddo
            enddo
         enddo
      enddo
      do k = 1, 45
         do j = 1, ntb_r1
            do i = 1, ntb_r
               tpi_qrfz(i,j,k) = 0.0d0
               tni_qrfz(i,j,k) = 0.0d0
               tpg_qrfz(i,j,k) = 0.0d0
               tnr_qrfz(i,j,k) = 0.0d0
            enddo
         enddo
         do i = 1, ntb_c
            tpi_qcfz(i,k) = 0.0d0
            tni_qcfz(i,k) = 0.0d0
         enddo
      enddo
      do j = 1, ntb_i1
         do i = 1, ntb_i
            tps_iaus(i,j) = 0.0d0
            tni_iaus(i,j) = 0.0d0
            tpi_ide(i,j) = 0.0d0
         enddo
      enddo
      do j = 1, nbc
         do i = 1, nbr
            t_Efrw(i,j) = 0.0
         enddo
         do i = 1, nbs
            t_Efsw(i,j) = 0.0
         enddo
      enddo
      do k = 1, ntb_r
         do j = 1, ntb_r1
            do i = 1, nbr
               tnr_rev(i,j,k) = 0.0d0
            enddo
         enddo
      enddo
      CALL wrf_debug(150, 'CREATING MICROPHYSICS LOOKUP TABLES ... ')
      WRITE (wrf_err_message, '(a, f5.2, a, f5.2, a, f5.2, a, f5.2)') &
          ' using: mu_c=',mu_c,' mu_i=',mu_i,' mu_r=',mu_r,' mu_g=',mu_g
      CALL wrf_debug(150, wrf_err_message)
!..Collision efficiency between rain/snow and cloud water.
      CALL wrf_debug(200, '  creating qc collision eff tables')
      call table_Efrw
      call table_Efsw
!..Drop evaporation.
!     CALL wrf_debug(200, '  creating rain evap table')
!     call table_dropEvap
!..Initialize various constants for computing radar reflectivity.
      xam_r = am_r
      xbm_r = bm_r
      xmu_r = mu_r
      xam_s = am_s
      xbm_s = bm_s
      xmu_s = mu_s
      xam_g = am_g
      xbm_g = bm_g
      xmu_g = mu_g
      call radar_init
      if (iiwarm) then
        CALL wrf_debug(200, '!!WARM MICROPHYSICS ONLY!!')
        CALL wrf_debug(200, '!!!!!NO ICE PROCESSES!!!!!')
      end if
      if (.not. iiwarm) then


        !bloss: Restrict lookup table file to qr_acr_qs and qr_acr_qg tables.
        !  This will keep the biggest lookup tables from being re-computed when
        !  changing the cloud droplet number concentration or cloud ice properties.
        INQUIRE(FILE='thompson_lookup_tables_r1.bin', & ! rename as _r1 to avoid confusion w/old files.
           EXIST=read_lookup_tables_from_file)

        if(read_lookup_tables_from_file) then
          CALL wrf_debug(200, '  TRYING TO GET LOOKUP TABLES FROM FILE')
          OPEN(UNIT=92,FILE='thompson_lookup_tables_r1.bin', &
               FORM='unformatted')

          ! read lookup table dimensions from file
          read(92) nbins_test, &
               nbr_test, nbs_test, nbg_test, &
               ntb_r_test, ntb_s_test, ntb_g_test, &
               ntb_g1_test, ntb_r1_test, ntb_t_test

        !check that the table dimensions in the file match those above.
        if(nbins.ne.nbins_test) read_lookup_tables_from_file = .false.
!!$        if(nbc.ne.nbc_test) read_lookup_tables_from_file = .false.
!!$        if(nbi.ne.nbi_test) read_lookup_tables_from_file = .false.
        if(nbr.ne.nbr_test) read_lookup_tables_from_file = .false.
        if(nbs.ne.nbs_test) read_lookup_tables_from_file = .false.
        if(nbg.ne.nbg_test) read_lookup_tables_from_file = .false.
!!$        if(ntb_c.ne.ntb_c_test) read_lookup_tables_from_file = .false.
!!$        if(ntb_i.ne.ntb_i_test) read_lookup_tables_from_file = .false.
        if(ntb_r.ne.ntb_r_test) read_lookup_tables_from_file = .false.
        if(ntb_s.ne.ntb_s_test) read_lookup_tables_from_file = .false.
        if(ntb_g.ne.ntb_g_test) read_lookup_tables_from_file = .false.
        if(ntb_g1.ne.ntb_g1_test) read_lookup_tables_from_file = .false.
        if(ntb_r1.ne.ntb_r1_test) read_lookup_tables_from_file = .false.
!!$        if(ntb_i1.ne.ntb_i1_test) read_lookup_tables_from_file = .false.
        if(ntb_t.ne.ntb_t_test) read_lookup_tables_from_file = .false.
        
        ! read microphysical parameters from file.
        read(92) mu_r_test, am_r_test, bm_r_test, av_r_test, bv_r_test, D0r_test, &
             mu_s_test, Kap0_test, Kap1_test, Lam0_test, Lam1_test, & ! extra snow parameters
             am_s_test, bm_s_test, av_s_test, bv_s_test, fv_s_test, D0s_test, &
             mu_g_test, am_g_test, bm_g_test, av_g_test, bv_g_test, D0g_test

        !check that the microphysical parameters used to compute the lookup table
        !   in the file match those above.  If not, recompute the table.
        if(mu_r.ne.mu_r_test) read_lookup_tables_from_file = .false. ! rain parameters
        if(am_r.ne.am_r_test) read_lookup_tables_from_file = .false.
        if(bm_r.ne.bm_r_test) read_lookup_tables_from_file = .false.
        if(av_r.ne.av_r_test) read_lookup_tables_from_file = .false.
        if(bv_r.ne.bv_r_test) read_lookup_tables_from_file = .false.
        if(D0r.ne.D0r_test) read_lookup_tables_from_file = .false.

        if(mu_s.ne.mu_s_test) read_lookup_tables_from_file = .false. ! snow parameters
        if(Kap0.ne.Kap0_test) read_lookup_tables_from_file = .false.
        if(Kap1.ne.Kap1_test) read_lookup_tables_from_file = .false.
        if(Lam0.ne.Lam0_test) read_lookup_tables_from_file = .false.
        if(Lam1.ne.Lam1_test) read_lookup_tables_from_file = .false.
        if(am_s.ne.am_s_test) read_lookup_tables_from_file = .false.
        if(bm_s.ne.bm_s_test) read_lookup_tables_from_file = .false.
        if(av_s.ne.av_s_test) read_lookup_tables_from_file = .false.
        if(bv_s.ne.bv_s_test) read_lookup_tables_from_file = .false.
        if(fv_s.ne.fv_s_test) read_lookup_tables_from_file = .false.
        if(D0s.ne.D0s_test) read_lookup_tables_from_file = .false.

        if(mu_g.ne.mu_g_test) read_lookup_tables_from_file = .false. ! graupel parameters
        if(am_g.ne.am_g_test) read_lookup_tables_from_file = .false.
        if(bm_g.ne.bm_g_test) read_lookup_tables_from_file = .false.
        if(av_g.ne.av_g_test) read_lookup_tables_from_file = .false.
        if(bv_g.ne.bv_g_test) read_lookup_tables_from_file = .false.
        if(D0g.ne.D0g_test) read_lookup_tables_from_file = .false.

        if(.NOT.read_lookup_tables_from_file) CLOSE(92)
      end if
      if(read_lookup_tables_from_file) then
        ! get tables computed in qr_acr_qg
        READ(92,IOSTAT=ioerror) tcg_racg, tmr_racg, tcr_gacr, tmg_gacr, &
             tnr_racg, tnr_gacr
        if(ioerror.ne.0) read_lookup_tables_from_file = .false.
        ! get tables computed in qr_acr_qs
        READ(92,IOSTAT=ioerror) tcs_racs1, tmr_racs1, tcs_racs2, tmr_racs2, &
             tcr_sacr1, tms_sacr1, tcr_sacr2, tms_sacr2, &
             tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2
        if(ioerror.ne.0) read_lookup_tables_from_file = .false.
        ! get tables computed in freezeH2O and qi_aut_qs
        READ(92,IOSTAT=ioerror) tpi_qcfz, tni_qcfz, tpi_qrfz, tpg_qrfz, &
             tni_qrfz, tnr_qrfz, tps_iaus, tni_iaus, tpi_ide
        if(ioerror.ne.0) read_lookup_tables_from_file = .false.
        CLOSE(92)
        if(read_lookup_tables_from_file) then
          CALL wrf_debug(200, '  READ RAIN-SNOW-GRAUPEL ACCRETION LOOKUP TABLES FROM FILE')
        end if
      end if
      if(.NOT.read_lookup_tables_from_file) then
!..Rain collecting graupel & graupel collecting rain.
      CALL wrf_debug(200, '  creating rain collecting graupel table')
      CALL wrf_debug(200, '  This may take a few minutes on a sinlge processor...')
      call qr_acr_qg
!..Rain collecting snow & snow collecting rain.
      CALL wrf_debug(200, '  creating rain collecting snow table')
      CALL wrf_debug(200, '  This may take a few minutes on a sinlge processor...')
      call qr_acr_qs

      !bloss: Write out microphysical lookup tables to file
      OPEN(UNIT=92,FILE='thompson_lookup_tables_r1.bin', &
           FORM='unformatted')

      ! write lookup table dimensions to file
      write(92) nbins, &
               nbr, nbs, nbg, &
               ntb_r, ntb_s, ntb_g, &
               ntb_g1, ntb_r1, ntb_t

      ! write microphysical parameters to file, so that we can check for consistency
      !   before deciding that we don't need to recompute the lookup tables.
      write(92) mu_r, am_r, bm_r, av_r, bv_r, D0r, &
             mu_s, Kap0, Kap1, Lam0, Lam1, & ! extra snow parameters
             am_s, bm_s, av_s, bv_s, fv_s, D0s, &
             mu_g, am_g, bm_g, av_g, bv_g, D0g

      WRITE(92) tcg_racg, tmr_racg, tcr_gacr, tmg_gacr, &
           tnr_racg, tnr_gacr ! computed in qr_acr_qg
      WRITE(92) tcs_racs1, tmr_racs1, tcs_racs2, tmr_racs2, &
           tcr_sacr1, tms_sacr1, tcr_sacr2, tms_sacr2, &
           tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2 ! computed in qr_acr_qs
      WRITE(92) tpi_qcfz, tni_qcfz, tpi_qrfz, tpg_qrfz, & 
           tni_qrfz, tnr_qrfz, & !computed in freezeH2O
           tps_iaus, tni_iaus, tpi_ide  ! computed in qi_aut_qs
      CLOSE(92)
      endif ! .NOT.read_lookup_tables_from_file

!..Cloud water and rain freezing (Bigg, 1953).
      CALL wrf_debug(200, '  creating freezing of water drops table')
      call freezeH2O

!..Conversion of some ice mass into snow category.
      CALL wrf_debug(200, '  creating ice converting to snow table')
      call qi_aut_qs

      endif ! .NOT.iiwarm
      CALL wrf_debug(150, ' ... DONE microphysical lookup tables')
      endif ! micro_init
      END SUBROUTINE thompiso_init
!+---+-----------------------------------------------------------------+
!
!+---+-----------------------------------------------------------------+
!..This is a wrapper routine designed to transfer values from 3D to 1D.
!+---+-----------------------------------------------------------------+
!bloss(2017-02-15): Add separate microphysical outputs for
!                  deposition, sublimation, rain evaporation tendencies
!                  rather than single net deposition tendency used previously.
!                  These changes are tagged below with bloss(dep/sub/revp).
      SUBROUTINE mp_gt_driver_iso(qv, qc, qr, qi, qs, qg, &
                              hdo_qv, hdo_qc, hdo_qr, hdo_qi, hdo_qs, hdo_qg, &
                              o18_qv, o18_qc, o18_qr, o18_qi, o18_qs, o18_qg, &
                              ni, nr, &
                              th, pii, p, dz, dt_in, itimestep, &
                              RAINNC, RAINNCV, &
                              SNOWNC, SNOWNCV, &
                              GRAUPELNC, GRAUPELNCV, SR, &
                              HDO_RAINNC, HDO_RAINNCV, &
                              HDO_SNOWNC, HDO_SNOWNCV, &
                              HDO_GRAUPELNC, HDO_GRAUPELNCV, &
                              O18_RAINNC, O18_RAINNCV, &
                              O18_SNOWNC, O18_SNOWNCV, &
                              O18_GRAUPELNC, O18_GRAUPELNCV, &
                          QLIQ_AUTO, QLIQ_COLL, QICE_AUTO, QICE_COLL, &
                          QICE_RIM, QVAP_CLD, QPCP_SED, &
                          QPCP_DEP, QPCP_SUB, QPCP_REVP, & !bloss(dep/sub/revp)
                          HDO_LIQ_AUTO, HDO_LIQ_COLL, HDO_ICE_AUTO, HDO_ICE_COLL, &
                          HDO_ICE_RIM, HDO_VAP_CLD, HDO_PCP_SED, &
                          HDO_PCP_DEP, HDO_PCP_SUB, HDO_PCP_REVP, & !bloss(dep/sub/revp)
                          O18_LIQ_AUTO, O18_LIQ_COLL, O18_ICE_AUTO, O18_ICE_COLL, &
                          O18_ICE_RIM, O18_VAP_CLD, O18_PCP_SED, &
                          O18_PCP_DEP, O18_PCP_SUB, O18_PCP_REVP, & !bloss(dep/sub/revp)
                              refl_10cm, diagflag, do_radar_ref,      &
                              re_cloud, re_ice, re_snow,              &
                              has_reqc, has_reqi, has_reqs,           &
                              ids,ide, jds,jde, kds,kde, &             ! domain dims
                              ims,ime, jms,jme, kms,kme, &             ! memory dims
                              its,ite, jts,jte, kts,kte)               ! tile dims
      implicit none
!..Subroutine arguments
      INTEGER, INTENT(IN):: ids,ide, jds,jde, kds,kde, &
                            ims,ime, jms,jme, kms,kme, &
                            its,ite, jts,jte, kts,kte
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
                          qv, qc, qr, qi, qs, qg, ni, nr, th
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
                          hdo_qv, hdo_qc, hdo_qr, hdo_qi, hdo_qs, hdo_qg, &
                          o18_qv, o18_qc, o18_qr, o18_qi, o18_qs, o18_qg
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
                          re_cloud, re_ice, re_snow
      INTEGER, INTENT(IN):: has_reqc, has_reqi, has_reqs
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
                          pii, p, dz
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT):: &
                          RAINNC, RAINNCV, SR
      REAL, DIMENSION(ims:ime, jms:jme), OPTIONAL, INTENT(INOUT)::      &
                          SNOWNC, SNOWNCV, GRAUPELNC, GRAUPELNCV
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT):: &
                          HDO_RAINNC, HDO_RAINNCV, HDO_SNOWNC, HDO_SNOWNCV, &
                          HDO_GRAUPELNC, HDO_GRAUPELNCV, &
                          O18_RAINNC, O18_RAINNCV, O18_SNOWNC, O18_SNOWNCV, &
                          O18_GRAUPELNC, O18_GRAUPELNCV

      ! Accumulated aggregate microphysical tendencies
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), OPTIONAL, INTENT(INOUT)::      &
                          QLIQ_AUTO, QLIQ_COLL, QICE_AUTO, QICE_COLL, &
                          QICE_RIM, QVAP_CLD, QPCP_SED, &
                          QPCP_DEP, QPCP_SUB, QPCP_REVP, & !bloss(dep/sub/revp)
                          HDO_LIQ_AUTO, HDO_LIQ_COLL, HDO_ICE_AUTO, HDO_ICE_COLL, &
                          HDO_ICE_RIM, HDO_VAP_CLD, HDO_PCP_SED, &
                          HDO_PCP_DEP, HDO_PCP_SUB, HDO_PCP_REVP, & !bloss(dep/sub/revp)
                          O18_LIQ_AUTO, O18_LIQ_COLL, O18_ICE_AUTO, O18_ICE_COLL, &
                          O18_ICE_RIM, O18_VAP_CLD, O18_PCP_SED, &
                          O18_PCP_DEP, O18_PCP_SUB, O18_PCP_REVP !bloss(dep/sub/revp)

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT)::       &
                          refl_10cm
      REAL, INTENT(IN):: dt_in
      INTEGER, INTENT(IN):: itimestep
!..Local variables
      REAL, DIMENSION(kts:kte):: &
                          qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d, &
                          nr1d, nc1d, t1d, p1d, dz1d, dBZ
      REAL, DIMENSION(kts:kte):: re_qc1d, re_qi1d, re_qs1d
      ! isotope masses
      INTEGER, PARAMETER :: niso_max = 3
      INTEGER :: niso, nn
      INTEGER, DIMENSION(niso_max) :: iso_index
      REAL, DIMENSION(kts:kte,niso_max) :: &
           iso_qv1d, iso_qc1d, iso_qi1d, iso_qr1d, iso_qs1d, iso_qg1d
      REAL, DIMENSION(niso_max) :: iso_pptrain, iso_pptsnow, iso_pptgraul, iso_pptice

      ! aggregate tendencies
      REAL, DIMENSION(kts:kte) :: &
                          qliq_auto1d, qliq_coll1d, qice_auto1d, qice_coll1d, &
                          qice_rim1d, qvap_cld1d, qpcp_sed1d, &
                          qpcp_dep1d, qpcp_sub1d, qpcp_revp1d !bloss(dep/sub/revp)
      REAL, DIMENSION(kts:kte,niso_max) :: &
                          iso_qliq_auto1d, iso_qliq_coll1d, iso_qice_auto1d, iso_qice_coll1d, &
                          iso_qice_rim1d, iso_qvap_cld1d, iso_qpcp_sed1d, &
                          iso_qpcp_dep1d, iso_qpcp_sub1d, iso_qpcp_revp1d !bloss(dep/sub/revp)
           
      LOGICAL :: accumulated_hdo, accumulated_o18
      REAL, DIMENSION(its:ite, jts:jte):: pcp_ra, pcp_sn, pcp_gr, pcp_ic
      REAL:: dt, pptrain, pptsnow, pptgraul, pptice
      REAL:: qc_max, qr_max, qs_max, qi_max, qg_max, ni_max, nr_max
      INTEGER:: i, j, k
      INTEGER:: imax_qc,imax_qr,imax_qi,imax_qs,imax_qg,imax_ni,imax_nr
      INTEGER:: jmax_qc,jmax_qr,jmax_qi,jmax_qs,jmax_qg,jmax_ni,jmax_nr
      INTEGER:: kmax_qc,kmax_qr,kmax_qi,kmax_qs,kmax_qg,kmax_ni,kmax_nr
      INTEGER:: i_start, j_start, i_end, j_end
      LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
      INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref
!+---+
      i_start = its
      j_start = jts
      i_end   = MIN(ite, ide-1)
      j_end   = MIN(jte, jde-1)
!..For idealized testing by developer.
!     if ( (ide-ids+1).gt.4 .and. (jde-jds+1).lt.4 .and.                &
!          ids.eq.its.and.ide.eq.ite.and.jds.eq.jts.and.jde.eq.jte) then
!        i_start = its + 2
!        i_end   = ite
!        j_start = jts
!        j_end   = jte
!     endif
      dt = dt_in
   
      qc_max = 0.
      qr_max = 0.
      qs_max = 0.
      qi_max = 0.
      qg_max = 0
      ni_max = 0.
      nr_max = 0.
      imax_qc = 0
      imax_qr = 0
      imax_qi = 0
      imax_qs = 0
      imax_qg = 0
      imax_ni = 0
      imax_nr = 0
      jmax_qc = 0
      jmax_qr = 0
      jmax_qi = 0
      jmax_qs = 0
      jmax_qg = 0
      jmax_ni = 0
      jmax_nr = 0
      kmax_qc = 0
      kmax_qr = 0
      kmax_qi = 0
      kmax_qs = 0
      kmax_qg = 0
      kmax_ni = 0
      kmax_nr = 0
      do i = 1, 256
         mp_debug(i:i) = char(0)
      enddo
      j_loop:  do j = j_start, j_end
      i_loop:  do i = i_start, i_end
         pptrain = 0.
         pptsnow = 0.
         pptgraul = 0.
         pptice = 0.
         RAINNCV(i,j) = 0.
         IF ( PRESENT (snowncv) ) THEN
            SNOWNCV(i,j) = 0.
         ENDIF
         IF ( PRESENT (graupelncv) ) THEN
            GRAUPELNCV(i,j) = 0.
         ENDIF
         SR(i,j) = 0.
         do k = kts, kte
            t1d(k) = th(i,k,j)*pii(i,k,j)
            p1d(k) = p(i,k,j)
            dz1d(k) = dz(i,k,j)
            qv1d(k) = qv(i,k,j)
            qc1d(k) = qc(i,k,j)
            qi1d(k) = qi(i,k,j)
            qr1d(k) = qr(i,k,j)
            qs1d(k) = qs(i,k,j)
            qg1d(k) = qg(i,k,j)
            ni1d(k) = ni(i,k,j)
            nr1d(k) = nr(i,k,j)
         enddo
         niso = 2
         iso_index(1) = 19 ! HDO
         do k = kts, kte
           iso_qv1d(k,1) = hdo_qv(i,k,j) ! dummy initialization !!!TODO!!! FIX LATER
           iso_qc1d(k,1) = hdo_qc(i,k,j)
           iso_qi1d(k,1) = hdo_qi(i,k,j)
           iso_qr1d(k,1) = hdo_qr(i,k,j)
           iso_qs1d(k,1) = hdo_qs(i,k,j)
           iso_qg1d(k,1) = hdo_qg(i,k,j)
         end do
         iso_index(2) = 20 ! H2O18
         do k = kts, kte
           iso_qv1d(k,2) = o18_qv(i,k,j) ! dummy initialization !!!TODO!!! FIX LATER
           iso_qc1d(k,2) = o18_qc(i,k,j)
           iso_qi1d(k,2) = o18_qi(i,k,j)
           iso_qr1d(k,2) = o18_qr(i,k,j)
           iso_qs1d(k,2) = o18_qs(i,k,j)
           iso_qg1d(k,2) = o18_qg(i,k,j)
         end do
         iso_pptrain(:) = 0.
         iso_pptsnow(:) = 0.
         iso_pptgraul(:) = 0.
         iso_pptice(:) = 0.
         call mp_thompiso(qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d, &
                      nr1d, t1d, p1d, dz1d, &
                      pptrain, pptsnow, pptgraul, pptice, &
                          iso_qv1d, iso_qc1d, iso_qi1d, iso_qr1d, iso_qs1d, iso_qg1d, & 
                          iso_pptrain, iso_pptsnow, iso_pptgraul, iso_pptice, &
                          qliq_auto1d, qliq_coll1d, qice_auto1d, qice_coll1d, qice_rim1d, qvap_cld1d, qpcp_sed1d, &
                          qpcp_dep1d, qpcp_sub1d, qpcp_revp1d, & !bloss(dep/sub/revp)
                          iso_qliq_auto1d, iso_qliq_coll1d, iso_qice_auto1d, iso_qice_coll1d, &
                          iso_qice_rim1d, iso_qvap_cld1d, iso_qpcp_sed1d, &
                          iso_qpcp_dep1d, iso_qpcp_sub1d, iso_qpcp_revp1d, & !bloss(dep/sub/revp)
                      kts, kte, dt, i, j &
                          , niso, iso_index &
                          )
         pcp_ra(i,j) = pptrain
         pcp_sn(i,j) = pptsnow
         pcp_gr(i,j) = pptgraul
         pcp_ic(i,j) = pptice
         RAINNCV(i,j) = pptrain + pptsnow + pptgraul + pptice
         RAINNC(i,j) = RAINNC(i,j) + pptrain + pptsnow + pptgraul + pptice
         IF ( PRESENT(snowncv) .AND. PRESENT(snownc) ) THEN
            SNOWNCV(i,j) = pptsnow + pptice
            SNOWNC(i,j) = SNOWNC(i,j) + pptsnow + pptice
         ENDIF
         IF ( PRESENT(graupelncv) .AND. PRESENT(graupelnc) ) THEN
            GRAUPELNCV(i,j) = pptgraul
            GRAUPELNC(i,j) = GRAUPELNC(i,j) + pptgraul
         ENDIF
         SR(i,j) = (pptsnow + pptgraul + pptice)/(RAINNCV(i,j)+1.e-12)
         do k = kts, kte
            nc1d(k) = Nt_c
         enddo
         
         !bloss(2014-10): update aggregate tendencies, mutliply be timestep,
         !  so that tendency becomes increment.  Units are now kg/kg.
         if (PRESENT(qliq_auto)) &
              qliq_auto(i,kts:kte,j) = qliq_auto(i,kts:kte,j) + dt*qliq_auto1d(kts:kte)
         if (PRESENT(qliq_coll)) &
              qliq_coll(i,kts:kte,j) = qliq_coll(i,kts:kte,j) + dt*qliq_coll1d(kts:kte)
         if (PRESENT(qice_auto)) &
              qice_auto(i,kts:kte,j) = qice_auto(i,kts:kte,j) + dt*qice_auto1d(kts:kte)
         if (PRESENT(qice_coll)) &
              qice_coll(i,kts:kte,j) = qice_coll(i,kts:kte,j) + dt*qice_coll1d(kts:kte)
         if (PRESENT(qice_rim)) &
              qice_rim(i,kts:kte,j) = qice_rim(i,kts:kte,j) + dt*qice_rim1d(kts:kte)
         if (PRESENT(qvap_cld)) &
              qvap_cld(i,kts:kte,j) = qvap_cld(i,kts:kte,j) + dt*qvap_cld1d(kts:kte)
         if (PRESENT(qpcp_sed)) &
              qpcp_sed(i,kts:kte,j) = qpcp_sed(i,kts:kte,j) + dt*qpcp_sed1d(kts:kte)
         if (PRESENT(qpcp_dep)) &
              qpcp_dep(i,kts:kte,j) = qpcp_dep(i,kts:kte,j) + dt*qpcp_dep1d(kts:kte) !bloss(dep/sub/revp)
         if (PRESENT(qpcp_sub)) &
              qpcp_sub(i,kts:kte,j) = qpcp_sub(i,kts:kte,j) + dt*qpcp_sub1d(kts:kte) !bloss(dep/sub/revp)
         if (PRESENT(qpcp_revp)) &
              qpcp_revp(i,kts:kte,j) = qpcp_revp(i,kts:kte,j) + dt*qpcp_revp1d(kts:kte) !bloss(dep/sub/revp)

         do k = kts, kte
            qv(i,k,j) = qv1d(k)
            qc(i,k,j) = qc1d(k)
            qi(i,k,j) = qi1d(k)
            qr(i,k,j) = qr1d(k)
            qs(i,k,j) = qs1d(k)
            qg(i,k,j) = qg1d(k)
            ni(i,k,j) = ni1d(k)
            nr(i,k,j) = nr1d(k)
            th(i,k,j) = t1d(k)/pii(i,k,j)

            if (qc1d(k) .gt. qc_max) then
             imax_qc = i
             jmax_qc = j
             kmax_qc = k
             qc_max = qc1d(k)
            elseif (qc1d(k) .lt. 0.0) then
             write(mp_debug,*) 'WARNING, negative qc ', qc1d(k),        &
                        ' at i,j,k=', i,j,k
             CALL wrf_debug(150, mp_debug)
            endif
            if (qr1d(k) .gt. qr_max) then
             imax_qr = i
             jmax_qr = j
             kmax_qr = k
             qr_max = qr1d(k)
            elseif (qr1d(k) .lt. 0.0) then
             write(mp_debug,*) 'WARNING, negative qr ', qr1d(k),        &
                        ' at i,j,k=', i,j,k
             CALL wrf_debug(150, mp_debug)
            endif
            if (nr1d(k) .gt. nr_max) then
             imax_nr = i
             jmax_nr = j
             kmax_nr = k
             nr_max = nr1d(k)
            elseif (nr1d(k) .lt. 0.0) then
             write(mp_debug,*) 'WARNING, negative nr ', nr1d(k),        &
                        ' at i,j,k=', i,j,k
             CALL wrf_debug(150, mp_debug)
            endif
            if (qs1d(k) .gt. qs_max) then
             imax_qs = i
             jmax_qs = j
             kmax_qs = k
             qs_max = qs1d(k)
            elseif (qs1d(k) .lt. 0.0) then
             write(mp_debug,*) 'WARNING, negative qs ', qs1d(k),        &
                        ' at i,j,k=', i,j,k
             CALL wrf_debug(150, mp_debug)
            endif
            if (qi1d(k) .gt. qi_max) then
             imax_qi = i
             jmax_qi = j
             kmax_qi = k
             qi_max = qi1d(k)
            elseif (qi1d(k) .lt. 0.0) then
             write(mp_debug,*) 'WARNING, negative qi ', qi1d(k),        &
                        ' at i,j,k=', i,j,k
             CALL wrf_debug(150, mp_debug)
            endif
            if (qg1d(k) .gt. qg_max) then
             imax_qg = i
             jmax_qg = j
             kmax_qg = k
             qg_max = qg1d(k)
            elseif (qg1d(k) .lt. 0.0) then
             write(mp_debug,*) 'WARNING, negative qg ', qg1d(k),        &
                        ' at i,j,k=', i,j,k
             CALL wrf_debug(150, mp_debug)
            endif
            if (ni1d(k) .gt. ni_max) then
             imax_ni = i
             jmax_ni = j
             kmax_ni = k
             ni_max = ni1d(k)
            elseif (ni1d(k) .lt. 0.0) then
             write(mp_debug,*) 'WARNING, negative ni ', ni1d(k),        &
                        ' at i,j,k=', i,j,k
             CALL wrf_debug(150, mp_debug)
            endif
            if (qv1d(k) .lt. 0.0) then
             write(mp_debug,*) 'WARNING, negative qv ', qv1d(k),        &
                        ' at i,j,k=', i,j,k
             CALL wrf_debug(150, mp_debug)
             if (k.lt.kte-2 .and. k.gt.kts+1) then
                write(mp_debug,*) '   below and above are: ', qv(i,k-1,j), qv(i,k+1,j)
                CALL wrf_debug(150, mp_debug)
                qv(i,k,j) = MAX(1.E-7, 0.5*(qv(i,k-1,j) + qv(i,k+1,j)))
             else
                qv(i,k,j) = 1.E-7
             endif
            endif
         enddo
         IF ( PRESENT (diagflag) ) THEN
         if (diagflag .and. do_radar_ref == 1) then
          call calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d,       &
                      t1d, p1d, dBZ, kts, kte, i, j)
          do k = kts, kte
             refl_10cm(i,k,j) = MAX(-35., dBZ(k))
          enddo
         endif
         ENDIF
         IF ( (has_reqc.ne.0) .and. (has_reqi.ne.0) .and. (has_reqs.ne.0) ) THEN
          do k = kts, kte
             re_qc1d(k) = 2.51E-6
             re_qi1d(k) = 5.01E-6
             re_qs1d(k) = 25.01E-6
          enddo
          call calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d,  &
                      re_qc1d, re_qi1d, re_qs1d, kts, kte)
          do k = kts, kte
             re_cloud(i,k,j) = MAX(2.51E-6, MIN(re_qc1d(k), 50.E-6))
             re_ice(i,k,j)   = MAX(5.01E-6, MIN(re_qi1d(k), 125.E-6))
             re_snow(i,k,j)  = MAX(25.E-6, MIN(re_qs1d(k), 999.E-6))
          enddo
         ENDIF
         accumulated_hdo = .false.
         accumulated_o18 = .false.
         do nn = 1,niso
           SELECT CASE(iso_index(nn))
           CASE(19)
             if(.NOT.accumulated_hdo) then
               HDO_RAINNCV(i,j) = iso_pptrain(nn) &
                    + iso_pptsnow(nn) + iso_pptgraul(nn) + iso_pptice(nn)
               HDO_RAINNC(i,j) = HDO_RAINNC(i,j) + iso_pptrain(nn) &
                    + iso_pptsnow(nn) + iso_pptgraul(nn) + iso_pptice(nn)
               HDO_SNOWNCV(i,j) = iso_pptsnow(nn) + iso_pptice(nn)
               HDO_SNOWNC(i,j) = HDO_SNOWNC(i,j) + iso_pptsnow(nn) + iso_pptice(nn)
               HDO_GRAUPELNCV(i,j) = iso_pptgraul(nn)
               HDO_GRAUPELNC(i,j) = HDO_GRAUPELNC(i,j) + iso_pptgraul(nn)
               do k = kts,kte
                 hdo_qv(i,k,j) = iso_qv1d(k,nn)
                 hdo_qc(i,k,j) = iso_qc1d(k,nn)
                 hdo_qi(i,k,j) = iso_qi1d(k,nn)
                 hdo_qr(i,k,j) = iso_qr1d(k,nn)
                 hdo_qs(i,k,j) = iso_qs1d(k,nn)
                 hdo_qg(i,k,j) = iso_qg1d(k,nn)
               end do
               accumulated_hdo = .true.
             end if

             ! accumulate aggregate microphysical HDO tendencies, multiply by dt
             !  so that they become increments.
             if (PRESENT(hdo_liq_auto)) &
                  hdo_liq_auto(i,kts:kte,j) = hdo_liq_auto(i,kts:kte,j) + dt*iso_qliq_auto1d(kts:kte,nn)
             if (PRESENT(hdo_liq_coll)) &
                  hdo_liq_coll(i,kts:kte,j) = hdo_liq_coll(i,kts:kte,j) + dt*iso_qliq_coll1d(kts:kte,nn)
             if (PRESENT(hdo_ice_auto)) &
                  hdo_ice_auto(i,kts:kte,j) = hdo_ice_auto(i,kts:kte,j) + dt*iso_qice_auto1d(kts:kte,nn)
             if (PRESENT(hdo_ice_coll)) &
                  hdo_ice_coll(i,kts:kte,j) = hdo_ice_coll(i,kts:kte,j) + dt*iso_qice_coll1d(kts:kte,nn)
             if (PRESENT(hdo_ice_rim)) &
                  hdo_ice_rim(i,kts:kte,j) = hdo_ice_rim(i,kts:kte,j) + dt*iso_qice_rim1d(kts:kte,nn)
             if (PRESENT(hdo_vap_cld)) &
                  hdo_vap_cld(i,kts:kte,j) = hdo_vap_cld(i,kts:kte,j) + dt*iso_qvap_cld1d(kts:kte,nn)
             if (PRESENT(hdo_pcp_sed)) &
                  hdo_pcp_sed(i,kts:kte,j) = hdo_pcp_sed(i,kts:kte,j) + dt*iso_qpcp_sed1d(kts:kte,nn)
             if (PRESENT(hdo_pcp_dep)) &
                  hdo_pcp_dep(i,kts:kte,j) = hdo_pcp_dep(i,kts:kte,j) + dt*iso_qpcp_dep1d(kts:kte,nn) !bloss(dep/sub/revp)
             if (PRESENT(hdo_pcp_sub)) &
                  hdo_pcp_sub(i,kts:kte,j) = hdo_pcp_sub(i,kts:kte,j) + dt*iso_qpcp_sub1d(kts:kte,nn) !bloss(dep/sub/revp)
             if (PRESENT(hdo_pcp_revp)) &
                  hdo_pcp_revp(i,kts:kte,j) = hdo_pcp_revp(i,kts:kte,j) + dt*iso_qpcp_revp1d(kts:kte,nn) !bloss(dep/sub/revp)

           CASE(20)
             if(.NOT.accumulated_o18) then
               O18_RAINNCV(i,j) = iso_pptrain(nn) &
                    + iso_pptsnow(nn) + iso_pptgraul(nn) + iso_pptice(nn)
               O18_RAINNC(i,j) = O18_RAINNC(i,j) + iso_pptrain(nn) &
                    + iso_pptsnow(nn) + iso_pptgraul(nn) + iso_pptice(nn)
               O18_SNOWNCV(i,j) = iso_pptsnow(nn) + iso_pptice(nn)
               O18_SNOWNC(i,j) = O18_SNOWNC(i,j) + iso_pptsnow(nn) + iso_pptice(nn)
               O18_GRAUPELNCV(i,j) = iso_pptgraul(nn)
               O18_GRAUPELNC(i,j) = O18_GRAUPELNC(i,j) + iso_pptgraul(nn)
               do k = kts,kte
                 o18_qv(i,k,j) = iso_qv1d(k,nn)
                 o18_qc(i,k,j) = iso_qc1d(k,nn)
                 o18_qi(i,k,j) = iso_qi1d(k,nn)
                 o18_qr(i,k,j) = iso_qr1d(k,nn)
                 o18_qs(i,k,j) = iso_qs1d(k,nn)
                 o18_qg(i,k,j) = iso_qg1d(k,nn)
               end do
               accumulated_o18 = .true.
             end if
               
             ! accumulate aggregate microphysicalH2 O18 tendencies, multiply by dt
             !  so that they become increments.
             if (PRESENT(o18_liq_auto)) &
                  o18_liq_auto(i,kts:kte,j) = o18_liq_auto(i,kts:kte,j) + dt*iso_qliq_auto1d(kts:kte,nn)
             if (PRESENT(o18_liq_coll)) &
                  o18_liq_coll(i,kts:kte,j) = o18_liq_coll(i,kts:kte,j) + dt*iso_qliq_coll1d(kts:kte,nn)
             if (PRESENT(o18_ice_auto)) &
                  o18_ice_auto(i,kts:kte,j) = o18_ice_auto(i,kts:kte,j) + dt*iso_qice_auto1d(kts:kte,nn)
             if (PRESENT(o18_ice_coll)) &
                  o18_ice_coll(i,kts:kte,j) = o18_ice_coll(i,kts:kte,j) + dt*iso_qice_coll1d(kts:kte,nn)
             if (PRESENT(o18_ice_rim)) &
                  o18_ice_rim(i,kts:kte,j) = o18_ice_rim(i,kts:kte,j) + dt*iso_qice_rim1d(kts:kte,nn)
             if (PRESENT(o18_vap_cld)) &
                  o18_vap_cld(i,kts:kte,j) = o18_vap_cld(i,kts:kte,j) + dt*iso_qvap_cld1d(kts:kte,nn)
             if (PRESENT(o18_pcp_sed)) &
                  o18_pcp_sed(i,kts:kte,j) = o18_pcp_sed(i,kts:kte,j) + dt*iso_qpcp_sed1d(kts:kte,nn)
             if (PRESENT(o18_pcp_dep)) &
                  o18_pcp_dep(i,kts:kte,j) = o18_pcp_dep(i,kts:kte,j) + dt*iso_qpcp_dep1d(kts:kte,nn) !bloss(dep/sub/revp)
             if (PRESENT(o18_pcp_sub)) &
                  o18_pcp_sub(i,kts:kte,j) = o18_pcp_sub(i,kts:kte,j) + dt*iso_qpcp_sub1d(kts:kte,nn) !bloss(dep/sub/revp)
             if (PRESENT(o18_pcp_revp)) &
                  o18_pcp_revp(i,kts:kte,j) = o18_pcp_revp(i,kts:kte,j) + dt*iso_qpcp_revp1d(kts:kte,nn) !bloss(dep/sub/revp)

           CASE DEFAULT
             write(*,*) 'BAD ISO_INDEX'
             STOP 'IN module_mp_thompiso'
           END SELECT
         END DO
      enddo i_loop
      enddo j_loop
! DEBUG - GT
      write(mp_debug,'(a,7(a,e13.6,1x,a,i3,a,i3,a,i3,a,1x))') 'MP-GT:', &
         'qc: ', qc_max, '(', imax_qc, ',', jmax_qc, ',', kmax_qc, ')', &
         'qr: ', qr_max, '(', imax_qr, ',', jmax_qr, ',', kmax_qr, ')', &
         'qi: ', qi_max, '(', imax_qi, ',', jmax_qi, ',', kmax_qi, ')', &
         'qs: ', qs_max, '(', imax_qs, ',', jmax_qs, ',', kmax_qs, ')', &
         'qg: ', qg_max, '(', imax_qg, ',', jmax_qg, ',', kmax_qg, ')', &
         'ni: ', ni_max, '(', imax_ni, ',', jmax_ni, ',', kmax_ni, ')', &
         'nr: ', nr_max, '(', imax_nr, ',', jmax_nr, ',', kmax_nr, ')'
      CALL wrf_debug(150, mp_debug)
!!$      write(mp_debug,'(a,6(a,e13.6,1x,a,i3,a,i3,a,i3,a,1x))') 'MP-GT-ISO HDO MAX:', &
!!$         'hdo_qv: ', hdo_qv_max, '(', imax_hdo_qv, ',', jmax_hdo_qv, ',', kmax_hdo_qv, ')', &
!!$         'hdo_qc: ', hdo_qc_max, '(', imax_hdo_qc, ',', jmax_hdo_qc, ',', kmax_hdo_qc, ')', &
!!$         'hdo_qr: ', hdo_qr_max, '(', imax_hdo_qr, ',', jmax_hdo_qr, ',', kmax_hdo_qr, ')', &
!!$         'hdo_qi: ', hdo_qi_max, '(', imax_hdo_qi, ',', jmax_hdo_qi, ',', kmax_hdo_qi, ')', &
!!$         'hdo_qs: ', hdo_qs_max, '(', imax_hdo_qs, ',', jmax_hdo_qs, ',', kmax_hdo_qs, ')', &
!!$         'hdo_qg: ', hdo_qg_max, '(', imax_hdo_qg, ',', jmax_hdo_qg, ',', kmax_hdo_qg, ')'
!!$      CALL wrf_debug(150, mp_debug)
!!$      write(mp_debug,'(a,6(a,e13.6,1x,a,i3,a,i3,a,i3,a,1x))') 'MP-GT-ISO O18 MAX:', &
!!$         'o18_qv: ', o18_qv_max, '(', imax_o18_qv, ',', jmax_o18_qv, ',', kmax_o18_qv, ')', &
!!$         'o18_qc: ', o18_qc_max, '(', imax_o18_qc, ',', jmax_o18_qc, ',', kmax_o18_qc, ')', &
!!$         'o18_qr: ', o18_qr_max, '(', imax_o18_qr, ',', jmax_o18_qr, ',', kmax_o18_qr, ')', &
!!$         'o18_qi: ', o18_qi_max, '(', imax_o18_qi, ',', jmax_o18_qi, ',', kmax_o18_qi, ')', &
!!$         'o18_qs: ', o18_qs_max, '(', imax_o18_qs, ',', jmax_o18_qs, ',', kmax_o18_qs, ')', &
!!$         'o18_qg: ', o18_qg_max, '(', imax_o18_qg, ',', jmax_o18_qg, ',', kmax_o18_qg, ')'
!!$      CALL wrf_debug(150, mp_debug)
! END DEBUG - GT
      do i = 1, 256
         mp_debug(i:i) = char(0)
      enddo
      END SUBROUTINE mp_gt_driver_iso
!+---+-----------------------------------------------------------------+
!
!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
!.. This subroutine computes the moisture tendencies of water vapor,
!.. cloud droplets, rain, cloud ice (pristine), snow, and graupel.
!.. Previously this code was based on Reisner et al (1998), but few of
!.. those pieces remain.  A complete description is now found in
!.. Thompson et al. (2004, 2008).
!+---+-----------------------------------------------------------------+
!
      subroutine mp_thompiso (qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d, &
                          nr1d, t1d, p1d, dzq, &
                          pptrain, pptsnow, pptgraul, pptice, &
                          iso_qv1d, iso_qc1d, iso_qi1d, iso_qr1d, iso_qs1d, iso_qg1d, & 
                          iso_pptrain, iso_pptsnow, iso_pptgraul, iso_pptice, &
                          qliq_auto, qliq_coll, qice_auto, qice_coll, qice_rim, qvap_cld, qpcp_sed, &
                          qpcp_dep, qpcp_sub, qpcp_revp, & !bloss(dep/sub/revp)
                          iso_qliq_auto, iso_qliq_coll, iso_qice_auto, iso_qice_coll, &
                          iso_qice_rim, iso_qvap_cld, iso_qpcp_sed, &
                          iso_qpcp_dep, iso_qpcp_sub, iso_qpcp_revp, & !bloss(dep/sub/revp)
                          kts, kte, dt, ii, jj &
                          , niso, iso_index &
                          )
      implicit none
!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte, ii, jj
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: &
                          qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d, &
                          nr1d, t1d, p1d
      REAL, DIMENSION(kts:kte), INTENT(IN):: dzq
      REAL, INTENT(INOUT):: pptrain, pptsnow, pptgraul, pptice
      REAL, INTENT(IN):: dt
      LOGICAL, PARAMETER :: dosedimentation = .true.
      ! isotope masses
      INTEGER, INTENT(in) :: niso
      INTEGER, DIMENSION(niso), INTENT(in) :: iso_index
      REAL, DIMENSION(kts:kte,niso), INTENT(INOUT) :: &
           iso_qv1d, iso_qc1d, iso_qi1d, iso_qr1d, iso_qs1d, iso_qg1d
      REAL, DIMENSION(niso), INTENT(INOUT):: &
           iso_pptrain, iso_pptsnow, iso_pptgraul, iso_pptice
      ! aggregate tendencies
      REAL, DIMENSION(kts:kte), INTENT(OUT) :: &
           qliq_auto, qliq_coll, qice_auto, qice_coll, qice_rim, qvap_cld, qpcp_sed, &
           qpcp_dep, qpcp_sub, qpcp_revp !bloss(dep/sub/revp)
      REAL, DIMENSION(kts:kte,niso), INTENT(OUT) :: &
           iso_qliq_auto, iso_qliq_coll, iso_qice_auto, iso_qice_coll, &
           iso_qice_rim, iso_qvap_cld, iso_qpcp_sed, &
           iso_qpcp_dep, iso_qpcp_sub, iso_qpcp_revp !bloss(dep/sub/revp)

!..Local variables
!bloss: Note on units.
!  q*ten is in kg/kg/s.
!  all mass processes (e.g., prr_rcg) are in kg/m3/s
!  EXCEPT FOR prw_vcd, which is in kg/kg/s.
!  The individual process rates are generally divided by 
!  rho when they are put together into q*ten.
      REAL, DIMENSION(kts:kte):: tten, qvten, qcten, qiten, &
           qrten, qsten, qgten, niten, nrten
      DOUBLE PRECISION, DIMENSION(kts:kte):: prw_vcd ! kg/kg/s <-- ONLY THIS ONE!!
      DOUBLE PRECISION, DIMENSION(kts:kte):: prr_wau, prr_rcw, prr_rcs, &
           prr_rcg, prr_sml, prr_gml, &
           prr_rci, prv_rev,          &  ! mass rates in kg/m3/s.
           pnr_wau, pnr_rcs, pnr_rcg, &
           pnr_rci, pnr_sml, pnr_gml, &
           pnr_rev, pnr_rcr, pnr_rfz
      DOUBLE PRECISION, DIMENSION(kts:kte):: pri_inu, pni_inu, pri_ihm, &
           pni_ihm, pri_wfz, pni_wfz, &
           pri_rfz, pni_rfz, pri_ide, &
           pni_ide, pri_rci, pni_rci, &
           pni_sci, pni_iau
      DOUBLE PRECISION, DIMENSION(kts:kte):: prs_iau, prs_sci, prs_rcs, &
           prs_scw, prs_sde, prs_ihm, &
           prs_ide
      DOUBLE PRECISION, DIMENSION(kts:kte):: prg_scw, prg_rfz, prg_gde, &
           prg_gcw, prg_rci, prg_rcs, prg_scr, &
           prg_rcg, prg_ihm
      DOUBLE PRECISION, PARAMETER:: zeroD0 = 0.0d0
      REAL, DIMENSION(kts:kte):: temp, pres, qv
      REAL, DIMENSION(kts:kte):: rc, ri, rr, rs, rg, ni, nr
      REAL, DIMENSION(kts:kte):: rho, rhof, rhof2
      REAL, DIMENSION(kts:kte):: qvs, qvsi, delQvs
      REAL, DIMENSION(kts:kte):: satw, sati, ssatw, ssati
      REAL, DIMENSION(kts:kte):: diffu, visco, vsc2, &
           tcond, lvap, ocp, lvt2
      DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilamg, N0_r, N0_g
      REAL, DIMENSION(kts:kte):: mvd_r, mvd_c
      REAL, DIMENSION(kts:kte):: smob, smo2, smo1, smo0, &
           smoc, smod, smoe, smof
      REAL, DIMENSION(kts:kte):: sed_r, sed_s, sed_g, sed_i, sed_n
      REAL:: rgvm, delta_tp, orho, lfus2
      REAL, DIMENSION(4):: onstep
      DOUBLE PRECISION:: N0_exp, N0_min, lam_exp, lamc, lamr, lamg
      DOUBLE PRECISION:: lami, ilami
      REAL:: xDc, Dc_b, Dc_g, xDi, xDr, xDs, xDg, Ds_m, Dg_m
      DOUBLE PRECISION:: Dr_star
      REAL:: zeta1, zeta, taud, tau
      REAL:: stoke_r, stoke_s, stoke_g, stoke_i
      REAL:: vti, vtr, vts, vtg
      REAL, DIMENSION(kts:kte+1):: vtik, vtnik, vtrk, vtnrk, vtsk, vtgk
      REAL, DIMENSION(kts:kte):: vts_boost
      REAL:: Mrat, ils1, ils2, t1_vts, t2_vts, t3_vts, t4_vts, C_snow
      REAL:: a_, b_, loga_, A1, A2, tf
      REAL:: tempc, tc0, r_mvd1, r_mvd2, xkrat
      REAL:: xnc, xri, xni, xmi, oxmi, xrc, xrr, xnr
      REAL:: xsat, rate_max, sump, ratio
      REAL:: clap, fcd, dfcd
      REAL:: otemp, rvs, rvs_p, rvs_pp, gamsc, alphsc, t1_evap, t1_subl
      REAL:: r_frac, g_frac
      REAL:: Ef_rw, Ef_sw, Ef_gw, Ef_rr
      REAL:: dtsave, odts, odt, odzq
      REAL:: xslw1, ygra1, zans1, eva_factor
      REAL :: updated_temp !bloss
      INTEGER:: i, k, k2, n, nn, nstep, k_0, kbot, IT, iexfrq
      INTEGER, DIMENSION(4):: ksed1
      INTEGER:: nir, nis, nig, nii, nic
      INTEGER:: idx_tc, idx_t, idx_s, idx_g1, idx_g, idx_r1, idx_r,     &
                idx_i1, idx_i, idx_c, idx, idx_d
      LOGICAL:: melti, no_micro
      LOGICAL, DIMENSION(kts:kte):: L_qc, L_qi, L_qr, L_qs, L_qg
      LOGICAL:: debug_flag
      !bloss: Duplicate a bunch of processes and tendencies for isotopologues
      REAL, DIMENSION(kts:kte,niso):: &
           iso_qvten, iso_qcten, iso_qiten, iso_qrten, iso_qsten, iso_qgten
      DOUBLE PRECISION, DIMENSION(kts:kte,niso):: &
           iso_prw_vcd, & ! vapor deposition
           iso_prr_wau, iso_prr_rcw, iso_prr_rcs, iso_prr_rcg, & ! rain processes
           iso_prr_sml, iso_prr_gml, iso_prr_rci, iso_prv_rev, &
           iso_pri_inu, iso_pri_ihm, &
           iso_pri_wfz, iso_pri_rfz, iso_pri_ide, iso_pri_rci, & ! cloud ice proc.
           iso_prs_iau, iso_prs_sci, iso_prs_rcs, &
           iso_prs_scw, iso_prs_sde, iso_prs_ihm, iso_prs_ide, & ! snow processes
           iso_prg_scw, iso_prg_rfz, iso_prg_gde, & ! graupel processes
           iso_prg_gcw, iso_prg_rci, iso_prg_rcs, iso_prg_scr, &
           iso_prg_rcg, iso_prg_ihm
      ! isotope ratios (e.g., isor_v = iso_qv/qv)
      REAL, DIMENSION(kts:kte,niso):: &
           isor_v, isor_c, isor_i, isor_r, isor_s, isor_g
      ! isotope densities (e.g, rho*iso_qv)
      REAL, DIMENSION(kts:kte,niso):: &
           iso_rv, iso_rc, iso_ri, iso_rr, iso_rs, iso_rg
      ! Updated isotopic mass mixing ratio of vapor
      REAL, DIMENSION(kts:kte,niso):: iso_qv
      REAL, DIMENSION(kts:kte,niso):: iso_sed_r, iso_sed_s, iso_sed_g, iso_sed_i
      REAL :: rv_ice_surface, Sice_tilde, tabs_ice_surface, &
           ventilation_factor_heavy, ventilation_factor_light, &
           ventilation_ratio, alpha_kinetic, &
           iso_qc_new, qcond, qvapor, iso_qtot, &
           iso_rate_max, diffu_heavy, &
           rv_rain_surface, Sliq_tilde, tabs_rain_surface, &
           iso_qv_equil, iso_qr_equil
      REAL, DIMENSION(niso) :: alpha_equil_liq, alpha_equil_rain, &
           alpha_equil_ice, iso_Sc3

      !bloss: turns off isotopic fractionation in the microphysics
      !        !!!THIS IS INTENDED MAINLY FOR TESTING PURPOSES!!!
      LOGICAL, PARAMETER :: fractionate = .true. 
!+---+
      debug_flag = .false.
!     if (ii.eq.315 .and. jj.eq.2) debug_flag = .true.
      no_micro = .true.
      dtsave = dt
      odt = 1./dt
      odts = 1./dtsave
      iexfrq = 1
!+---+-----------------------------------------------------------------+
!.. Source/sink terms.  First 2 chars: "pr" represents source/sink of
!.. mass while "pn" represents source/sink of number.  Next char is one
!.. of "v" for water vapor, "r" for rain, "i" for cloud ice, "w" for
!.. cloud water, "s" for snow, and "g" for graupel.  Next chars
!.. represent processes: "de" for sublimation/deposition, "ev" for
!.. evaporation, "fz" for freezing, "ml" for melting, "au" for
!.. autoconversion, "nu" for ice nucleation, "hm" for Hallet/Mossop
!.. secondary ice production, and "c" for collection followed by the
!.. character for the species being collected.  ALL of these terms are
!.. positive (except for deposition/sublimation terms which can switch
!.. signs based on super/subsaturation) and are treated as negatives
!.. where necessary in the tendency equations.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         tten(k) = 0.
         qvten(k) = 0.
         qcten(k) = 0.
         qiten(k) = 0.
         qrten(k) = 0.
         qsten(k) = 0.
         qgten(k) = 0.
         niten(k) = 0.
         nrten(k) = 0.
         prw_vcd(k) = 0.
         prv_rev(k) = 0.
         prr_wau(k) = 0.
         prr_rcw(k) = 0.
         prr_rcs(k) = 0.
         prr_rcg(k) = 0.
         prr_sml(k) = 0.
         prr_gml(k) = 0.
         prr_rci(k) = 0.
         pnr_wau(k) = 0.
         pnr_rcs(k) = 0.
         pnr_rcg(k) = 0.
         pnr_rci(k) = 0.
         pnr_sml(k) = 0.
         pnr_gml(k) = 0.
         pnr_rev(k) = 0.
         pnr_rcr(k) = 0.
         pnr_rfz(k) = 0.
         pri_inu(k) = 0.
         pni_inu(k) = 0.
         pri_ihm(k) = 0.
         pni_ihm(k) = 0.
         pri_wfz(k) = 0.
         pni_wfz(k) = 0.
         pri_rfz(k) = 0.
         pni_rfz(k) = 0.
         pri_ide(k) = 0.
         pni_ide(k) = 0.
         pri_rci(k) = 0.
         pni_rci(k) = 0.
         pni_sci(k) = 0.
         pni_iau(k) = 0.
         prs_iau(k) = 0.
         prs_sci(k) = 0.
         prs_rcs(k) = 0.
         prs_scw(k) = 0.
         prs_sde(k) = 0.
         prs_ihm(k) = 0.
         prs_ide(k) = 0.
         prg_scw(k) = 0.
         prg_rfz(k) = 0.
         prg_gde(k) = 0.
         prg_gcw(k) = 0.
!bloss         prg_rci(k) = 0.
         prg_rcs(k) = 0.
         prg_scr(k) = 0.
         prg_rcg(k) = 0.
         prg_ihm(k) = 0.
      enddo
      do nn = 1,niso
        do k = kts, kte
          iso_qvten(k,nn) = 0.
          iso_qcten(k,nn) = 0.
          iso_qiten(k,nn) = 0.
          iso_qrten(k,nn) = 0.
          iso_qsten(k,nn) = 0.
          iso_qgten(k,nn) = 0.
          iso_prw_vcd(k,nn) = 0.
          iso_prv_rev(k,nn) = 0.
          iso_prr_wau(k,nn) = 0.
          iso_prr_rcw(k,nn) = 0.
          iso_prr_rcs(k,nn) = 0.
          iso_prr_rcg(k,nn) = 0.
          iso_prr_sml(k,nn) = 0.
          iso_prr_gml(k,nn) = 0.
          iso_prr_rci(k,nn) = 0.
          iso_pri_inu(k,nn) = 0.
          iso_pri_ihm(k,nn) = 0.
          iso_pri_wfz(k,nn) = 0.
          iso_pri_rfz(k,nn) = 0.
          iso_pri_ide(k,nn) = 0.
          iso_pri_rci(k,nn) = 0.
          iso_prs_iau(k,nn) = 0.
          iso_prs_sci(k,nn) = 0.
          iso_prs_rcs(k,nn) = 0.
          iso_prs_scw(k,nn) = 0.
          iso_prs_sde(k,nn) = 0.
          iso_prs_ihm(k,nn) = 0.
          iso_prs_ide(k,nn) = 0.
          iso_prg_scw(k,nn) = 0.
          iso_prg_rfz(k,nn) = 0.
          iso_prg_gde(k,nn) = 0.
          iso_prg_gcw(k,nn) = 0.
!bloss          iso_prg_rci(k,nn) = 0.
          iso_prg_rcs(k,nn) = 0.
          iso_prg_scr(k,nn) = 0.
          iso_prg_rcg(k,nn) = 0.
          iso_prg_ihm(k,nn) = 0.
          isor_v(k,nn) = 0.
          isor_c(k,nn) = 0.
          isor_i(k,nn) = 0.
          isor_r(k,nn) = 0.
          isor_s(k,nn) = 0.
          isor_g(k,nn) = 0.
        enddo
      enddo
!..Schmidt number to one-third used numerous times.
      do nn = 1,niso
        iso_Sc3(nn) = (Sc*Drat_light_over_heavy(iso_index(nn)))**(1./3.)
      end do

      !bloss(2014-10): Initialize aggregate tendencies to zero
      qliq_auto(:) = 0.
      qliq_coll(:) = 0.
      qice_auto(:) = 0.
      qice_coll(:) = 0.
      qice_rim(:) = 0.
      qvap_cld(:) = 0.
      qpcp_sed(:) = 0.
      qpcp_dep(:) = 0. !bloss(dep/sub/revp)
      qpcp_sub(:) = 0. !bloss(dep/sub/revp)
      qpcp_revp(:) = 0. !bloss(dep/sub/revp)
      iso_qliq_coll(:,:) = 0.
      iso_qice_auto(:,:) = 0.
      iso_qice_coll(:,:) = 0.
      iso_qice_rim(:,:) = 0.
      iso_qvap_cld(:,:) = 0.
      iso_qpcp_sed(:,:) = 0.
      iso_qpcp_dep(:,:) = 0. !bloss(dep/sub/revp)
      iso_qpcp_sub(:,:) = 0. !bloss(dep/sub/revp)
      iso_qpcp_revp(:,:) = 0. !bloss(dep/sub/revp)

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         if (qc1d(k) .gt. R1) then
            no_micro = .false.
            rc(k) = qc1d(k)*rho(k)
            L_qc(k) = .true.
         else
            qc1d(k) = 0.0
            rc(k) = R1
            L_qc(k) = .false.
         endif
         if (qi1d(k) .gt. R1) then
            no_micro = .false.
            ri(k) = qi1d(k)*rho(k)
            ni(k) = MAX(R2, ni1d(k)*rho(k))
            L_qi(k) = .true.
            lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
            ilami = 1./lami
            xDi = (bm_i + mu_i + 1.) * ilami
            if (xDi.lt. 20.E-6) then
             lami = cie(2)/20.E-6
             ni(k) = MIN(250.D3, cig(1)*oig2*ri(k)/am_i*lami**bm_i)
            elseif (xDi.gt. 300.E-6) then
             lami = cie(2)/300.E-6
             ni(k) = cig(1)*oig2*ri(k)/am_i*lami**bm_i
            endif
         else
            qi1d(k) = 0.0
            ni1d(k) = 0.0
            ri(k) = R1
            ni(k) = R2
            L_qi(k) = .false.
         endif
         if (qr1d(k) .gt. R1) then
            no_micro = .false.
            rr(k) = qr1d(k)*rho(k)
            nr(k) = MAX(R2, nr1d(k)*rho(k))
            L_qr(k) = .true.
            lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
            mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
            if (mvd_r(k) .gt. 2.5E-3) then
               mvd_r(k) = 2.5E-3
               lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
               nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
            elseif (mvd_r(k) .lt. D0r*0.75) then
               mvd_r(k) = D0r*0.75
               lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
               nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
            endif
         else
            qr1d(k) = 0.0
            nr1d(k) = 0.0
            rr(k) = R1
            nr(k) = R2
            L_qr(k) = .false.
         endif
         if (qs1d(k) .gt. R1) then
            no_micro = .false.
            rs(k) = qs1d(k)*rho(k)
            L_qs(k) = .true.
         else
            qs1d(k) = 0.0
            rs(k) = R1
            L_qs(k) = .false.
         endif
         if (qg1d(k) .gt. R1) then
            no_micro = .false.
            rg(k) = qg1d(k)*rho(k)
            L_qg(k) = .true.
         else
            qg1d(k) = 0.0
            rg(k) = R1
            L_qg(k) = .false.
         endif
      enddo
      ! initialize to unreasonable value.  This should make obvious
      ! any spurious use of this ratio (e.g., multiplying the cloud
      ! liquid isotope ratio by a process rate when no cloud liquid exists)
      isor_v(:,:) = 1.e8
      isor_c(:,:) = 1.e10
      isor_i(:,:) = 1.e12
      isor_r(:,:) = 1.e14
      isor_s(:,:) = 1.e16
      isor_g(:,:) = 1.e18
      do nn = 1,niso
        do k = kts, kte
          isor_v(k,nn) = iso_qv1d(k,nn)/qv(k) ! updated
          if(L_qc(k)) isor_c(k,nn) = iso_qc1d(k,nn)/qc1d(k)
          if(L_qi(k)) isor_i(k,nn) = iso_qi1d(k,nn)/qi1d(k)
          if(L_qr(k)) isor_r(k,nn) = iso_qr1d(k,nn)/qr1d(k)
          if(L_qs(k)) isor_s(k,nn) = iso_qs1d(k,nn)/qs1d(k)
          if(L_qg(k)) isor_g(k,nn) = iso_qg1d(k,nn)/qg1d(k)
        end do
      end do
!+---+-----------------------------------------------------------------+
!..Derive various thermodynamic variables frequently used.
!.. Saturation vapor pressure (mixing ratio) over liquid/ice comes from
!.. Flatau et al. 1992; enthalpy (latent heat) of vaporization from
!.. Bohren & Albrecht 1998; others from Pruppacher & Klett 1978.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         tempc = temp(k) - 273.15
         rhof(k) = SQRT(RHO_NOT/rho(k))
         rhof2(k) = SQRT(rhof(k))
         qvs(k) = rslf(pres(k), temp(k))
         delQvs(k) = MAX(0.0, rslf(pres(k), 273.15)-qv(k))
         if (tempc .le. 0.0) then
          qvsi(k) = rsif(pres(k), temp(k))
         else
          qvsi(k) = qvs(k)
         endif
         satw(k) = qv(k)/qvs(k)
         sati(k) = qv(k)/qvsi(k)
         ssatw(k) = satw(k) - 1.
         ssati(k) = sati(k) - 1.
         if (abs(ssatw(k)).lt. eps) ssatw(k) = 0.0
         if (abs(ssati(k)).lt. eps) ssati(k) = 0.0
         if (no_micro .and. ssati(k).gt. 0.0) no_micro = .false.
         diffu(k) = 2.11E-5*(temp(k)/273.15)**1.94 * (101325./pres(k))
         if (tempc .ge. 0.0) then
            visco(k) = (1.718+0.0049*tempc)*1.0E-5
         else
            visco(k) = (1.718+0.0049*tempc-1.2E-5*tempc*tempc)*1.0E-5
         endif
         ocp(k) = 1./(Cp*(1.+0.887*qv(k)))
         vsc2(k) = SQRT(rho(k)/visco(k))
         lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
         tcond(k) = (5.69 + 0.0168*tempc)*1.0E-5 * 418.936
      enddo
!+---+-----------------------------------------------------------------+
!..If no existing hydrometeor species and no chance to initiate ice or
!.. condense cloud water, just exit quickly!
!+---+-----------------------------------------------------------------+
      if (no_micro) return
!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope, and useful moments for snow.
!+---+-----------------------------------------------------------------+
      if (.not. iiwarm) then
      do k = kts, kte
         if (.not. L_qs(k)) CYCLE
         tc0 = MIN(-0.1, temp(k)-273.15)
         smob(k) = rs(k)*oams
!..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
!.. then we must compute actual 2nd moment and use as reference.
         if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
            smo2(k) = smob(k)
         else
            loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
               + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
               + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
               + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
               + sa(10)*bm_s*bm_s*bm_s
            a_ = 10.0**loga_
            b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
               + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
               + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
               + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
               + sb(10)*bm_s*bm_s*bm_s
            smo2(k) = (smob(k)/a_)**(1./b_)
         endif
!..Calculate 0th moment.  Represents snow number concentration.
         loga_ = sa(1) + sa(2)*tc0 + sa(5)*tc0*tc0 + sa(9)*tc0*tc0*tc0
         a_ = 10.0**loga_
         b_ = sb(1) + sb(2)*tc0 + sb(5)*tc0*tc0 + sb(9)*tc0*tc0*tc0
         smo0(k) = a_ * smo2(k)**b_
!..Calculate 1st moment.  Useful for depositional growth and melting.
         loga_ = sa(1) + sa(2)*tc0 + sa(3) &
               + sa(4)*tc0 + sa(5)*tc0*tc0 &
               + sa(6) + sa(7)*tc0*tc0 &
               + sa(8)*tc0 + sa(9)*tc0*tc0*tc0 &
               + sa(10)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3) + sb(4)*tc0 &
              + sb(5)*tc0*tc0 + sb(6) &
              + sb(7)*tc0*tc0 + sb(8)*tc0 &
              + sb(9)*tc0*tc0*tc0 + sb(10)
         smo1(k) = a_ * smo2(k)**b_
!..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
               + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
               + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
               + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(1)*cse(1)*cse(1)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
              + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
              + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
         smoc(k) = a_ * smo2(k)**b_
!..Calculate bv_s+2 (th) moment.  Useful for riming.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(13) &
               + sa(4)*tc0*cse(13) + sa(5)*tc0*tc0 &
               + sa(6)*cse(13)*cse(13) + sa(7)*tc0*tc0*cse(13) &
               + sa(8)*tc0*cse(13)*cse(13) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(13)*cse(13)*cse(13)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(13) + sb(4)*tc0*cse(13) &
              + sb(5)*tc0*tc0 + sb(6)*cse(13)*cse(13) &
              + sb(7)*tc0*tc0*cse(13) + sb(8)*tc0*cse(13)*cse(13) &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(13)*cse(13)*cse(13)
         smoe(k) = a_ * smo2(k)**b_
!..Calculate 1+(bv_s+1)/2 (th) moment.  Useful for depositional growth.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(16) &
               + sa(4)*tc0*cse(16) + sa(5)*tc0*tc0 &
               + sa(6)*cse(16)*cse(16) + sa(7)*tc0*tc0*cse(16) &
               + sa(8)*tc0*cse(16)*cse(16) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(16)*cse(16)*cse(16)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(16) + sb(4)*tc0*cse(16) &
              + sb(5)*tc0*tc0 + sb(6)*cse(16)*cse(16) &
              + sb(7)*tc0*tc0*cse(16) + sb(8)*tc0*cse(16)*cse(16) &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(16)*cse(16)*cse(16)
         smof(k) = a_ * smo2(k)**b_
      enddo
!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for graupel.
!+---+-----------------------------------------------------------------+
      N0_min = gonv_max
      do k = kte, kts, -1
         if (temp(k).lt.270.65 .and. L_qr(k) .and. mvd_r(k).gt.100.E-6) then
            xslw1 = 4.01 + alog10(mvd_r(k))
         else
            xslw1 = 0.01
         endif
         ygra1 = 4.31 + alog10(max(5.E-5, rg(k)))
         zans1 = 3.1 + (100./(300.*xslw1*ygra1/(10./xslw1+1.+0.25*ygra1)+30.+10.*ygra1))
         N0_exp = 10.**(zans1)
         N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
         N0_min = MIN(N0_exp, N0_min)
         N0_exp = N0_min
         lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
         lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
         ilamg(k) = 1./lamg
         N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
!+---+-----------------------------------------------------------------+
!     if( debug_flag .and. k.lt.42) then
!        if (k.eq.41) write(mp_debug,*) 'DEBUG-GT:   K,   zans1,      rc,        rr,         rg,        N0_g'
!        if (k.eq.41) CALL wrf_debug(0, mp_debug)
!        write(mp_debug, 'a, i2, 1x, f6.3, 1x, 4(1x,e13.6,1x)')         &
!                   '  GT ', k, zans1, rc(k), rr(k), rg(k), N0_g(k)
!        CALL wrf_debug(0, mp_debug)
!     endif
!+---+-----------------------------------------------------------------+
      enddo
      endif
!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for rain.
!+---+-----------------------------------------------------------------+
      do k = kte, kts, -1
         lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
         ilamr(k) = 1./lamr
         mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
         N0_r(k) = nr(k)*org2*lamr**cre(2)
      enddo
!+---+-----------------------------------------------------------------+
!..Compute warm-rain process terms (except evap done later).
!+---+-----------------------------------------------------------------+
      do k = kts, kte
!..Rain self-collection follows Seifert, 1994 and drop break-up
!.. follows Verlinde and Cotton, 1993.                                        RAIN2M
!.. Reduced leading coeff because higher values overly force to singular
!.. mean size of rain.
         if (L_qr(k) .and. mvd_r(k).gt. D0r) then
!-GT      Ef_rr = 1.0
!-GT      if (mvd_r(k) .gt. 1500.0E-6) then
             Ef_rr = 2.0 - EXP(2300.0*(mvd_r(k)-1600.0E-6))
!-GT      endif
          pnr_rcr(k) = Ef_rr * 0.5*nr(k)*rr(k)
         endif
         mvd_c(k) = D0c
         if (.not. L_qc(k)) CYCLE
         xDc = MAX(D0c*1.E6, ((rc(k)/(am_r*Nt_c))**obmr) * 1.E6)
         lamc = (Nt_c*am_r* ccg(2) * ocg1 / rc(k))**obmr
         mvd_c(k) = (3.0+mu_c+0.672) / lamc
!..Autoconversion follows Berry & Reinhardt (1974) with characteristic
!.. diameters correctly computed from gamma distrib of cloud droplets.
         if (rc(k).gt. 0.01e-3) then
          Dc_g = ((ccg(3)*ocg2)**obmr / lamc) * 1.E6
          Dc_b = (xDc*xDc*xDc*Dc_g*Dc_g*Dc_g - xDc*xDc*xDc*xDc*xDc*xDc) &
                 **(1./6.)
          zeta1 = 0.5*((6.25E-6*xDc*Dc_b*Dc_b*Dc_b - 0.4) &
                     + abs(6.25E-6*xDc*Dc_b*Dc_b*Dc_b - 0.4))
          zeta = 0.027*rc(k)*zeta1
          taud = 0.5*((0.5*Dc_b - 7.5) + abs(0.5*Dc_b - 7.5)) + R1
          tau  = 3.72/(rc(k)*taud)
          prr_wau(k) = zeta/tau
          prr_wau(k) = MIN(DBLE(rc(k)*odts), prr_wau(k))
          pnr_wau(k) = prr_wau(k) / (am_r*mu_c*D0r*D0r*D0r)              ! RAIN2M
         endif
!..Rain collecting cloud water.  In CE, assume Dc<<Dr and vtc=~0.
         if (L_qr(k) .and. mvd_r(k).gt. D0r .and. mvd_c(k).gt. D0c) then
          lamr = 1./ilamr(k)
          idx = 1 + INT(nbr*DLOG(mvd_r(k)/Dr(1))/DLOG(Dr(nbr)/Dr(1)))
          idx = MIN(idx, nbr)
          Ef_rw = t_Efrw(idx, INT(mvd_c(k)*1.E6))
          prr_rcw(k) = rhof(k)*t1_qr_qc*Ef_rw*rc(k)*N0_r(k) &
                         *((lamr+fv_r)**(-cre(9)))
          prr_rcw(k) = MIN(DBLE(rc(k)*odts), prr_rcw(k))
         endif
      enddo
!+---+-----------------------------------------------------------------+
!..Compute all frozen hydrometeor species' process terms.
!+---+-----------------------------------------------------------------+
      if (.not. iiwarm) then
      do k = kts, kte
         vts_boost(k) = 1.5
!..Temperature lookup table indexes.
         tempc = temp(k) - 273.15
         idx_tc = MAX(1, MIN(NINT(-tempc), 45) )
         idx_t = INT( (tempc-2.5)/5. ) - 1
         idx_t = MAX(1, -idx_t)
         idx_t = MIN(idx_t, ntb_t)
         IT = MAX(1, MIN(NINT(-tempc), 31) )
!..Cloud water lookup table index.
         if (rc(k).gt. r_c(1)) then
          nic = NINT(ALOG10(rc(k)))
          do nn = nic-1, nic+1
             n = nn
             if ( (rc(k)/10.**nn).ge.1.0 .and. &
                  (rc(k)/10.**nn).lt.10.0) goto 141
          enddo
 141      continue
          idx_c = INT(rc(k)/10.**n) + 10*(n-nic2) - (n-nic2)
          idx_c = MAX(1, MIN(idx_c, ntb_c))
         else
          idx_c = 1
         endif
!..Cloud ice lookup table indexes.
         if (ri(k).gt. r_i(1)) then
          nii = NINT(ALOG10(ri(k)))
          do nn = nii-1, nii+1
             n = nn
             if ( (ri(k)/10.**nn).ge.1.0 .and. &
                  (ri(k)/10.**nn).lt.10.0) goto 142
          enddo
 142      continue
          idx_i = INT(ri(k)/10.**n) + 10*(n-nii2) - (n-nii2)
          idx_i = MAX(1, MIN(idx_i, ntb_i))
         else
          idx_i = 1
         endif
         if (ni(k).gt. Nt_i(1)) then
          nii = NINT(ALOG10(ni(k)))
          do nn = nii-1, nii+1
             n = nn
             if ( (ni(k)/10.**nn).ge.1.0 .and. &
                  (ni(k)/10.**nn).lt.10.0) goto 143
          enddo
 143      continue
          idx_i1 = INT(ni(k)/10.**n) + 10*(n-nii3) - (n-nii3)
          idx_i1 = MAX(1, MIN(idx_i1, ntb_i1))
         else
          idx_i1 = 1
         endif
!..Rain lookup table indexes.
         if (rr(k).gt. r_r(1)) then
          nir = NINT(ALOG10(rr(k)))
          do nn = nir-1, nir+1
             n = nn
             if ( (rr(k)/10.**nn).ge.1.0 .and. &
                  (rr(k)/10.**nn).lt.10.0) goto 144
          enddo
 144      continue
          idx_r = INT(rr(k)/10.**n) + 10*(n-nir2) - (n-nir2)
          idx_r = MAX(1, MIN(idx_r, ntb_r))
          lamr = 1./ilamr(k)
          lam_exp = lamr * (crg(3)*org2*org1)**bm_r
          N0_exp = org1*rr(k)/am_r * lam_exp**cre(1)
          nir = NINT(DLOG10(N0_exp))
          do nn = nir-1, nir+1
             n = nn
             if ( (N0_exp/10.**nn).ge.1.0 .and. &
                  (N0_exp/10.**nn).lt.10.0) goto 145
          enddo
 145      continue
          idx_r1 = INT(N0_exp/10.**n) + 10*(n-nir3) - (n-nir3)
          idx_r1 = MAX(1, MIN(idx_r1, ntb_r1))
         else
          idx_r = 1
          idx_r1 = ntb_r1
         endif
!..Snow lookup table index.
         if (rs(k).gt. r_s(1)) then
          nis = NINT(ALOG10(rs(k)))
          do nn = nis-1, nis+1
             n = nn
             if ( (rs(k)/10.**nn).ge.1.0 .and. &
                  (rs(k)/10.**nn).lt.10.0) goto 146
          enddo
 146      continue
          idx_s = INT(rs(k)/10.**n) + 10*(n-nis2) - (n-nis2)
          idx_s = MAX(1, MIN(idx_s, ntb_s))
         else
          idx_s = 1
         endif
!..Graupel lookup table index.
         if (rg(k).gt. r_g(1)) then
          nig = NINT(ALOG10(rg(k)))
          do nn = nig-1, nig+1
             n = nn
             if ( (rg(k)/10.**nn).ge.1.0 .and. &
                  (rg(k)/10.**nn).lt.10.0) goto 147
          enddo
 147      continue
          idx_g = INT(rg(k)/10.**n) + 10*(n-nig2) - (n-nig2)
          idx_g = MAX(1, MIN(idx_g, ntb_g))
          lamg = 1./ilamg(k)
          lam_exp = lamg * (cgg(3)*ogg2*ogg1)**bm_g
          N0_exp = ogg1*rg(k)/am_g * lam_exp**cge(1)
          nig = NINT(DLOG10(N0_exp))
          do nn = nig-1, nig+1
             n = nn
             if ( (N0_exp/10.**nn).ge.1.0 .and. &
                  (N0_exp/10.**nn).lt.10.0) goto 148
          enddo
 148      continue
          idx_g1 = INT(N0_exp/10.**n) + 10*(n-nig3) - (n-nig3)
          idx_g1 = MAX(1, MIN(idx_g1, ntb_g1))
         else
          idx_g = 1
          idx_g1 = ntb_g1
         endif
!..Deposition/sublimation prefactor (from Srivastava & Coen 1992).
         otemp = 1./temp(k)
         rvs = rho(k)*qvsi(k)
         rvs_p = rvs*otemp*(lsub*otemp*oRv - 1.)
         rvs_pp = rvs * ( otemp*(lsub*otemp*oRv - 1.) &
                         *otemp*(lsub*otemp*oRv - 1.) &
                         + (-2.*lsub*otemp*otemp*otemp*oRv) &
                         + otemp*otemp)
         gamsc = lsub*diffu(k)/tcond(k) * rvs_p
         alphsc = 0.5*(gamsc/(1.+gamsc))*(gamsc/(1.+gamsc)) &
                    * rvs_pp/rvs_p * rvs/rvs_p
         alphsc = MAX(1.E-9, alphsc)
         xsat = ssati(k)
         if (abs(xsat).lt. 1.E-9) xsat=0.
         t1_subl = 4.*PI*( 1.0 - alphsc*xsat &
                + 2.*alphsc*alphsc*xsat*xsat &
                - 5.*alphsc*alphsc*alphsc*xsat*xsat*xsat ) &
                / (1.+gamsc)
         if(ssati(k).gt.0.) then
           ! if supersaturated, figure out the effective saturation
           ! vapor density at the ice surface computed from Srivastava
           ! & Coen (1992).
           rv_ice_surface = rho(k)*qv(k) &
                - (t1_subl/4./pi)*ssati(k)*rvs
           Sice_tilde = (rho(k)*qv(k))/rv_ice_surface
           ! One can back out the surface temperature of the ice from
           !   Srivastava & Coen (1992, eqn. 5/12).  This temperature
           !   is the one relevant for deposition onto the ice
           !   surface, so that we use it (rather than the air
           !   temperature) to compute the equilibrium fractionation
           !   factor.  This is probably a small effect, but seems
           !   like a more faithful representation of the isotope
           !   physics, so that we use it here.
           tabs_ice_surface = temp(k) &
                + (t1_subl/4./pi)*ssati(k)*rvs*gamsc/rvs_p
           do nn = 1,niso
             alpha_equil_ice(nn) = &
                  alfaI_equilibrium(tabs_ice_surface,iso_index(nn))
           end do
           !bloss: commentary 1: interestingly, the surface temperature
           ! would be the same for snow, ice and graupel, assuming
           ! that the ventilation factors for heat and standard water
           ! have the same dependence on the particle Reynolds
           ! number.  (Is this really true?  It's assumed in
           ! Srivastava & Coen, 1992, below equation 7 if I
           ! understand properly.)  Perhaps, the way to think about 
           ! this is the extra vapor being deposited onto graupel due to
           ! ventilation is complemented by additional transfer of
           ! heat away from the graupel, so that its surface
           ! temperature will be more or less identical to that of a
           ! tiny ice particle with a negligible fall speed.
           ! As a result, the equilibrium fractionation factor is
           ! identical for the various types of ice.  The kinetic
           ! fractionation factors will differ due to the ratio of
           ! the standard ventilation factor to that of the heavy
           ! isotopologues.  
           !bloss: commentary 2: Note that the Thompson microphysics 
           ! neglects the change in surface temperature due to the
           ! accretion of liquid water when computing deposition,
           ! which would transfer additional heat to
           ! the ice particle and tend to reduce the deposition rate
           ! onto the ice.  Ideally, one would solve a version of
           ! equation 2 in Srivastava & Coen that included the latent
           ! heating of the ice/snow/graupel due to accretion of
           ! liquid water.
         else
           ! if subsaturated, no fractionation
           alpha_equil_ice(:) = 1.
          end if
!..Snow collecting cloud water.  In CE, assume Dc<<Ds and vtc=~0.
         if (L_qc(k) .and. mvd_c(k).gt. D0c) then
          xDs = 0.0
          if (L_qs(k)) xDs = smoc(k) / smob(k)
          if (xDs .gt. D0s) then
           idx = 1 + INT(nbs*DLOG(xDs/Ds(1))/DLOG(Ds(nbs)/Ds(1)))
           idx = MIN(idx, nbs)
           Ef_sw = t_Efsw(idx, INT(mvd_c(k)*1.E6))
           prs_scw(k) = rhof(k)*t1_qs_qc*Ef_sw*rc(k)*smoe(k)
          endif
!..Graupel collecting cloud water.  In CE, assume Dc<<Dg and vtc=~0.
          if (rg(k).ge. r_g(1) .and. mvd_c(k).gt. D0c) then
           xDg = (bm_g + mu_g + 1.) * ilamg(k)
           vtg = rhof(k)*av_g*cgg(6)*ogg3 * ilamg(k)**bv_g
           stoke_g = mvd_c(k)*mvd_c(k)*vtg*rho_w/(9.*visco(k)*xDg)
           if (xDg.gt. D0g) then
            if (stoke_g.ge.0.4 .and. stoke_g.le.10.) then
             Ef_gw = 0.55*ALOG10(2.51*stoke_g)
            elseif (stoke_g.lt.0.4) then
             Ef_gw = 0.0
            elseif (stoke_g.gt.10) then
             Ef_gw = 0.77
            endif
            prg_gcw(k) = rhof(k)*t1_qg_qc*Ef_gw*rc(k)*N0_g(k) &
                          *ilamg(k)**cge(9)
           endif
          endif
         endif
!..Rain collecting snow.  Cannot assume Wisner (1972) approximation
!.. or Mizuno (1990) approach so we solve the CE explicitly and store
!.. results in lookup table.
         if (rr(k).ge. r_r(1)) then
          if (rs(k).ge. r_s(1)) then
            !bloss: re-arrange rain-snow collision processes
            !  so that each process has a single source and
            !  single destination.  Requires addition of
            !  prg_scr, which is the conversion of rain to
            !  graupel following a collision with snow.
            ! snow collects rain --> snow (above or below 0C).
            prs_rcs(k) = tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
                         + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r)
            prs_rcs(k) = MIN(DBLE(rr(k)*odts), prs_rcs(k))
            ! rain collects snow --> rain (T>0C), graupel (T<0C)
            prr_rcs(k) = tcs_racs1(idx_s,idx_t,idx_r1,idx_r)           &
                        + tms_sacr1(idx_s,idx_t,idx_r1,idx_r)
            prr_rcs(k) = MIN(DBLE(rs(k)*odts), prr_rcs(k))
            ! loss of rain number due to prs_rcs
            pnr_rcs(k) = tnr_racs2(idx_s,idx_t,idx_r1,idx_r)            &   ! RAIN2M
                         + tnr_sacr2(idx_s,idx_t,idx_r1,idx_r)
           if (temp(k).lt.T_0) then
             ! prr_rcs --> prg_rcs (make graupel, not rain)
             !   source hydrometeor is snow.
             prg_rcs(k) = prr_rcs(k)
             prr_rcs(k) = 0.
             ! now account for rain conversion to graupel
             !  after collision with snow.  
             !  Source hydrometeor: rain.
             prg_scr(k) = tmr_racs1(idx_s,idx_t,idx_r1,idx_r) &
                           + tcr_sacr1(idx_s,idx_t,idx_r1,idx_r)
             prg_scr(k) = MIN(DBLE(rr(k)*odts), prg_scr(k))
             pnr_rcs(k) = pnr_rcs(k) &
                         + tnr_racs1(idx_s,idx_t,idx_r1,idx_r) &   ! RAIN2M
                         + tnr_sacr1(idx_s,idx_t,idx_r1,idx_r)
           endif
           pnr_rcs(k) = MIN(DBLE(nr(k)*odts), pnr_rcs(k))
          endif
!..Rain collecting graupel.  Cannot assume Wisner (1972) approximation
!.. or Mizuno (1990) approach so we solve the CE explicitly and store
!.. results in lookup table.
          if (rg(k).ge. r_g(1)) then
           if (temp(k).lt.T_0) then
            prg_rcg(k) = tmr_racg(idx_g1,idx_g,idx_r1,idx_r) &
                         + tcr_gacr(idx_g1,idx_g,idx_r1,idx_r)
            prg_rcg(k) = MIN(DBLE(rr(k)*odts), prg_rcg(k))
!bloss:  Make prg_rcg(k) a sink in rain tendency, so that each process
!           has a single source hydrometeor.
!bloss            prr_rcg(k) = -prg_rcg(k)
            pnr_rcg(k) = tnr_racg(idx_g1,idx_g,idx_r1,idx_r)            &   ! RAIN2M
                         + tnr_gacr(idx_g1,idx_g,idx_r1,idx_r)
            pnr_rcg(k) = MIN(DBLE(nr(k)*odts), pnr_rcg(k))
           else
            prr_rcg(k) = tcg_racg(idx_g1,idx_g,idx_r1,idx_r)
            prr_rcg(k) = MIN(DBLE(rg(k)*odts), prr_rcg(k))
!bloss:  Make prr_rcg(k) a sink in the graupel tendency, so that each process
!           has a single source hydrometeor.
!bloss            prg_rcg(k) = -prr_rcg(k)
           endif
          endif
         endif
!+---+-----------------------------------------------------------------+
!..Next IF block handles only those processes below 0C.
!+---+-----------------------------------------------------------------+
         if (temp(k).lt.T_0) then
          vts_boost(k) = 1.0
          rate_max = (qv(k)-qvsi(k))*rho(k)*odts*0.999
!..Freezing of water drops into graupel/cloud ice (Bigg 1953).
          if (rr(k).gt. r_r(1)) then
           prg_rfz(k) = tpg_qrfz(idx_r,idx_r1,idx_tc)*odts
           pri_rfz(k) = tpi_qrfz(idx_r,idx_r1,idx_tc)*odts
           pni_rfz(k) = tni_qrfz(idx_r,idx_r1,idx_tc)*odts
           pnr_rfz(k) = tnr_qrfz(idx_r,idx_r1,idx_tc)*odts                 ! RAIN2M
           pnr_rfz(k) = MIN(DBLE(nr(k)*odts), pnr_rfz(k))
          elseif (rr(k).gt. R1 .and. temp(k).lt.HGFR) then
           pri_rfz(k) = rr(k)*odts
           pnr_rfz(k) = nr(k)*odts                                         ! RAIN2M
           pni_rfz(k) = pnr_rfz(k)
          endif
          if (rc(k).gt. r_c(1)) then
           pri_wfz(k) = tpi_qcfz(idx_c,idx_tc)*odts
           pri_wfz(k) = MIN(DBLE(rc(k)*odts), pri_wfz(k))
           pni_wfz(k) = tni_qcfz(idx_c,idx_tc)*odts
           pni_wfz(k) = MIN(DBLE(Nt_c*odts), pri_wfz(k)/(2.*xm0i), &
                                pni_wfz(k))
          elseif (rc(k).gt. R1 .and. temp(k).lt.HGFR) then
           pri_wfz(k) = rc(k)*odts
           pni_wfz(k) = MIN(DBLE(Nt_c*odts), pri_wfz(k)/(2.*xm0i), &
                                pni_wfz(k))
          endif
!..Nucleate ice from deposition & condensation freezing (Cooper 1986)
!.. but only if water sat and T<-12C or 25%+ ice supersaturated.
          if ( (ssati(k).ge. 0.25) .or. (ssatw(k).gt. eps &
                                .and. temp(k).lt.261.15) ) then
           xnc = MIN(250.E3, TNO*EXP(ATO*(T_0-temp(k))))
           xni = ni(k) + (pni_rfz(k)+pni_wfz(k))*dtsave
           pni_inu(k) = 0.5*(xnc-xni + abs(xnc-xni))*odts
           pri_inu(k) = MIN(DBLE(rate_max), xm0i*pni_inu(k))
           pni_inu(k) = pri_inu(k)/xm0i
          endif
!..Deposition/sublimation of cloud ice (Srivastava & Coen 1992).
          if (L_qi(k)) then
           lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
           ilami = 1./lami
           xDi = MAX(DBLE(D0i), (bm_i + mu_i + 1.) * ilami)
           xmi = am_i*xDi**bm_i
           oxmi = 1./xmi
           pri_ide(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                  *oig1*cig(5)*ni(k)*ilami
           if (pri_ide(k) .lt. 0.0) then
            pri_ide(k) = MAX(DBLE(-ri(k)*odts), pri_ide(k), DBLE(rate_max))
            pni_ide(k) = pri_ide(k)*oxmi
            pni_ide(k) = MAX(DBLE(-ni(k)*odts), pni_ide(k))
           else
            pri_ide(k) = MIN(pri_ide(k), DBLE(rate_max))
            prs_ide(k) = (1.0D0-tpi_ide(idx_i,idx_i1))*pri_ide(k)
            pri_ide(k) = tpi_ide(idx_i,idx_i1)*pri_ide(k)
           endif
!..Some cloud ice needs to move into the snow category.  Use lookup
!.. table that resulted from explicit bin representation of distrib.
           if ( (idx_i.eq. ntb_i) .or. (xDi.gt. 5.0*D0s) ) then
            prs_iau(k) = ri(k)*.99*odts
            pni_iau(k) = ni(k)*.95*odts
           elseif (xDi.lt. 0.1*D0s) then
            prs_iau(k) = 0.
            pni_iau(k) = 0.
           else
            prs_iau(k) = tps_iaus(idx_i,idx_i1)*odts
            prs_iau(k) = MIN(DBLE(ri(k)*.99*odts), prs_iau(k))
            pni_iau(k) = tni_iaus(idx_i,idx_i1)*odts
            pni_iau(k) = MIN(DBLE(ni(k)*.95*odts), pni_iau(k))
           endif
          endif
!..Deposition/sublimation of snow/graupel follows Srivastava & Coen
!.. (1992).
          if (L_qs(k)) then
           C_snow = C_sqrd + (tempc+15.)*(C_cube-C_sqrd)/(-30.+15.)
           C_snow = MAX(C_sqrd, MIN(C_snow, C_cube))
           prs_sde(k) = C_snow*t1_subl*diffu(k)*ssati(k)*rvs &
                        * (t1_qs_sd*smo1(k) &
                         + t2_qs_sd*rhof2(k)*vsc2(k)*smof(k))
           if (prs_sde(k).lt. 0.) then
            prs_sde(k) = MAX(DBLE(-rs(k)*odts), prs_sde(k), DBLE(rate_max))
           else
            prs_sde(k) = MIN(prs_sde(k), DBLE(rate_max))
           endif
          endif
          if (L_qg(k) .and. ssati(k).lt. -eps) then
           prg_gde(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
               * N0_g(k) * (t1_qg_sd*ilamg(k)**cge(10) &
               + t2_qg_sd*vsc2(k)*rhof2(k)*ilamg(k)**cge(11))
           if (prg_gde(k).lt. 0.) then
            prg_gde(k) = MAX(DBLE(-rg(k)*odts), prg_gde(k), DBLE(rate_max))
           else
            prg_gde(k) = MIN(prg_gde(k), DBLE(rate_max))
           endif
          endif
!..Snow collecting cloud ice.  In CE, assume Di<<Ds and vti=~0.
          if (L_qi(k)) then
           lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
           ilami = 1./lami
           xDi = MAX(DBLE(D0i), (bm_i + mu_i + 1.) * ilami)
           xmi = am_i*xDi**bm_i
           oxmi = 1./xmi
           if (rs(k).ge. r_s(1)) then
            prs_sci(k) = t1_qs_qi*rhof(k)*Ef_si*ri(k)*smoe(k)
            pni_sci(k) = prs_sci(k) * oxmi
           endif
!..Rain collecting cloud ice.  In CE, assume Di<<Dr and vti=~0.
           if (rr(k).ge. r_r(1) .and. mvd_r(k).gt. 4.*xDi) then
            lamr = 1./ilamr(k)
            pri_rci(k) = rhof(k)*t1_qr_qi*Ef_ri*ri(k)*N0_r(k) &
                           *((lamr+fv_r)**(-cre(9)))
            pnr_rci(k) = rhof(k)*t1_qr_qi*Ef_ri*ni(k)*N0_r(k)           &   ! RAIN2M
                           *((lamr+fv_r)**(-cre(9)))
            pni_rci(k) = pri_rci(k) * oxmi
            prr_rci(k) = rhof(k)*t2_qr_qi*Ef_ri*ni(k)*N0_r(k) &
                           *((lamr+fv_r)**(-cre(8)))
            prr_rci(k) = MIN(DBLE(rr(k)*odts), prr_rci(k))
!bloss: use individual pri_rci and prr_rci terms in graupel tendency.
!         This simplifies both the limiting of process rates to
!         prevent over-depletion of rain and/or cloud ice and also
!         the implmentation of water isotopologues.  In both cases,
!         the individual processes are preferred because prg_rci
!         has two sources, so that it must be re-computed after the 
!         individual pri_rci and prr_rci processes are limited.
!         Better to eliminate prg_rci and use the individual rates
!         directly in the graupel tendency.
!bloss            prg_rci(k) = pri_rci(k) + prr_rci(k)
           endif
          endif
!..Ice multiplication from rime-splinters (Hallet & Mossop 1974).
          if (prg_gcw(k).gt. eps .and. tempc.gt.-8.0) then
            !bloss: Should we use some prediction of graupel surface
            ! temperature here, rather than air temperature???
            !  My impression from HM74 and from Pruppacher & Klett is
            ! that the graupel surface  temperature is the important
            ! thing.  One could compute this in a manner similar to
            ! that for rain surface temperature during evaporation,
            ! with latent heating (due to vapor deposition when
            ! prg_gde>0 and water freezing with various prg_***
            ! terms) balanced by convective cooling of the graupel
            ! particles. 
           tf = 0.
           if (tempc.ge.-5.0 .and. tempc.lt.-3.0) then
            tf = 0.5*(-3.0 - tempc)
           elseif (tempc.gt.-8.0 .and. tempc.lt.-5.0) then
            tf = 0.33333333*(8.0 + tempc)
           endif
           pni_ihm(k) = 3.5E8*tf*prg_gcw(k)
           pri_ihm(k) = xm0i*pni_ihm(k)
           prs_ihm(k) = prs_scw(k)/(prs_scw(k)+prg_gcw(k)) &
                          * pri_ihm(k)
           prg_ihm(k) = prg_gcw(k)/(prs_scw(k)+prg_gcw(k)) &
                          * pri_ihm(k)
          endif
!..A portion of rimed snow converts to graupel but some remains snow.
!.. Interp from 5 to 75% as riming factor increases from 5.0 to 30.0
!.. 0.028 came from (.75-.05)/(30.-5.).  This remains ad-hoc and should
!.. be revisited.
          if (prs_scw(k).gt.5.0*prs_sde(k) .and. &
                         prs_sde(k).gt.eps) then
           r_frac = MIN(30.0D0, prs_scw(k)/prs_sde(k))
           g_frac = MIN(0.75, 0.05 + (r_frac-5.)*.028)
           vts_boost(k) = MIN(1.5, 1.1 + (r_frac-5.)*.016)
           prg_scw(k) = g_frac*prs_scw(k)
           prs_scw(k) = (1. - g_frac)*prs_scw(k)
          endif
         else
!..Melt snow and graupel and enhance from collisions with liquid.
!.. We also need to sublimate snow and graupel if subsaturated.
          if (L_qs(k)) then
           prr_sml(k) = (tempc*tcond(k)-lvap0*diffu(k)*delQvs(k))       &
                      * (t1_qs_me*smo1(k) + t2_qs_me*rhof2(k)*vsc2(k)*smof(k))
           prr_sml(k) = prr_sml(k) + 4218.*olfus*tempc &
                                   * (prr_rcs(k)+prs_scw(k))
           prr_sml(k) = MIN(DBLE(rs(k)*odts), MAX(0.D0, prr_sml(k)))
           pnr_sml(k) = smo0(k)/rs(k)*prr_sml(k) * 10.0**(-0.75*tempc)      ! RAIN2M
           pnr_sml(k) = MIN(DBLE(smo0(k)*odts), pnr_sml(k))
           if (tempc.gt.3.5 .or. rs(k).lt.0.005E-3) pnr_sml(k)=0.0
           if (ssati(k).lt. 0.) then
            prs_sde(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                         * (t1_qs_sd*smo1(k) &
                          + t2_qs_sd*rhof2(k)*vsc2(k)*smof(k))
            prs_sde(k) = MAX(DBLE(-rs(k)*odts), prs_sde(k))
           endif
          endif
          if (L_qg(k)) then
           prr_gml(k) = (tempc*tcond(k)-lvap0*diffu(k)*delQvs(k))       &
                      * N0_g(k)*(t1_qg_me*ilamg(k)**cge(10)             &
                      + t2_qg_me*rhof2(k)*vsc2(k)*ilamg(k)**cge(11))
!-GT       prr_gml(k) = prr_gml(k) + 4218.*olfus*tempc &
!-GT                               * (prr_rcg(k)+prg_gcw(k))
           prr_gml(k) = MIN(DBLE(rg(k)*odts), MAX(0.D0, prr_gml(k)))
           pnr_gml(k) = N0_g(k)*cgg(2)*ilamg(k)**cge(2) / rg(k)         &   ! RAIN2M
                      * prr_gml(k) * 10.0**(-.25*tempc)
           if (tempc.gt.7.5 .or. rg(k).lt.0.005E-3) pnr_gml(k)=0.0
           if (ssati(k).lt. 0.) then
            prg_gde(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                * N0_g(k) * (t1_qg_sd*ilamg(k)**cge(10) &
                + t2_qg_sd*vsc2(k)*rhof2(k)*ilamg(k)**cge(11))
            prg_gde(k) = MAX(DBLE(-rg(k)*odts), prg_gde(k))
           endif
          endif
!.. This change will be required if users run adaptive time step that
!.. results in delta-t that is generally too long to allow cloud water
!.. collection by snow/graupel above melting temperature.
!.. Credit to Bjorn-Egil Nygaard for discovering.
          if (dt .gt. 120.) then
             prr_rcw(k)=prr_rcw(k)+prs_scw(k)+prg_gcw(k)
             prs_scw(k)=0.
             prg_gcw(k)=0.
          endif
         endif
          !bloss: handle fractionation of heavy isotopes due to
          ! deposition onto cloud ice, snow and graupel here.
          ! We assume that there is no fractionation due to
          ! sublimation, which neglects layering of isotopic
          ! composition onto ice, but could be not a bad
          ! approximation when the ice is heavily rimed. (??)
          ! Note: ventilation factor for cloud ice is one in Thompson
          ! scheme, so that the heavy isotope ventilation factor will
          ! also be one.
          if(pri_inu(k).gt.0.) then
            do nn = 1,niso
              ventilation_ratio = 1.
              alpha_kinetic = alfaK( alpha_equil_ice(nn), Sice_tilde, &
                   ventilation_ratio, iso_index(nn) )
              iso_pri_inu(k,nn) = &
                 alpha_equil_ice(nn)*alpha_kinetic*isor_v(k,nn)*pri_inu(k)
            end do
          else
            iso_pri_inu(k,:) = isor_i(k,:)*pri_inu(k) ! should never happen.
          end if
            
          if(pri_ide(k).gt.0.) then
            do nn = 1,niso
              ventilation_ratio = 1.
              alpha_kinetic = alfaK( alpha_equil_ice(nn), Sice_tilde, &
                   ventilation_ratio, iso_index(nn) )
              iso_pri_ide(k,nn) = &
                   alpha_equil_ice(nn)*alpha_kinetic*isor_v(k,nn)*pri_ide(k)
            end do
          else
            iso_pri_ide(k,:) = isor_i(k,:)*pri_ide(k) 
          end if
            
          if(prs_ide(k).gt.0.) then
            do nn = 1,niso
              ventilation_ratio = 1.
              alpha_kinetic = alfaK( alpha_equil_ice(nn), Sice_tilde, &
                   ventilation_ratio, iso_index(nn) )
              iso_prs_ide(k,nn) = &
                   alpha_equil_ice(nn)*alpha_kinetic*isor_v(k,nn)*prs_ide(k)
            end do
          else
            iso_prs_ide(k,:) = isor_i(k,:)*prs_ide(k) 
          end if
            
          if(prs_sde(k).gt.0.) then
          ! The ventilation factors for the light/standard isotope of
           ! snow are taken from the expression for prs_sde above.
            ventilation_factor_light = &
                 t1_qs_sd + t2_qs_sd*rhof2(k)*vsc2(k)*smof(k)/smo1(k)
            do nn = 1,niso
              ventilation_factor_heavy = t1_qs_sd &
                   + (iso_Sc3(nn)/Sc3)*t2_qs_sd &
                       *rhof2(k)*vsc2(k)*smof(k)/smo1(k)
              ventilation_ratio = MAX(1.,ventilation_factor_light) &
                   /MAX(1.,ventilation_factor_heavy)
              alpha_kinetic = alfaK( alpha_equil_ice(nn), Sice_tilde, &
                   ventilation_ratio, iso_index(nn) )
              iso_prs_sde(k,nn) = &
                   alpha_equil_ice(nn)*alpha_kinetic*isor_v(k,nn)*prs_sde(k)
            end do
          else
            ! sublimation, no fractionation.
            iso_prs_sde(k,:) = isor_s(k,:)*prs_sde(k) 
          end if
            
          if(prg_gde(k).gt.0.) then
            !bloss: Apparently, deposition onto graupel is prohibited in this
            !  microphysics, so that this code will likely be unused.  It seems
            !  useful to include the code, though, in case we (or someone else) 
            !  might want to enable deposition onto graupel in the future.
          ! The ventilation factors for the light/standard isotope of
           ! snow are taken from the expression for prs_gde above.
            !  Note that the ventilation factor is scaled by cgg(10)
            ! within t1_qg_sd and by (1/lambda)^(cge(10) in the
            ! expression below.  We divide by those to get the raw
            ! ventilation factor.
            ventilation_factor_light = &
                 (t1_qg_sd*ilamg(k)**cge(10) &
                 + t2_qg_sd*vsc2(k)*rhof2(k)*ilamg(k)**cge(11)) &
                 / (cgg(10)*ilamg(k)**cge(10))
            do nn = 1,niso
              ventilation_factor_heavy = &
                   (t1_qg_sd*ilamg(k)**cge(10) &
                   + (iso_Sc3(nn)/Sc3)*t2_qg_sd*vsc2(k)*rhof2(k)*ilamg(k)**cge(11)) &
                   / (cgg(10)*ilamg(k)**cge(10))
              ventilation_ratio = MAX(1.,ventilation_factor_light) &
                   /MAX(1.,ventilation_factor_heavy)
              alpha_kinetic = alfaK( alpha_equil_ice(nn), Sice_tilde, &
                   ventilation_ratio, iso_index(nn) )
              iso_prg_gde(k,nn) = &
                   alpha_equil_ice(nn)*alpha_kinetic*isor_v(k,nn)*prg_gde(k)
            end do
          else
            ! sublimation, no fractionation.
            iso_prg_gde(k,:) = isor_g(k,:)*prg_gde(k) 
          end if
      enddo
      endif
!+---+-----------------------------------------------------------------+
!..Ensure we do not deplete more hydrometeor species than exists.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
!..If ice supersaturated, ensure sum of depos growth terms does not
!.. deplete more vapor than possibly exists.  If subsaturated, limit
!.. sum of sublimation terms such that vapor does not reproduce ice
!.. supersat again.
         sump = pri_inu(k) + pri_ide(k) + prs_ide(k) &
              + prs_sde(k) + prg_gde(k)
         rate_max = (qv(k)-qvsi(k))*rho(k)*odts*0.999 !bloss: added *rho(k)
         if ( (sump.gt. eps .and. sump.gt. rate_max) .or. &
              (sump.lt. -eps .and. sump.lt. rate_max) ) then
          ratio = rate_max/sump
          pri_inu(k) = pri_inu(k) * ratio
          pri_ide(k) = pri_ide(k) * ratio
          pni_ide(k) = pni_ide(k) * ratio
          prs_ide(k) = prs_ide(k) * ratio
          prs_sde(k) = prs_sde(k) * ratio
          prg_gde(k) = prg_gde(k) * ratio
          !bloss: limit isotopic process rates as well,
          !  since these are directly proportional to those
          !  for the standard isotope.
          iso_pri_inu(k,:) = iso_pri_inu(k,:) * ratio
          iso_pri_ide(k,:) = iso_pri_ide(k,:) * ratio
          iso_prs_ide(k,:) = iso_prs_ide(k,:) * ratio
          iso_prs_sde(k,:) = iso_prs_sde(k,:) * ratio
          iso_prg_gde(k,:) = iso_prg_gde(k,:) * ratio
         endif
!..Cloud water conservation.
         sump = -prr_wau(k) - pri_wfz(k) - prr_rcw(k) &
                - prs_scw(k) - prg_scw(k) - prg_gcw(k)
         rate_max = -rc(k)*odts
         if (sump.lt. rate_max .and. L_qc(k)) then
          ratio = rate_max/sump
          prr_wau(k) = prr_wau(k) * ratio
          pri_wfz(k) = pri_wfz(k) * ratio
          prr_rcw(k) = prr_rcw(k) * ratio
          prs_scw(k) = prs_scw(k) * ratio
          prg_scw(k) = prg_scw(k) * ratio
          prg_gcw(k) = prg_gcw(k) * ratio
         endif
!..Cloud ice conservation.
         sump = pri_ide(k) - prs_iau(k) - prs_sci(k) &
                - pri_rci(k)
         rate_max = -ri(k)*odts
         if (sump.lt. rate_max .and. L_qi(k)) then
          ratio = rate_max/sump
          pri_ide(k) = pri_ide(k) * ratio
          prs_iau(k) = prs_iau(k) * ratio
          prs_sci(k) = prs_sci(k) * ratio
          pri_rci(k) = pri_rci(k) * ratio
          !bloss: limit isotopic process rates as well,
          !  since these are directly proportional to those
          !  for the standard isotope.
          iso_pri_ide(k,:) = iso_pri_ide(k,:) * ratio
         endif
!..Rain conservation.
         sump = -prg_rfz(k) - pri_rfz(k) - prr_rci(k) &
                - prs_rcs(k) - prg_scr(k) - prg_rcg(k) !bloss
         rate_max = -rr(k)*odts
         if (sump.lt. rate_max .and. L_qr(k)) then
          ratio = rate_max/sump
          prg_rfz(k) = prg_rfz(k) * ratio
          pri_rfz(k) = pri_rfz(k) * ratio
          prr_rci(k) = prr_rci(k) * ratio
          prs_rcs(k) = prs_rcs(k) * ratio !rain --> snow via collision
          prg_scr(k) = prg_scr(k) * ratio !rain --> graupel after collision with snow
          prg_rcg(k) = prg_rcg(k) * ratio !bloss: prg_rcg is sink for rain
         endif
!..Snow conservation.
         sump = prs_sde(k) - prs_ihm(k) - prr_sml(k) &
                - prr_rcs(k) - prg_rcs(k) !blos
         rate_max = -rs(k)*odts
         if (sump.lt. rate_max .and. L_qs(k)) then
          ratio = rate_max/sump
          prs_sde(k) = prs_sde(k) * ratio
          prs_ihm(k) = prs_ihm(k) * ratio
          prr_sml(k) = prr_sml(k) * ratio
          prr_rcs(k) = prr_rcs(k) * ratio ! snow --> rain via collision
          prg_rcs(k) = prg_rcs(k) * ratio ! snow --> graupel after collision w/rain
          !bloss: limit isotopic process rates as well,
          !  since these are directly proportional to those
          !  for the standard isotope.
          iso_prs_sde(k,:) = iso_prs_sde(k,:) * ratio
         endif
!..Graupel conservation.
         sump = prg_gde(k) - prg_ihm(k) - prr_gml(k) &
              - prr_rcg(k) !bloss: prg_rcg --> prr_rcg
         rate_max = -rg(k)*odts
         if (sump.lt. rate_max .and. L_qg(k)) then
          ratio = rate_max/sump
          prg_gde(k) = prg_gde(k) * ratio
          prg_ihm(k) = prg_ihm(k) * ratio
          prr_gml(k) = prr_gml(k) * ratio
          prr_rcg(k) = prr_rcg(k) * ratio !bloss: prr_rcg is sink for graupel
          !bloss: limit isotopic process rates as well,
          !  since these are directly proportional to those
          !  for the standard isotope.
          iso_prg_gde(k,:) = iso_prg_gde(k,:) * ratio
         endif
!..Re-enforce proper mass conservation for subsequent elements in case
!.. any of the above terms were altered.  Thanks P. Blossey. 2009Sep28
         pri_ihm(k) = prs_ihm(k) + prg_ihm(k)
!bloss: not needed after separating processes into prr_rcg (graup --> rain)
!          and prg_rcg (rain --> graupel)
!bloss         ratio = MIN( ABS(prr_rcg(k)), ABS(prg_rcg(k)) )
!bloss         prr_rcg(k) = ratio * SIGN(1.0, SNGL(prr_rcg(k)))
!bloss         prg_rcg(k) = -prr_rcg(k)
!!$         if (temp(k).gt.T_0) then
!!$            ratio = MIN( ABS(prr_rcs(k)), ABS(prs_rcs(k)) )
!!$            prr_rcs(k) = ratio * SIGN(1.0, SNGL(prr_rcs(k)))
!!$            prs_rcs(k) = -prr_rcs(k)
!!$         endif
      enddo
!+---+-----------------------------------------------------------------+
!..Calculate tendencies of all species but constrain the number of ice
!.. to reasonable values.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         orho = 1./rho(k)
         lfus2 = lsub - lvap(k)
!..Water vapor tendency
         qvten(k) = qvten(k) + (-pri_inu(k) - pri_ide(k) &
                      - prs_ide(k) - prs_sde(k) - prg_gde(k)) &
                      * orho
!..Cloud water tendency
         qcten(k) = qcten(k) + (-prr_wau(k) - pri_wfz(k) &
                      - prr_rcw(k) - prs_scw(k) - prg_scw(k) &
                      - prg_gcw(k)) &
                      * orho
!..Cloud ice mixing ratio tendency
         qiten(k) = qiten(k) + (pri_inu(k) + pri_ihm(k) &
                      + pri_wfz(k) + pri_rfz(k) + pri_ide(k) &
                      - prs_iau(k) - prs_sci(k) - pri_rci(k)) &
                      * orho
!..Cloud ice number tendency.
         niten(k) = niten(k) + (pni_inu(k) + pni_ihm(k) &
                      + pni_wfz(k) + pni_rfz(k) + pni_ide(k) &
                      - pni_iau(k) - pni_sci(k) - pni_rci(k)) &
                      * orho
!..Cloud ice mass/number balance; keep mass-wt mean size between
!.. 20 and 300 microns.  Also no more than 250 xtals per liter.
         xri=MAX(R1,(qi1d(k) + qiten(k)*dtsave)*rho(k))
         xni=MAX(R2,(ni1d(k) + niten(k)*dtsave)*rho(k))
         if (xri.gt. R1) then
           lami = (am_i*cig(2)*oig1*xni/xri)**obmi
           ilami = 1./lami
           xDi = (bm_i + mu_i + 1.) * ilami
           if (xDi.lt. 20.E-6) then
            lami = cie(2)/20.E-6
            xni = MIN(250.D3, cig(1)*oig2*xri/am_i*lami**bm_i)
            niten(k) = (xni-ni1d(k)*rho(k))*odts*orho
           elseif (xDi.gt. 300.E-6) then
            lami = cie(2)/300.E-6
            xni = cig(1)*oig2*xri/am_i*lami**bm_i
            niten(k) = (xni-ni1d(k)*rho(k))*odts*orho
           endif
         else
          niten(k) = -ni1d(k)*odts
         endif
         xni=MAX(0.,(ni1d(k) + niten(k)*dtsave)*rho(k))
         if (xni.gt.250.E3) &
                niten(k) = (250.E3-ni1d(k)*rho(k))*odts*orho
!..Rain tendency
         qrten(k) = qrten(k) + (prr_wau(k) + prr_rcw(k) &
                      + prr_sml(k) + prr_gml(k) &
                      + prr_rcs(k) - prs_rcs(k) - prg_scr(k) & !bloss: separate tendencies
                      + prr_rcg(k) - prg_rcg(k) & !bloss: separate tendencies
                      - prg_rfz(k) &
                      - pri_rfz(k) - prr_rci(k)) &
                      * orho
!..Rain number tendency
         nrten(k) = nrten(k) + (pnr_wau(k) + pnr_sml(k) + pnr_gml(k)    &
                      - (pnr_rfz(k) + pnr_rcr(k) + pnr_rcg(k)           &
                      + pnr_rcs(k) + pnr_rci(k)) )                      &
                      * orho
!..Rain mass/number balance; keep median volume diameter between
!.. 37 microns (D0r*0.75) and 2.5 mm.
         xrr=MAX(R1,(qr1d(k) + qrten(k)*dtsave)*rho(k))
         xnr=MAX(R2,(nr1d(k) + nrten(k)*dtsave)*rho(k))
         if (xrr.gt. R1) then
           lamr = (am_r*crg(3)*org2*xnr/xrr)**obmr
           mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
           if (mvd_r(k) .gt. 2.5E-3) then
              mvd_r(k) = 2.5E-3
              lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
              xnr = crg(2)*org3*xrr*lamr**bm_r / am_r
              nrten(k) = (xnr-nr1d(k)*rho(k))*odts*orho
           elseif (mvd_r(k) .lt. D0r*0.75) then
              mvd_r(k) = D0r*0.75
              lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
              xnr = crg(2)*org3*xrr*lamr**bm_r / am_r
              nrten(k) = (xnr-nr1d(k)*rho(k))*odts*orho
           endif
         else
           qrten(k) = -qr1d(k)*odts
           nrten(k) = -nr1d(k)*odts
         endif
!..Snow tendency
         qsten(k) = qsten(k) + (prs_iau(k) + prs_sde(k) &
                      + prs_sci(k) + prs_scw(k) &
                      + prs_rcs(k) - prr_rcs(k) - prg_rcs(k) & !bloss: separate tendencies
                      + prs_ide(k) - prs_ihm(k) - prr_sml(k)) &
                      * orho
!..Graupel tendency
         qgten(k) = qgten(k) + (prg_scw(k) + prg_rfz(k) &
                      + prg_gde(k) + prg_gcw(k) &
                      + prg_rcg(k) - prr_rcg(k) & !bloss: separate tendencies
                      + prr_rci(k) + pri_rci(k) & !bloss: separate tendencies
                      + prg_rcs(k) + prg_scr(k) - prg_ihm(k) &
                      - prr_gml(k)) &
                      * orho
!..Temperature tendency
         if (temp(k).lt.T_0) then
          tten(k) = tten(k) &
                    + ( lsub*ocp(k)*(pri_inu(k) + pri_ide(k) &
                                     + prs_ide(k) + prs_sde(k) &
                                     + prg_gde(k)) &
                     + lfus2*ocp(k)*(pri_wfz(k) + pri_rfz(k) &
                                     + prg_rfz(k) + prs_scw(k) &
                                     + prg_scw(k) + prg_gcw(k) &
                                     + prg_scr(k) + prs_rcs(k) & !bloss: prg_rcs --> prg_scr
                                     + prr_rci(k) + prg_rcg(k)) &
                       )*orho * (1-IFDRY)
         else
          tten(k) = tten(k) &
                    + ( lfus*ocp(k)*(-prr_sml(k) - prr_gml(k) &
                                     - prr_rcs(k) + prs_rcs(k) & !bloss: separate processes
                                     - prr_rcg(k)) &
                      + lsub*ocp(k)*(prs_sde(k) + prg_gde(k)) &
                       )*orho * (1-IFDRY)
         endif
      enddo
!+---+-----------------------------------------------------------------+
!..Define process rates for isotopic exchanges
!..Limit process rates to prevent over-depletion of isotopic mass
!..Define overall tendencies for isotopic mass
!+---+-----------------------------------------------------------------+
      do k = kts, kte
!=====================================================================
!-- Define process rates for isotopic mass here (except for special
!--         cases, which are handled above)
!=====================================================================
        ! sinks for cloud liquid -- NO FRACTIONATION
        iso_prr_wau(k,:) = prr_wau(k) * isor_c(k,:)
        iso_pri_wfz(k,:) = pri_wfz(k) * isor_c(k,:)
        iso_prr_rcw(k,:) = prr_rcw(k) * isor_c(k,:)
        iso_prs_scw(k,:) = prs_scw(k) * isor_c(k,:)
        iso_prg_scw(k,:) = prg_scw(k) * isor_c(k,:)
        iso_prg_gcw(k,:) = prg_gcw(k) * isor_c(k,:)
        ! sinks for cloud ice -- NO FRACTIONATION
        iso_prs_iau(k,:) = prs_iau(k) * isor_i(k,:)
        iso_prs_sci(k,:) = prs_sci(k) * isor_i(k,:)
        iso_pri_rci(k,:) = pri_rci(k) * isor_i(k,:)
        !NOTE: We are neglecting layering within ice species,
        !  so that the isotopic composition of vapor produced by sublimation
        !  is the mean composition of the ice species.
        iso_prg_rfz(k,:) = prg_rfz(k) * isor_r(k,:)
        iso_pri_rfz(k,:) = pri_rfz(k) * isor_r(k,:)
        iso_prr_rci(k,:) = prr_rci(k) * isor_r(k,:)
        iso_prr_rcs(k,:) = prr_rcs(k) * isor_s(k,:) ! now has single source: snow
        iso_prr_rcg(k,:) = prr_rcg(k) * isor_g(k,:) ! now has single source: graupel
        iso_prs_ihm(k,:) = prs_ihm(k) * isor_c(k,:)
        iso_prr_sml(k,:) = prr_sml(k) * isor_s(k,:)
        iso_prs_rcs(k,:) = prs_rcs(k) * isor_r(k,:) ! now has single source: rain
        iso_prg_ihm(k,:) = prg_ihm(k) * isor_c(k,:)
        iso_prr_gml(k,:) = prr_gml(k) * isor_g(k,:)
        iso_prg_rcg(k,:) = prg_rcg(k) * isor_r(k,:) ! now has single source: rain
        iso_prg_scr(k,:) = prg_scr(k) * isor_r(k,:) ! single source: rain
        iso_prg_rcs(k,:) = prg_rcs(k) * isor_s(k,:) ! now has single source: snow
        !========================================
        !-- THIS CODE DISABLES FRACTIONATION. ---
        !========================================
        if(.NOT.fractionate) then
          ! this code is useful for testing whether there is any
          ! spurious fractionation in the scheme.  Do not use in 
          ! any real simulations.
          iso_pri_inu(k,:) = pri_inu(k) * isor_v(k,:)
          iso_pri_ide(k,:) = pri_ide(k) * isor_v(k,:)
          iso_prs_ide(k,:) = prs_ide(k) * isor_v(k,:)
          iso_prs_sde(k,:) = prs_sde(k) * isor_v(k,:)
          iso_prg_gde(k,:) = prg_gde(k) * isor_v(k,:)
          if(pri_ide(k).lt.0.) iso_pri_ide(k,:) = pri_ide(k) * isor_i(k,:)
          if(prs_sde(k).lt.0.) iso_prs_sde(k,:) = prs_sde(k) * isor_s(k,:)
          if(prg_gde(k).lt.0.) iso_prg_gde(k,:) = prg_gde(k) * isor_g(k,:)
        end if
!=====================================================================
!-- LIMIT ISOTOPIC PROCESS RATES TO PREVENT OVER-DEPLETION OF ISOTOPES.
!=====================================================================
!bloss: For standard water, these processes are limited to prevent 
!.. them from over-depleting super- or subsaturation.  The heavy isotopic
!.. counterparts of these processes don't have the same constraint, so
!.. here we limit only to prevent the over-depletion of heavy vapor itself.
!.. Note that the standard process rates were limited before they were used
!.. to compute the isotopic 
!..Water vapor conservation.
        do nn = 1,niso
         sump = -iso_pri_inu(k,nn) - iso_pri_ide(k,nn) - iso_prs_ide(k,nn) &
              - iso_prs_sde(k,nn) - iso_prg_gde(k,nn)
         rate_max = -MAX(R1,rho(k)*iso_qv1d(k,nn))*odts
         if ( sump.lt. rate_max ) then
          ratio = rate_max/sump
          iso_pri_inu(k,nn) = iso_pri_inu(k,nn) * ratio
          iso_pri_ide(k,nn) = iso_pri_ide(k,nn) * ratio
          iso_prs_ide(k,nn) = iso_prs_ide(k,nn) * ratio
          iso_prs_sde(k,nn) = iso_prs_sde(k,nn) * ratio
          iso_prg_gde(k,nn) = iso_prg_gde(k,nn) * ratio
         endif
!..Cloud water conservation.
         sump = - iso_prr_wau(k,nn) - iso_pri_wfz(k,nn) - iso_prr_rcw(k,nn) &
                - iso_prs_scw(k,nn) - iso_prg_scw(k,nn) - iso_prg_gcw(k,nn)
         rate_max = -MAX(R1,rho(k)*iso_qc1d(k,nn))*odts
         if (sump.lt. rate_max .and. L_qc(k)) then
          ratio = rate_max/sump
          iso_prr_wau(k,nn) = iso_prr_wau(k,nn) * ratio
          iso_pri_wfz(k,nn) = iso_pri_wfz(k,nn) * ratio
          iso_prr_rcw(k,nn) = iso_prr_rcw(k,nn) * ratio
          iso_prs_scw(k,nn) = iso_prs_scw(k,nn) * ratio
          iso_prg_scw(k,nn) = iso_prg_scw(k,nn) * ratio
          iso_prg_gcw(k,nn) = iso_prg_gcw(k,nn) * ratio
         endif
!..Cloud ice conservation.
         sump = iso_pri_ide(k,nn) - iso_prs_iau(k,nn) - iso_prs_sci(k,nn) &
                - iso_pri_rci(k,nn)
         rate_max = -MAX(R1,rho(k)*iso_qi1d(k,nn))*odts
         if (sump.lt. rate_max .and. L_qi(k)) then
          ratio = rate_max/sump
          iso_pri_ide(k,nn) = iso_pri_ide(k,nn) * ratio
          iso_prs_iau(k,nn) = iso_prs_iau(k,nn) * ratio
          iso_prs_sci(k,nn) = iso_prs_sci(k,nn) * ratio
          iso_pri_rci(k,nn) = iso_pri_rci(k,nn) * ratio
         endif
!..Rain conservation.
         sump = - iso_prg_rfz(k,nn) - iso_pri_rfz(k,nn) - iso_prr_rci(k,nn) &
                - iso_prs_rcs(k,nn) - iso_prg_scr(k,nn) - iso_prg_rcg(k,nn)
         rate_max = -MAX(R1,rho(k)*iso_qr1d(k,nn))*odts
         if (sump.lt. rate_max .and. L_qr(k)) then
          ratio = rate_max/sump
          iso_prg_rfz(k,nn) = iso_prg_rfz(k,nn) * ratio
          iso_pri_rfz(k,nn) = iso_pri_rfz(k,nn) * ratio
          iso_prr_rci(k,nn) = iso_prr_rci(k,nn) * ratio
          iso_prs_rcs(k,nn) = iso_prs_rcs(k,nn) * ratio !prs_rcs is sink for rain
          iso_prg_scr(k,nn) = iso_prg_scr(k,nn) * ratio !prg_scr is sink for rain
          iso_prg_rcg(k,nn) = iso_prg_rcg(k,nn) * ratio !prg_rcg is sink for rain
         endif
!..Snow conservation.
         sump = iso_prs_sde(k,nn) - iso_prs_ihm(k,nn) - iso_prr_sml(k,nn) &
                - iso_prr_rcs(k,nn) - iso_prg_rcs(k,nn)
         rate_max = -MAX(R1,rho(k)*iso_qs1d(k,nn))*odts
         if (sump.lt. rate_max .and. L_qs(k)) then
          ratio = rate_max/sump
          iso_prs_sde(k,nn) = iso_prs_sde(k,nn) * ratio
          iso_prs_ihm(k,nn) = iso_prs_ihm(k,nn) * ratio
          iso_prr_sml(k,nn) = iso_prr_sml(k,nn) * ratio
          iso_prr_rcs(k,nn) = iso_prr_rcs(k,nn) * ratio !snow --> rain
          iso_prg_rcs(k,nn) = iso_prg_rcs(k,nn) * ratio !snow --> graupel
         endif
!..Graupel conservation.
         sump = iso_prg_gde(k,nn) - iso_prg_ihm(k,nn) - iso_prr_gml(k,nn) &
              - iso_prr_rcg(k,nn)
         rate_max = -MAX(R1,rho(k)*iso_qg1d(k,nn))*odts
         if (sump.lt. rate_max .and. L_qg(k)) then
          ratio = rate_max/sump
          iso_prg_gde(k,nn) = iso_prg_gde(k,nn) * ratio
          iso_prg_ihm(k,nn) = iso_prg_ihm(k,nn) * ratio
          iso_prr_gml(k,nn) = iso_prr_gml(k,nn) * ratio
          iso_prr_rcg(k,nn) = iso_prr_rcg(k,nn) * ratio !prr_rcg is sink for graup
         endif
!..Re-enforce proper mass conservation for subsequent elements in case
!.. any of the above terms were altered.  Thanks P. Blossey. 2009Sep28
         iso_pri_ihm(k,nn) = iso_prs_ihm(k,nn) + iso_prg_ihm(k,nn)
       end do ! do nn = 1,niso
!!$         ratio = MIN( ABS(prr_rcg(k)), ABS(prg_rcg(k)) )
!!$         prr_rcg(k) = ratio * SIGN(1.0, SNGL(prr_rcg(k)))
!!$         prg_rcg(k) = -prr_rcg(k)
!!$         if (temp(k).gt.T_0) then
!!$            ratio = MIN( ABS(prr_rcs(k)), ABS(prs_rcs(k)) )
!!$            prr_rcs(k) = ratio * SIGN(1.0, SNGL(prr_rcs(k)))
!!$            prs_rcs(k) = -prr_rcs(k)
!!$         endif
     enddo ! do k = kts,kte
!=====================================================================
!..Calculate tendencies for masses of water isotopologues
!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
     do nn = 1,niso
      do k = kts, kte
         orho = 1./rho(k)
!..Water vapor tendency
         iso_qvten(k,nn) = iso_qvten(k,nn) + (-iso_pri_inu(k,nn) - iso_pri_ide(k,nn) &
                      - iso_prs_ide(k,nn) - iso_prs_sde(k,nn) - iso_prg_gde(k,nn)) &
                      * orho
!..Cloud water tendency
         iso_qcten(k,nn) = iso_qcten(k,nn) + (-iso_prr_wau(k,nn) - iso_pri_wfz(k,nn) &
                      - iso_prr_rcw(k,nn) - iso_prs_scw(k,nn) - iso_prg_scw(k,nn) &
                      - iso_prg_gcw(k,nn)) &
                      * orho
!..Cloud ice mixing ratio tendency
         iso_qiten(k,nn) = iso_qiten(k,nn) + (iso_pri_inu(k,nn) + iso_pri_ihm(k,nn) &
                      + iso_pri_wfz(k,nn) + iso_pri_rfz(k,nn) + iso_pri_ide(k,nn) &
                      - iso_prs_iau(k,nn) - iso_prs_sci(k,nn) - iso_pri_rci(k,nn)) &
                      * orho
!..Rain tendency
         iso_qrten(k,nn) = iso_qrten(k,nn) + (iso_prr_wau(k,nn) + iso_prr_rcw(k,nn) &
                      + iso_prr_sml(k,nn) + iso_prr_gml(k,nn) &
                      + iso_prr_rcs(k,nn) - iso_prs_rcs(k,nn) - iso_prg_scr(k,nn) &
                      + iso_prr_rcg(k,nn) - iso_prg_rcg(k,nn) &
                      - iso_prg_rfz(k,nn) &
                      - iso_pri_rfz(k,nn) - iso_prr_rci(k,nn)) &
                      * orho
!..Snow tendency
         iso_qsten(k,nn) = iso_qsten(k,nn) + (iso_prs_iau(k,nn) + iso_prs_sde(k,nn) &
                      + iso_prs_sci(k,nn) + iso_prs_scw(k,nn) &
                      + iso_prs_rcs(k,nn) - iso_prr_rcs(k,nn) - iso_prg_rcs(k,nn) &
                      + iso_prs_ide(k,nn) - iso_prs_ihm(k,nn) - iso_prr_sml(k,nn)) &
                      * orho
!..Graupel tendency
         iso_qgten(k,nn) = iso_qgten(k,nn) + (iso_prg_scw(k,nn) + iso_prg_rfz(k,nn) &
                      + iso_prg_gde(k,nn) + iso_prg_gcw(k,nn) &
                      + iso_prg_rcg(k,nn) - iso_prr_rcg(k,nn) &
                      + iso_prr_rci(k,nn) + iso_pri_rci(k,nn) & !bloss: separate tendencies
                      + iso_prg_rcs(k,nn) + iso_prg_scr(k,nn) &
                      - iso_prg_ihm(k,nn) - iso_prr_gml(k,nn)) &
                      * orho
      enddo ! do k = kts,kte
    end do ! do nn = 1,niso
!+---+-----------------------------------------------------------------+
!..Update variables for TAU+1 before condensation & sedimention.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k) + DT*tten(k)
         otemp = 1./temp(k)
         tempc = temp(k) - 273.15
         qv(k) = MAX(1.E-10, qv1d(k) + DT*qvten(k))
         do nn = 1,niso
           iso_qv(k,nn) = MAX(1.E-10*isor_v(k,nn), &
                iso_qv1d(k,nn) + DT*iso_qvten(k,nn) )
         end do
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         rhof(k) = SQRT(RHO_NOT/rho(k))
         rhof2(k) = SQRT(rhof(k))
         qvs(k) = rslf(pres(k), temp(k))
         ssatw(k) = qv(k)/qvs(k) - 1.
         if (abs(ssatw(k)).lt. eps) ssatw(k) = 0.0
         diffu(k) = 2.11E-5*(temp(k)/273.15)**1.94 * (101325./pres(k))
         if (tempc .ge. 0.0) then
            visco(k) = (1.718+0.0049*tempc)*1.0E-5
         else
            visco(k) = (1.718+0.0049*tempc-1.2E-5*tempc*tempc)*1.0E-5
         endif
         vsc2(k) = SQRT(rho(k)/visco(k))
         lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
         tcond(k) = (5.69 + 0.0168*tempc)*1.0E-5 * 418.936
         ocp(k) = 1./(Cp*(1.+0.887*qv(k)))
         lvt2(k)=lvap(k)*lvap(k)*ocp(k)*oRv*otemp*otemp
         if ((qc1d(k) + qcten(k)*DT) .gt. R1) then
            rc(k) = (qc1d(k) + qcten(k)*DT)*rho(k)
            iso_rc(k,:) = (iso_qc1d(k,:) + iso_qcten(k,:)*DT)*rho(k)
            L_qc(k) = .true.
         else
            rc(k) = R1
            L_qc(k) = .false.
            iso_rc(k,:) = R1
            ! push any extra isotopic cloud liquid mass into vapor
            iso_qvten(k,:) = iso_qvten(k,:) &
                 + MAX(0., iso_qc1d(k,:)*odt + iso_qcten(k,:) )
            iso_qcten(k,:) = -iso_qc1d(k,:)*odt
         endif
         if ((qi1d(k) + qiten(k)*DT) .gt. R1) then
            ri(k) = (qi1d(k) + qiten(k)*DT)*rho(k)
            iso_ri(k,:) = (iso_qi1d(k,:) + iso_qiten(k,:)*DT)*rho(k)
            ni(k) = MAX(R2, (ni1d(k) + niten(k)*DT)*rho(k))
            L_qi(k) = .true. 
         else
            ri(k) = R1
            ni(k) = R2
            L_qi(k) = .false.
            iso_ri(k,:) = R1
            ! push any extra isotopic cloud ice mass into vapor
            iso_qvten(k,:) = iso_qvten(k,:) &
                 + MAX(0., iso_qi1d(k,:)*odt + iso_qiten(k,:) )
            iso_qiten(k,:) = -iso_qi1d(k,:)*odt
         endif
         if ((qr1d(k) + qrten(k)*DT) .gt. R1) then
            rr(k) = (qr1d(k) + qrten(k)*DT)*rho(k)
            iso_rr(k,:) = (iso_qr1d(k,:) + iso_qrten(k,:)*DT)*rho(k)
            nr(k) = MAX(R2, (nr1d(k) + nrten(k)*DT)*rho(k))
            L_qr(k) = .true.
            lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
            mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
            if (mvd_r(k) .gt. 2.5E-3) then
               mvd_r(k) = 2.5E-3
               lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
               nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
            elseif (mvd_r(k) .lt. D0r*0.75) then
               mvd_r(k) = D0r*0.75
               lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
               nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
            endif
         else
            rr(k) = R1
            nr(k) = R2
            L_qr(k) = .false.
            iso_rr(k,:) = R1
            ! push any extra isotopic rain mass into vapor
            iso_qvten(k,:) = iso_qvten(k,:) &
                 + MAX(0., iso_qr1d(k,:)*odt + iso_qrten(k,:) )
            iso_qrten(k,:) = -iso_qr1d(k,:)*odt
         endif
               
         if ((qs1d(k) + qsten(k)*DT) .gt. R1) then
            rs(k) = (qs1d(k) + qsten(k)*DT)*rho(k)
            iso_rs(k,:) = (iso_qs1d(k,:) + iso_qsten(k,:)*DT)*rho(k)
            L_qs(k) = .true.
         else
            rs(k) = R1
            L_qs(k) = .false.
            iso_rs(k,:) = R1
            ! push any extra isotopic snow mass into vapor
            iso_qvten(k,:) = iso_qvten(k,:) &
                 + MAX(0., iso_qs1d(k,:)*odt + iso_qsten(k,:) )
            iso_qsten(k,:) = -iso_qs1d(k,:)*odt
         endif
         if ((qg1d(k) + qgten(k)*DT) .gt. R1) then
            rg(k) = (qg1d(k) + qgten(k)*DT)*rho(k)
            iso_rg(k,:) = (iso_qg1d(k,:) + iso_qgten(k,:)*DT)*rho(k)
            L_qg(k) = .true.
         else
            rg(k) = R1
            L_qg(k) = .false.
            iso_rg(k,:) = R1
            ! push any extra isotopic graupel mass into vapor
            iso_qvten(k,:) = iso_qvten(k,:) &
                 + MAX(0., iso_qg1d(k,:)*odt + iso_qgten(k,:) )
            iso_qgten(k,:) = -iso_qg1d(k,:)*odt
         endif
      enddo
      ! initialize to unreasonable value.  This should make obvious
      ! any spurious use of this ratio (e.g., multiplying the cloud
      ! liquid isotope ratio by a process rate when no cloud liquid exists)
      isor_v(:,:) = 1.e8
      isor_c(:,:) = 1.e10
      isor_i(:,:) = 1.e12
      isor_r(:,:) = 1.e14
      isor_s(:,:) = 1.e16
      isor_g(:,:) = 1.e18
      !Re-define isotope ratios based on updated mixing ratios
      do nn = 1,niso
        do k = kts, kte
          isor_v(k,nn) = iso_qv(k,nn)/qv(k)
          if(L_qc(k)) isor_c(k,nn) = iso_rc(k,nn)/rc(k)
          if(L_qi(k)) isor_i(k,nn) = iso_ri(k,nn)/ri(k)
          if(L_qr(k)) isor_r(k,nn) = iso_rr(k,nn)/rr(k)
          if(L_qs(k)) isor_s(k,nn) = iso_rs(k,nn)/rs(k)
          if(L_qg(k)) isor_g(k,nn) = iso_rg(k,nn)/rg(k)
        end do
      end do
!+---+-----------------------------------------------------------------+
!..With tendency-updated mixing ratios, recalculate snow moments and
!.. intercepts/slopes of graupel and rain.
!+---+-----------------------------------------------------------------+
      if (.not. iiwarm) then
      do k = kts, kte
         if (.not. L_qs(k)) CYCLE
         tc0 = MIN(-0.1, temp(k)-273.15)
         smob(k) = rs(k)*oams
!..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
!.. then we must compute actual 2nd moment and use as reference.
         if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
            smo2(k) = smob(k)
         else
            loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
               + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
               + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
               + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
               + sa(10)*bm_s*bm_s*bm_s
            a_ = 10.0**loga_
            b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
               + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
               + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
               + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
               + sb(10)*bm_s*bm_s*bm_s
            smo2(k) = (smob(k)/a_)**(1./b_)
         endif
!..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
               + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
               + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
               + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(1)*cse(1)*cse(1)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
              + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
              + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
         smoc(k) = a_ * smo2(k)**b_
!..Calculate bm_s+bv_s (th) moment.  Useful for sedimentation.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(14) &
               + sa(4)*tc0*cse(14) + sa(5)*tc0*tc0 &
               + sa(6)*cse(14)*cse(14) + sa(7)*tc0*tc0*cse(14) &
               + sa(8)*tc0*cse(14)*cse(14) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(14)*cse(14)*cse(14)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(14) + sb(4)*tc0*cse(14) &
              + sb(5)*tc0*tc0 + sb(6)*cse(14)*cse(14) &
              + sb(7)*tc0*tc0*cse(14) + sb(8)*tc0*cse(14)*cse(14) &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(14)*cse(14)*cse(14)
         smod(k) = a_ * smo2(k)**b_
      enddo
!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for graupel.
!+---+-----------------------------------------------------------------+
      N0_min = gonv_max
      do k = kte, kts, -1
         if (temp(k).lt.270.65 .and. L_qr(k) .and. mvd_r(k).gt.100.E-6) then
            xslw1 = 4.01 + alog10(mvd_r(k))
         else
            xslw1 = 0.01
         endif
         ygra1 = 4.31 + alog10(max(5.E-5, rg(k)))
         zans1 = 3.1 + (100./(300.*xslw1*ygra1/(10./xslw1+1.+0.25*ygra1)+30.+10.*ygra1))
         N0_exp = 10.**(zans1)
         N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
         N0_min = MIN(N0_exp, N0_min)
         N0_exp = N0_min
         lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
         lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
         ilamg(k) = 1./lamg
         N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
      enddo
      endif
!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for rain.
!+---+-----------------------------------------------------------------+
      do k = kte, kts, -1
         lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
         ilamr(k) = 1./lamr
         mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
         N0_r(k) = nr(k)*org2*lamr**cre(2)
      enddo
!+---+-----------------------------------------------------------------+
!..Cloud water condensation and evaporation.  Newly formulated using
!.. Newton-Raphson iterations (3 should suffice) as provided by B. Hall.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         if ( (ssatw(k).gt. eps) .or. (ssatw(k).lt. -eps .and. &
                   L_qc(k)) ) then
          clap = (qv(k)-qvs(k))/(1. + lvt2(k)*qvs(k))
          do n = 1, 3
             fcd = qvs(k)* EXP(lvt2(k)*clap) - qv(k) + clap
             dfcd = qvs(k)*lvt2(k)* EXP(lvt2(k)*clap) + 1.
             clap = clap - fcd/dfcd
          enddo
          !bloss: rho(k)*clap has the same units as rc(k).
          !  Previously, this has rc(k) + clap, which was inconsistent.
          xrc = rc(k) + rho(k)*clap
          if (xrc.gt. 0.0) then
             prw_vcd(k) = clap*odt
             L_qc(k) = .true.
          else
             prw_vcd(k) = -rc(k)/rho(k)*odts
             L_qc(k) = .false.
          endif
          qcten(k) = qcten(k) + prw_vcd(k)
          qvten(k) = qvten(k) - prw_vcd(k)
          tten(k) = tten(k) + lvap(k)*ocp(k)*prw_vcd(k)*(1-IFDRY)
          rc(k) = MAX(R1, (qc1d(k) + DT*qcten(k))*rho(k))
          qv(k) = MAX(1.E-10, qv1d(k) + DT*qvten(k))
          temp(k) = t1d(k) + DT*tten(k)
          rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
          qvs(k) = rslf(pres(k), temp(k))
          ssatw(k) = qv(k)/qvs(k) - 1.
          ! since this microphysics scheme assumes saturation 
          !  adjustment, allow the isotopes to equilibrate as well.
          !  This assumption prevents kinetic fractionation of vapor
          ! deposition onto or evaporation from small liquid droplets
          ! in super/subsaturated conditions.  However, I don't know
          ! how to allow kinetic fractionation in a way that's
          ! consistent with a saturation adjustment scheme that works
          ! instantaneously to remove any supersaturation with
          ! respect to liquid or to remove any cloud liquid in
          ! subsaturated conditions.
          ! That being said, perform this equilibration after the rain
          !   evaporation/equilibration process.  In test with rain
          !   and cloud equilibrating in a saturated environment, the
          !   placement of the equilibration after the evaporation
          !   improved the equilibration of cloud liquid.  Not sure
          !   entirely why this happened, but I think that it's useful
          !   to consider the vapor equation as the crucial one.  Vapor
          !   is exchanging with both cloud and rain.  The timescale of
          !   exchange with cloud is fast, so that you can assume that
          !   the isotope ratio of vapor is constantly in equilibrium
          !   with cloud during the adjustment process.  As a result,
          !   it makes sense to equilibrate the cloud and vapor after
          !   the rain exchange process.  The difference with the
          !   standard isotope is that rain won't exchange with vapor
          !   when vapor and cloud are in equilibrium (as they are
          !   after saturation adjustment), so that the standard water
          !   microphysics doesn't have to worry about this subtlety
          !   and can perform saturation adjustment with cloud and rain
          !   evaporation sequentially, since the processes don't
          !   interact in the same way they do for heavy isotopologues.
         endif
      enddo
!+---+-----------------------------------------------------------------+
!.. If still subsaturated, allow rain to evaporate, following
!.. Srivastava & Coen (1992).
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         if ( L_qr(k) ) then
           ! There is exchange of water isotopologues even if there
           !   is no evaporation of the standard isotope, so that
           !   we need to compute the isotopic tendency whenever
           !   rain is present.
          tempc = temp(k) - 273.15
          otemp = 1./temp(k)
          rhof(k) = SQRT(RHO_NOT/rho(k))
          rhof2(k) = SQRT(rhof(k))
          diffu(k) = 2.11E-5*(temp(k)/273.15)**1.94 * (101325./pres(k))
          if (tempc .ge. 0.0) then
             visco(k) = (1.718+0.0049*tempc)*1.0E-5
          else
             visco(k) = (1.718+0.0049*tempc-1.2E-5*tempc*tempc)*1.0E-5
          endif
          vsc2(k) = SQRT(rho(k)/visco(k))
          lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
          tcond(k) = (5.69 + 0.0168*tempc)*1.0E-5 * 418.936
          ocp(k) = 1./(Cp*(1.+0.887*qv(k)))
          rvs = rho(k)*qvs(k)
          rvs_p = rvs*otemp*(lvap(k)*otemp*oRv - 1.)
          rvs_pp = rvs * ( otemp*(lvap(k)*otemp*oRv - 1.) &
                          *otemp*(lvap(k)*otemp*oRv - 1.) &
                          + (-2.*lvap(k)*otemp*otemp*otemp*oRv) &
                          + otemp*otemp)
          gamsc = lvap(k)*diffu(k)/tcond(k) * rvs_p
          alphsc = 0.5*(gamsc/(1.+gamsc))*(gamsc/(1.+gamsc)) &
                     * rvs_pp/rvs_p * rvs/rvs_p
          alphsc = MAX(1.E-9, alphsc)
          xsat   = MIN(-1.E-9, ssatw(k))
          t1_evap = 2.*PI*( 1.0 - alphsc*xsat  &
                 + 2.*alphsc*alphsc*xsat*xsat  &
                 - 5.*alphsc*alphsc*alphsc*xsat*xsat*xsat ) &
                 / (1.+gamsc)
          lamr = 1./ilamr(k)
!..Rapidly eliminate near zero values when low humidity (<95%)
          if (qv(k)/qvs(k) .lt. 0.95 .AND. rr(k)/rho(k).le.1.E-8) then
          prv_rev(k) = rr(k)/rho(k)*odts
          !bloss: Since all rain is removed, remove all isotopic
          !  mass as well.
          iso_prv_rev(k,:) = iso_rr(k,:)/rho(k)*odts
          else
            eva_factor = 1. !bloss: use this here, and then again below.
!bloss          if ( (ssatw(k).lt. -eps) .and. (.not.(prw_vcd(k).gt. 0.)) ) then
            ! Only evaporate the standard isotope if subsaturated and
            !   there hasn't been any deposition onto cloud water.  This
            !   makes the evaporation treatment of the standard
            !   isotope consistent between the isotopic and non
            !   -isotopic versions of the code.
          prv_rev(k) = t1_evap*diffu(k)*(-ssatw(k))*N0_r(k)*rvs &
              * (t1_qr_ev*ilamr(k)**cre(10) &
              + t2_qr_ev*vsc2(k)*rhof2(k)*((lamr+0.5*fv_r)**(-cre(11))))
!bloss: Since we are allowing prv_rev to go in both directions in the
! isotopic code for consistency with the heavy isotope implementation,
! Change limiting so that we don't chop off deposition onto rain.
! Though deposition onto rain will likely be small since the cloud
! removed most all of the saturation excess above.
          if(qvs(k).gt.qv(k)) then
            rate_max = MIN((rr(k)/rho(k)*odts), ABS(qvs(k)-qv(k))*odts)
          else
            rate_max = rr(k)/rho(k)*odts
          end if
          ! compute eva_factor
          eva_factor = MIN(1.,DBLE(rate_max)*rho(k)/MAX(R1,prv_rev(k)))
          ! bloss:  NOTE that the factor of rho is being removed
          !  from prv_rev at this point!!! 
          prv_rev(k) = MIN(DBLE(rate_max), prv_rev(k)/rho(k))
!bloss         end if ! subsaturated and no deposition onto cloud liquid.
          if(fractionate) then
            ! wait until after the rain is updated to compute the
            ! isotopic exchange between rain and vapor.
            
            ! First, compute the effective saturation vapor density at the 
            !   surface of the rain, as computed from Srivastava & Coen (1992).
            rv_rain_surface = rho(k)*qv(k) - (t1_evap/2./pi)*ssatw(k)*rvs 
            ! bloss (random question): Why 2pi for rain and 4pi for ice??
            ! This is the effective supersaturation
            Sliq_tilde = (rho(k)*qv(k))/rv_rain_surface
            ! Next, back out the surface temperature of the rain from
            !   Srivastava & Coen (1992, eqn. 5/12).  This temperature
            !   is the one that is important for the equilibrium exchange
            !   of water isotopologues under the saturated conditions
            ! at the surface of the raindrop, so that we use it
            ! (rather than the air temperature) to compute the
            ! equilibrium fractionation factor.  This is probably a
            ! small effect, but seems like a more faithful
            ! representation of the isotope physics, so that we use
            ! it here. 
            tabs_rain_surface = temp(k) &
                 + (t1_evap/2./pi)*ssatw(k)*rvs*gamsc/rvs_p
            ! for fun (and diagnostic purposes), compute the
            ! ventilation factor for standard rain.
            ventilation_factor_light = (t1_qr_ev*ilamr(k)**cre(10) &
              + t2_qr_ev*vsc2(k)*rhof2(k)*((lamr+0.5*fv_r)**(-cre(11))))&
              /(crg(10)*ilamr(k)**cre(10))
            do nn = 1,niso
              ! compute equilibrium fractionation factor relevant to rain temperature
              alpha_equil_rain(nn) = &
                   alfaW_equilibrium(tabs_rain_surface,iso_index(nn))
              
              ! equilibrium fractionation factor for ambient air temp.
              updated_temp = temp(k)-DT*lvap(k)*ocp(k)*prv_rev(k)*(1-IFDRY)
              alpha_equil_liq(nn) = alfaW_equilibrium( &
                   updated_temp, iso_index(nn))
!!$              alpha_equil_liq(nn) = alfaW_equilibrium( &
!!$                   temp(k)-DT*lvap(k)*ocp(k)*prv_rev(k)*(1-IFDRY), & !updated temp
!!$                   iso_index(nn))
              ! compute ventilation factor for heavy isotope
              ventilation_factor_heavy = (t1_qr_ev*ilamr(k)**cre(10) &
                   + (iso_Sc3(nn)/Sc3)*t2_qr_ev*vsc2(k)*rhof2(k) &
                       *((lamr+0.5*fv_r)**(-cre(11)))) &
                   /(crg(10)*ilamr(k)**cre(10))
              ventilation_ratio = MAX(1.,ventilation_factor_light) &
                   /MAX(1.,ventilation_factor_heavy)
              if(rc(k).gt.R1) then
                ! compute isotope ratio for vapor in equilibrium with
                !   cloud liquid
                isor_v(k,nn) = (rho(k)*iso_qv(k,nn)+iso_rc(k,nn)) &
                     /(rho(k)*qv(k)+alpha_equil_liq(nn)*rc(k))
              else
                isor_v(k,nn) = iso_qv(k,nn)/qv(k)
              end if
              ! compute diffusivity of heavy isotopologue.
              diffu_heavy = diffu(k)/Drat_light_over_heavy(iso_index(nn))
              iso_prv_rev(k,nn) = 2.*pi*diffu_heavy*N0_r(k) &
                   * (t1_qr_ev*ilamr(k)**cre(10) &
                        + (iso_Sc3(nn)/Sc3)*t2_qr_ev*vsc2(k)*rhof2(k) &
                          *((lamr+0.5*fv_r)**(-cre(11)))) &
                   *( rv_rain_surface*isor_r(k,nn)/alpha_equil_rain(nn) &
                     - rho(k)*qv(k)*isor_v(k,nn) )
              ! normalize by rho to give units of kg/kg/s
              iso_prv_rev(k,nn) = iso_prv_rev(k,nn)/rho(k)
              ! multiply by eva_factor to be consistent with any limiting
              !   of prv_rev(k) above.
              iso_prv_rev(k,nn) = eva_factor*iso_prv_rev(k,nn)
              ! Make sure we don't remove more than the existing rain heavy
              !   isotopologue.
              iso_prv_rev(k,nn) = MAX(-iso_qv(k,nn)*odt, &
                   MIN(iso_rr(k,nn)/rho(k)*odts, iso_prv_rev(k,nn)) )
              !bloss: Need to limit to prevent over-equilibration.
              !  Find the mass mixing ratio of heavy vapor if it were in isotopic
              !    equilibrium with vapor.
              iso_qtot = ( rho(k)*iso_qv(k,nn) +iso_rc(k,nn) + iso_rr(k,nn) )/rho(k)
              qcond = ( rc(k) + rr(k) - DT*rho(k)*prv_rev(k) )/rho(k)
              qvapor = qv(k) + DT*prv_rev(k)
              iso_qr_equil = (rr(k) - DT*rho(k)*prv_rev(k) )/rho(k) &
                   * alpha_equil_liq(nn) * iso_qtot &
                   / ( qvapor + alpha_equil_liq(nn)*qcond )
             
              if( (iso_prv_rev(k,nn)*(iso_rr(k,nn)-rho(k)*iso_qr_equil) .gt. 0.) &
                   .OR.(rc(k).gt.R1)) then
                ! if the isotopic tendency due to evaporation is pushing the
                !   rain towards equilibrium, make sure that it doesn't overshoot
                !   the equilibrium.  This is mainly important when running with
                !   longer timesteps.
                iso_rate_max = odts*ABS(iso_rr(k,nn)/rho(k) - iso_qr_equil)
              iso_prv_rev(k,nn) = &
                   MAX(-iso_rate_max, MIN(iso_rate_max, iso_prv_rev(k,nn) ) )
            end if
              !bloss (optional approach): Figure out what the rain
              !  mixing ratio would be if it equilibrated jointly with
              !  vapor and cloud liquid with the vapor.  Then use that
              !  value to limit the isotopic exchange, so that we don't
              !   overshoot equilibrium in one time step.  The
              ! timescales for rain adjustment to equilibrium should
              ! be longer, so that this shouldn't be a big deal, I think.
            end do
          else
            !bloss: Evaporate without fractionation.
            !  NOTE: This is unphysical.  Use only for
            !  testing or senstivity study.
            iso_prv_rev(k,:) = prv_rev(k)*isor_r(k,:)
          end if
!..REDUCE the rain evaporation in same places as melting graupel occurs. 
!..Rationale: not much shedding of the water from the graupel so
!..likely that the water-coated graupel evaporating much slower than
!..if the water was immediately shed off.
          IF (prr_gml(k).gt.0.0) THEN
             eva_factor = MIN(1.0, 0.01+(0.99-0.01)*(tempc/20.0))
             prv_rev(k) = prv_rev(k)*eva_factor
             iso_prv_rev(k,:) = iso_prv_rev(k,:)*eva_factor
          ENDIF
          endif
          pnr_rev(k) = MIN(DBLE(nr(k)*0.99/rho(k)*odts),                &   ! RAIN2M
                       prv_rev(k) * nr(k)/rr(k))
          qrten(k) = qrten(k) - prv_rev(k)
          qvten(k) = qvten(k) + prv_rev(k)
          nrten(k) = nrten(k) - pnr_rev(k)
          tten(k) = tten(k) - lvap(k)*ocp(k)*prv_rev(k)*(1-IFDRY)
          rr(k) = MAX(R1, (qr1d(k) + DT*qrten(k))*rho(k))
          qv(k) = MAX(1.E-10, qv1d(k) + DT*qvten(k))
          nr(k) = MAX(R2, (nr1d(k) + DT*nrten(k))*rho(k))
          temp(k) = t1d(k) + DT*tten(k)
          do nn = 1,niso
            iso_qrten(k,nn) = iso_qrten(k,nn) - iso_prv_rev(k,nn)
            iso_qvten(k,nn) = iso_qvten(k,nn) + iso_prv_rev(k,nn)
            iso_rr(k,nn) = MAX(R1, (iso_qr1d(k,nn) + DT*iso_qrten(k,nn))*rho(k))
            iso_qv(k,nn) = MAX(1.E-10, &
                 iso_qv1d(k,nn) + DT*iso_qvten(k,nn))
          end do
         endif
         if(L_qc(k)) then
           ! equilibrate cloud liquid and vapor here, if any cloud liquid exists.
           do nn = 1,niso
             ! compute the equilibrium fractionation factor over
             ! liquid at the updated temperature
             alpha_equil_liq(nn) = alfaW_equilibrium(temp(k),iso_index(nn))
             ! compute the mass mixing ratio of heavy vapor when in
             !   isotopic equilibrium with cloud.  
             iso_qtot = ( rho(k)*iso_qv(k,nn) + iso_rc(k,nn) )/rho(k)
             qcond = ( rc(k) )/rho(k)
             qvapor = qv(k) 
             iso_qv_equil = &
                  iso_qtot * qvapor / ( qvapor + alpha_equil_liq(nn)*qcond )
             !  The isotopic tendency due to vapor deposition onto cloud
             !    takes the vapor-cloud liquid system into equilibrium in
             !    one time step.  The vapor and cloud liquid should be in 
             !    equilibrium at the end of the time step UNLESS cloud ice
             !    has been instantaneously sublimated after sedimentation.
             !    This process is down below (search for "xri") and will
             !    take the vapor and cloud liquid away from equilibrium.
             iso_prw_vcd(k,nn) = odt*( iso_qv(k,nn) - iso_qv_equil )
           end do
         else ! if no cloud, evaporate all heavy isotopes from cloud liquid
           iso_prw_vcd(k,:) = -(iso_qc1d(k,:)+DT*iso_qcten(k,:))*odt
         end if
         if(.NOT.fractionate) then ! mainly for testing purposes
           if(prw_vcd(k).gt.0.) then
             iso_prw_vcd(k,:) = prw_vcd(k)*isor_v(k,:)
           else
             iso_prw_vcd(k,:) = prw_vcd(k)*isor_c(k,:)
           end if 
         end if
         ! add equilibration tendency to overall vapor/cloud tendency
         iso_qcten(k,:) = iso_qcten(k,:) + iso_prw_vcd(k,:)
         iso_qvten(k,:) = iso_qvten(k,:) - iso_prw_vcd(k,:)
         ! update cloud and vapor isotopologues
         do nn = 1,niso
           iso_rc(k,nn) = MAX(R1, (iso_qc1d(k,nn) + DT*iso_qcten(k,nn))*rho(k))
           iso_qv(k,nn) = MAX(1.E-10, iso_qv1d(k,nn) + DT*iso_qvten(k,nn))
         end do
          rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
       enddo

       !bloss(2014-10): Store aggregate tendency of precipitating species before
       !   sedimentation.  We'll subtract this away from the tendency afterwards
       !   to get the sedimentation tendency.  This is for rain, snow and graupel together.
       qpcp_sed(kts:kte) = qrten(kts:kte) + qsten(kts:kte) + qgten(kts:kte)
       iso_qpcp_sed(kts:kte,:) = iso_qrten(kts:kte,:) + iso_qsten(kts:kte,:) + iso_qgten(kts:kte,:)

      !bloss: flag to turn off sedimentation, mostly for testing purposes
      !  would also allow microphysics to be run as a parcel model.
      if(dosedimentation) then 
!+---+-----------------------------------------------------------------+
!..Find max terminal fallspeed (distribution mass-weighted mean
!.. velocity) and use it to determine if we need to split the timestep
!.. (var nstep>1).  Either way, only bother to do sedimentation below
!.. 1st level that contains any sedimenting particles (k=ksed1 on down).
!.. New in v3.0+ is computing separate for rain, ice, snow, and
!.. graupel species thus making code faster with credit to J. Schmidt.
!+---+-----------------------------------------------------------------+
      nstep = 0
      onstep(:) = 1.0
      ksed1(:) = 1
      do k = kte+1, kts, -1
         vtrk(k) = 0.
         vtnrk(k) = 0.
         vtik(k) = 0.
         vtnik(k) = 0.
         vtsk(k) = 0.
         vtgk(k) = 0.
      enddo
      do k = kte, kts, -1
         vtr = 0.
         rhof(k) = SQRT(RHO_NOT/rho(k))
         if (rr(k).gt. R1) then
          lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
          vtr = rhof(k)*av_r*crg(6)*org3 * lamr**cre(3)                 &
                      *((lamr+fv_r)**(-cre(6)))
          vtrk(k) = vtr
! First below is technically correct:
!         vtr = rhof(k)*av_r*crg(5)*org2 * lamr**cre(2)                 &
!                     *((lamr+fv_r)**(-cre(5)))
! Test: make number fall faster (but still slower than mass)
! Goal: less prominent size sorting
          vtr = rhof(k)*av_r*crg(7)/crg(12) * lamr**cre(12)             &
                      *((lamr+fv_r)**(-cre(7)))
          vtnrk(k) = vtr
         else
          vtrk(k) = vtrk(k+1)
          vtnrk(k) = vtnrk(k+1)
         endif
         if (MAX(vtrk(k),vtnrk(k)) .gt. 1.E-3) then
            ksed1(1) = MAX(ksed1(1), k)
            delta_tp = dzq(k)/(MAX(vtrk(k),vtnrk(k)))
            nstep = MAX(nstep, INT(DT/delta_tp + 1.))
         endif
      enddo
      if (ksed1(1) .eq. kte) ksed1(1) = kte-1
      if (nstep .gt. 0) onstep(1) = 1./REAL(nstep)
!+---+-----------------------------------------------------------------+
      if (.not. iiwarm) then
       nstep = 0
       do k = kte, kts, -1
          vti = 0.
          if (ri(k).gt. R1) then
           lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
           ilami = 1./lami
           vti = rhof(k)*av_i*cig(3)*oig2 * ilami**bv_i
           vtik(k) = vti
! First below is technically correct:
!          vti = rhof(k)*av_i*cig(4)*oig1 * ilami**bv_i
! Goal: less prominent size sorting
           vti = rhof(k)*av_i*cig(6)/cig(7) * ilami**bv_i
           vtnik(k) = vti
          else
           vtik(k) = vtik(k+1)
           vtnik(k) = vtnik(k+1)
          endif
          if (vtik(k) .gt. 1.E-3) then
             ksed1(2) = MAX(ksed1(2), k)
             delta_tp = dzq(k)/vtik(k)
             nstep = MAX(nstep, INT(DT/delta_tp + 1.))
          endif
       enddo
       if (ksed1(2) .eq. kte) ksed1(2) = kte-1
       if (nstep .gt. 0) onstep(2) = 1./REAL(nstep)
!+---+-----------------------------------------------------------------+
       nstep = 0
       do k = kte, kts, -1
          vts = 0.
          if (rs(k).gt. R1) then
           xDs = smoc(k) / smob(k)
           Mrat = 1./xDs
           ils1 = 1./(Mrat*Lam0 + fv_s)
           ils2 = 1./(Mrat*Lam1 + fv_s)
           t1_vts = Kap0*csg(4)*ils1**cse(4)
           t2_vts = Kap1*Mrat**mu_s*csg(10)*ils2**cse(10)
           ils1 = 1./(Mrat*Lam0)
           ils2 = 1./(Mrat*Lam1)
           t3_vts = Kap0*csg(1)*ils1**cse(1)
           t4_vts = Kap1*Mrat**mu_s*csg(7)*ils2**cse(7)
           vts = rhof(k)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)
           if (temp(k).gt. T_0) then
            vtsk(k) = MAX(vts*vts_boost(k), vtrk(k))
           else
            vtsk(k) = vts*vts_boost(k)
           endif
          else
            vtsk(k) = vtsk(k+1)
          endif
          if (vtsk(k) .gt. 1.E-3) then
             ksed1(3) = MAX(ksed1(3), k)
             delta_tp = dzq(k)/vtsk(k)
             nstep = MAX(nstep, INT(DT/delta_tp + 1.))
          endif
       enddo
       if (ksed1(3) .eq. kte) ksed1(3) = kte-1
       if (nstep .gt. 0) onstep(3) = 1./REAL(nstep)
!+---+-----------------------------------------------------------------+
       nstep = 0
       do k = kte, kts, -1
          vtg = 0.
          if (rg(k).gt. R1) then
           vtg = rhof(k)*av_g*cgg(6)*ogg3 * ilamg(k)**bv_g
           if (temp(k).gt. T_0) then
            vtgk(k) = MAX(vtg, vtrk(k))
           else
            vtgk(k) = vtg
           endif
          else
            vtgk(k) = vtgk(k+1)
          endif
          if (vtgk(k) .gt. 1.E-3) then
             ksed1(4) = MAX(ksed1(4), k)
             delta_tp = dzq(k)/vtgk(k)
             nstep = MAX(nstep, INT(DT/delta_tp + 1.))
          endif
       enddo
       if (ksed1(4) .eq. kte) ksed1(4) = kte-1
       if (nstep .gt. 0) onstep(4) = 1./REAL(nstep)
      endif
!+---+-----------------------------------------------------------------+
!..Sedimentation of mixing ratio is the integral of v(D)*m(D)*N(D)*dD,
!.. whereas neglect m(D) term for number concentration.  Therefore,
!.. cloud ice has proper differential sedimentation.
!.. New in v3.0+ is computing separate for rain, ice, snow, and
!.. graupel species thus making code faster with credit to J. Schmidt.
!+---+-----------------------------------------------------------------+
      nstep = NINT(1./onstep(1))
      do n = 1, nstep
         do k = kte, kts, -1
            sed_r(k) = vtrk(k)*rr(k)
            sed_n(k) = vtnrk(k)*nr(k)
            iso_sed_r(k,:) = vtrk(k)*iso_rr(k,:)
         enddo
         k = kte
         odzq = 1./dzq(k)
         orho = 1./rho(k)
         qrten(k) = qrten(k) - sed_r(k)*odzq*onstep(1)*orho
         nrten(k) = nrten(k) - sed_n(k)*odzq*onstep(1)*orho
         rr(k) = MAX(R1, rr(k) - sed_r(k)*odzq*DT*onstep(1))
         nr(k) = MAX(R2, nr(k) - sed_n(k)*odzq*DT*onstep(1))
         do nn = 1,niso
           iso_qrten(k,nn) = iso_qrten(k,nn) &
                - iso_sed_r(k,nn)*odzq*onstep(1)*orho
           iso_rr(k,nn) = MAX( R1, &
                iso_rr(k,nn) - iso_sed_r(k,nn)*odzq*DT*onstep(1) )
         end do
         do k = ksed1(1), kts, -1
            odzq = 1./dzq(k)
            orho = 1./rho(k)
            qrten(k) = qrten(k) + (sed_r(k+1)-sed_r(k)) &
                                               *odzq*onstep(1)*orho
            nrten(k) = nrten(k) + (sed_n(k+1)-sed_n(k)) &
                                               *odzq*onstep(1)*orho
            rr(k) = MAX(R1, rr(k) + (sed_r(k+1)-sed_r(k)) &
                                           *odzq*DT*onstep(1))
            nr(k) = MAX(R2, nr(k) + (sed_n(k+1)-sed_n(k)) &
                                           *odzq*DT*onstep(1))
            do nn = 1,niso
              iso_qrten(k,nn) = iso_qrten(k,nn) &
                   + (iso_sed_r(k+1,nn)-iso_sed_r(k,nn)) &
                      *odzq*onstep(1)*orho
              iso_rr(k,nn) = MAX(R1, &
                   iso_rr(k,nn) + (iso_sed_r(k+1,nn)-iso_sed_r(k,nn)) &
                                  *odzq*DT*onstep(1) )
            end do
         enddo
         if (rr(kts).gt.R1*10.) then
         pptrain = pptrain + sed_r(kts)*DT*onstep(1)
         iso_pptrain(:) = iso_pptrain(:) + iso_sed_r(kts,:)*DT*onstep(1)
         end if
      enddo
!+---+-----------------------------------------------------------------+
      nstep = NINT(1./onstep(2))
      do n = 1, nstep
         do k = kte, kts, -1
            sed_i(k) = vtik(k)*ri(k)
            sed_n(k) = vtnik(k)*ni(k)
            iso_sed_i(k,:) = vtik(k)*iso_ri(k,:)
         enddo
         k = kte
         odzq = 1./dzq(k)
         orho = 1./rho(k)
         qiten(k) = qiten(k) - sed_i(k)*odzq*onstep(2)*orho
         niten(k) = niten(k) - sed_n(k)*odzq*onstep(2)*orho
         ri(k) = MAX(R1, ri(k) - sed_i(k)*odzq*DT*onstep(2))
         ni(k) = MAX(R2, ni(k) - sed_n(k)*odzq*DT*onstep(2))
         iso_qiten(k,:) = iso_qiten(k,:) - iso_sed_i(k,:)*odzq*onstep(2)*orho
         iso_ri(k,:) = iso_ri(k,:) - iso_sed_i(k,:)*odzq*DT*onstep(2)
         do k = ksed1(2), kts, -1
            odzq = 1./dzq(k)
            orho = 1./rho(k)
            qiten(k) = qiten(k) + (sed_i(k+1)-sed_i(k)) &
                                               *odzq*onstep(2)*orho
            niten(k) = niten(k) + (sed_n(k+1)-sed_n(k)) &
                                               *odzq*onstep(2)*orho
            ri(k) = MAX(R1, ri(k) + (sed_i(k+1)-sed_i(k)) &
                                           *odzq*DT*onstep(2))
            ni(k) = MAX(R2, ni(k) + (sed_n(k+1)-sed_n(k)) &
                                           *odzq*DT*onstep(2))
            iso_qiten(k,:) = iso_qiten(k,:) + (iso_sed_i(k+1,:)-iso_sed_i(k,:)) &
                                               *odzq*onstep(2)*orho
            iso_ri(k,:) = iso_ri(k,:) + (iso_sed_i(k+1,:)-iso_sed_i(k,:)) &
                                           *odzq*DT*onstep(2)
            do nn = 1,niso
              iso_ri(k,nn) = MAX(R1, iso_ri(k,nn))
            end do
         enddo
         if (ri(kts).gt.R1*10.) then
         pptice = pptice + sed_i(kts)*DT*onstep(2)
         iso_pptice(:) = iso_pptice(:) + iso_sed_i(kts,:)*DT*onstep(2)
         end if
      enddo
!+---+-----------------------------------------------------------------+
      nstep = NINT(1./onstep(3))
      do n = 1, nstep
         do k = kte, kts, -1
            sed_s(k) = vtsk(k)*rs(k)
            iso_sed_s(k,:) = vtsk(k)*iso_rs(k,:)
        enddo
         k = kte
         odzq = 1./dzq(k)
         orho = 1./rho(k)
         qsten(k) = qsten(k) - sed_s(k)*odzq*onstep(3)*orho
         rs(k) = MAX(R1, rs(k) - sed_s(k)*odzq*DT*onstep(3))
         iso_qsten(k,:) = iso_qsten(k,:) - iso_sed_s(k,:)*odzq*onstep(3)*orho
         iso_rs(k,:) = iso_rs(k,:) - iso_sed_s(k,:)*odzq*DT*onstep(3)
         do k = ksed1(3), kts, -1
            odzq = 1./dzq(k)
            orho = 1./rho(k)
            qsten(k) = qsten(k) + (sed_s(k+1)-sed_s(k)) &
                                               *odzq*onstep(3)*orho
            rs(k) = MAX(R1, rs(k) + (sed_s(k+1)-sed_s(k)) &
                                           *odzq*DT*onstep(3))
            iso_qsten(k,:) = iso_qsten(k,:) + (iso_sed_s(k+1,:)-iso_sed_s(k,:)) &
                                               *odzq*onstep(3)*orho
            iso_rs(k,:) = iso_rs(k,:) + (iso_sed_s(k+1,:)-iso_sed_s(k,:)) &
                                           *odzq*DT*onstep(3)
            do nn = 1,niso
              iso_rs(k,nn) = MAX(R1, iso_rs(k,nn))
            end do
         enddo
         if (rs(kts).gt.R1*10.) then
         pptsnow = pptsnow + sed_s(kts)*DT*onstep(3)
         iso_pptsnow(:) = iso_pptsnow(:) + iso_sed_s(kts,:)*DT*onstep(3)
         end if
      enddo
!+---+-----------------------------------------------------------------+
      nstep = NINT(1./onstep(4))
      do n = 1, nstep
         do k = kte, kts, -1
            sed_g(k) = vtgk(k)*rg(k)
            iso_sed_g(k,:) = vtgk(k)*iso_rg(k,:)
         enddo
         k = kte
         odzq = 1./dzq(k)
         orho = 1./rho(k)
         qgten(k) = qgten(k) - sed_g(k)*odzq*onstep(4)*orho
         rg(k) = MAX(R1, rg(k) - sed_g(k)*odzq*DT*onstep(4))
         iso_qgten(k,:) = iso_qgten(k,:) - iso_sed_g(k,:)*odzq*onstep(4)*orho
         iso_rg(k,:) = iso_rg(k,:) - iso_sed_g(k,:)*odzq*DT*onstep(4)
         do k = ksed1(4), kts, -1
            odzq = 1./dzq(k)
            orho = 1./rho(k)
            qgten(k) = qgten(k) + (sed_g(k+1)-sed_g(k)) &
                                               *odzq*onstep(4)*orho
            rg(k) = MAX(R1, rg(k) + (sed_g(k+1)-sed_g(k)) &
                                           *odzq*DT*onstep(4))
            iso_qgten(k,:) = iso_qgten(k,:) + (iso_sed_g(k+1,:)-iso_sed_g(k,:)) &
                                               *odzq*onstep(4)*orho
            iso_rg(k,:) = iso_rg(k,:) + (iso_sed_g(k+1,:)-iso_sed_g(k,:)) &
                                           *odzq*DT*onstep(4)
            do nn = 1,niso
              iso_rg(k,nn) = MAX(R1, iso_rg(k,nn))
            end do
         enddo
         if (rg(kts).gt.R1*10.) then
         pptgraul = pptgraul + sed_g(kts)*DT*onstep(4)
         iso_pptgraul(:) = iso_pptgraul(:) + iso_sed_g(kts,:)*DT*onstep(4)
         end if
      enddo
    end if

       !bloss(2014-10): Store aggregate tendency of precipitating species before
       !   sedimentation.  We'll subtract this away from the tendency afterwards
       !   to get the sedimentation tendency.  This is for rain, snow and graupel together.
       qpcp_sed(kts:kte) = qrten(kts:kte) + qsten(kts:kte) + qgten(kts:kte) &
            - qpcp_sed(kts:kte)
       iso_qpcp_sed(kts:kte,:) = iso_qrten(kts:kte,:) + iso_qsten(kts:kte,:) + iso_qgten(kts:kte,:) &
            - iso_qpcp_sed(kts:kte,:)

!+---+-----------------------------------------------------------------+
!.. Instantly melt any cloud ice into cloud water if above 0C and
!.. instantly freeze any cloud water found below HGFR.
!+---+-----------------------------------------------------------------+
      if (.not. iiwarm) then
      do k = kts, kte
         xri = MAX(0.0, qi1d(k) + qiten(k)*DT)
         if ( (temp(k).gt. T_0) .and. (xri.gt. 0.0) ) then
          qcten(k) = qcten(k) + xri*odt
          qiten(k) = qiten(k) - xri*odt
          niten(k) = -ni1d(k)*odt
          tten(k) = tten(k) - lfus*ocp(k)*xri*odt*(1-IFDRY)
          do nn = 1,niso
            xri = MAX(0.0, iso_qi1d(k,nn) + iso_qiten(k,nn)*DT)
            iso_qcten(k,nn) = iso_qcten(k,nn) + xri*odt
            iso_qiten(k,nn) = iso_qiten(k,nn) - xri*odt
          end do
         endif
         xrc = MAX(0.0, qc1d(k) + qcten(k)*DT)
         if ( (temp(k).lt. HGFR) .and. (xrc.gt. 0.0) ) then
          lfus2 = lsub - lvap(k)
          qiten(k) = qiten(k) + xrc*odt
          niten(k) = niten(k) + xrc/xm0i * odt
          qcten(k) = qcten(k) - xrc*odt
          tten(k) = tten(k) + lfus2*ocp(k)*xrc*odt*(1-IFDRY)
          do nn = 1,niso
            xrc = MAX(0.0, iso_qc1d(k,nn) + iso_qcten(k,nn)*DT)
            iso_qiten(k,nn) = iso_qiten(k,nn) + xrc*odt
            iso_qcten(k,nn) = iso_qcten(k,nn) - xrc*odt
          end do
         endif
      enddo
      endif
!+---+-----------------------------------------------------------------+
!.. All tendencies computed, apply and pass back final values to parent.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         t1d(k)  = t1d(k) + tten(k)*DT
         qv1d(k) = MAX(1.E-10, qv1d(k) + qvten(k)*DT)
         do nn = 1,niso
           iso_qv1d(k,nn) = MAX(1.E-10, iso_qv1d(k,nn) + iso_qvten(k,nn)*DT)
         end do
         qc1d(k) = qc1d(k) + qcten(k)*DT
         iso_qc1d(k,:) = iso_qc1d(k,:) + iso_qcten(k,:)*DT
         if (qc1d(k) .le. R1) then
           qc1d(k) = 0.0
           iso_qv1d(k,:) = iso_qv1d(k,:) + iso_qc1d(k,:)
           iso_qc1d(k,:) = 0.
         end if
         qi1d(k) = qi1d(k) + qiten(k)*DT
         ni1d(k) = MAX(R2/rho(k), ni1d(k) + niten(k)*DT)
         iso_qi1d(k,:) = iso_qi1d(k,:) + iso_qiten(k,:)*DT
         if (qi1d(k) .le. R1) then
           qi1d(k) = 0.0
           ni1d(k) = 0.0
           iso_qv1d(k,:) = iso_qv1d(k,:) + iso_qi1d(k,:)
           iso_qi1d(k,:) = 0.
         else
           lami = (am_i*cig(2)*oig1*ni1d(k)/qi1d(k))**obmi
           ilami = 1./lami
           xDi = (bm_i + mu_i + 1.) * ilami
           if (xDi.lt. 20.E-6) then
            lami = cie(2)/20.E-6
           elseif (xDi.gt. 300.E-6) then
            lami = cie(2)/300.E-6
           endif
           ni1d(k) = MIN(cig(1)*oig2*qi1d(k)/am_i*lami**bm_i,           &
                         250.D3/rho(k))
         endif
         qr1d(k) = qr1d(k) + qrten(k)*DT
         nr1d(k) = MAX(R2/rho(k), nr1d(k) + nrten(k)*DT)
         iso_qr1d(k,:) = iso_qr1d(k,:) + iso_qrten(k,:)*DT
         if (qr1d(k) .le. R1) then
           qr1d(k) = 0.0
           nr1d(k) = 0.0
           iso_qv1d(k,:) = iso_qv1d(k,:) + iso_qr1d(k,:)
           iso_qr1d(k,:) = 0.
         else
           lamr = (am_r*crg(3)*org2*nr1d(k)/qr1d(k))**obmr
           mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
           if (mvd_r(k) .gt. 2.5E-3) then
              mvd_r(k) = 2.5E-3
           elseif (mvd_r(k) .lt. D0r*0.75) then
              mvd_r(k) = D0r*0.75
           endif
           lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
           nr1d(k) = crg(2)*org3*qr1d(k)*lamr**bm_r / am_r
         endif
         qs1d(k) = qs1d(k) + qsten(k)*DT
         iso_qs1d(k,:) = iso_qs1d(k,:) + iso_qsten(k,:)*DT
         if (qs1d(k) .le. R1) then
           qs1d(k) = 0.0
           iso_qv1d(k,:) = iso_qv1d(k,:) + iso_qs1d(k,:)
           iso_qs1d(k,:) = 0.
         end if
         qg1d(k) = qg1d(k) + qgten(k)*DT
         iso_qg1d(k,:) = iso_qg1d(k,:) + iso_qgten(k,:)*DT
         if (qg1d(k) .le. R1) then
           qg1d(k) = 0.0
           iso_qv1d(k,:) = iso_qv1d(k,:) + iso_qg1d(k,:)
           iso_qg1d(k,:) = 0.
         end if
      enddo

      !bloss(2014-10): note at this point that all p**_*** tendencies are scaled by
      ! rho except for prw_vcd, prv_rev and prw_iml.
      !  Save then as kg/kg/s here.

      !bloss(2014-10): save aggregated tendencies for precipitating species 
      ! (rain/snow/graupel) to compute budgets for contributions to surface precip.
      ! I realize that cloud ice may precipitate to the surface, but I'm assuming that
      ! this will be a small term.


      !bloss(dep/sub/revp): Separate tendencies for deposition, sublimation and rain evaporation
      !  All based on the sign of deposition/sublimation for standard water.
      where(prs_ide(:) > 0.)
        qpcp_dep(:) = qpcp_dep(:) + prs_ide(:)/rho(:)
      elsewhere
        qpcp_sub(:) = qpcp_sub(:) + prs_ide(:)/rho(:)
      end where

      where(prs_sde(:) > 0.)
        qpcp_dep(:) = qpcp_dep(:) + prs_sde(:)/rho(:)
      elsewhere
        qpcp_sub(:) = qpcp_sub(:) + prs_sde(:)/rho(:)
      end where

      where(prg_gde(:) > 0.)
        qpcp_dep(:) = qpcp_dep(:) + prg_gde(:)/rho(:)
      elsewhere
        qpcp_sub(:) = qpcp_sub(:) + prg_gde(:)/rho(:)
      end where

      qpcp_revp(:) = - prv_rev(:)
        
      qvap_cld(:) = prw_vcd(:) ! vapor--> cloud liquid, may help us identify Bergeron process
      qliq_auto(:) = prr_wau(:) / rho(:) ! cloud liq --> rain autoconversion
      qliq_coll(:) = prr_rcw(:) / rho(:) ! collection of cloud liq by rain
      qice_auto(:) = prs_iau(:) / rho(:) ! cloud ice --> snow autoconversion
      qice_coll(:) = prs_sci(:) / rho(:) ! collection of cloud ice by snow
      qice_rim(:) = ( prs_scw(:) + prg_scw(:) + prg_gcw(:) &
           + pri_rci(:) - pri_rfz(:) - pri_ihm(:) ) / rho(:)

      !bloss(2014-10): save aggregated tendencies to compute budgets for liquid/ice precipitation
      do nn = 1,niso
        !bloss(dep/sub/revp): Separate tendencies for deposition, sublimation and rain evaporation
        !  All based on the sign of deposition/sublimation for standard water.
        where(prs_ide(:) > 0.) 
          iso_qpcp_dep(:,nn) = iso_qpcp_dep(:,nn) + iso_prs_ide(:,nn)/rho(:)
        elsewhere
          iso_qpcp_sub(:,nn) = iso_qpcp_sub(:,nn) + iso_prs_ide(:,nn)/rho(:)
        end where

        where(prs_sde(:) > 0.)
          iso_qpcp_dep(:,nn) = iso_qpcp_dep(:,nn) + iso_prs_sde(:,nn)/rho(:)
        elsewhere
          iso_qpcp_sub(:,nn) = iso_qpcp_sub(:,nn) + iso_prs_sde(:,nn)/rho(:)
        end where

        where(prg_gde(:) > 0.)
          iso_qpcp_dep(:,nn) = iso_qpcp_dep(:,nn) + iso_prg_gde(:,nn)/rho(:)
        elsewhere
          iso_qpcp_sub(:,nn) = iso_qpcp_sub(:,nn) + iso_prg_gde(:,nn)/rho(:)
        end where

        iso_qpcp_revp(:,nn) = - iso_prv_rev(:,nn)

        iso_qvap_cld(:,nn) = iso_prw_vcd(:,nn)
        iso_qliq_auto(:,nn) = iso_prr_wau(:,nn) / rho(:) ! cloud liq --> rain autoconversion
        iso_qliq_coll(:,nn) = iso_prr_rcw(:,nn) / rho(:) ! collection of cloud liq by rain
        iso_qice_auto(:,nn) = iso_prs_iau(:,nn) / rho(:) ! cloud ice --> snow autoconversion
        iso_qice_coll(:,nn) = iso_prs_sci(:,nn) / rho(:) ! collection of cloud ice by snow
        iso_qice_rim(:,nn) = ( iso_prs_scw(:,nn) + iso_prg_scw(:,nn) + iso_prg_gcw(:,nn) &
             + iso_pri_rci(:,nn) - iso_pri_rfz(:,nn) - iso_pri_ihm(:,nn) ) / rho(:)
      end do

      end subroutine mp_thompiso
!+---+-----------------------------------------------------------------+
!
!+---+-----------------------------------------------------------------+
!..Creation of the lookup tables and support functions found below here.
!+---+-----------------------------------------------------------------+
!..Rain collecting graupel (and inverse).  Explicit CE integration.
!+---+-----------------------------------------------------------------+
      subroutine qr_acr_qg
      implicit none
!..Local variables
      INTEGER:: i, j, k, m, n, n2
      INTEGER:: km, km_s, km_e
      DOUBLE PRECISION, DIMENSION(nbg):: vg, N_g
      DOUBLE PRECISION, DIMENSION(nbr):: vr, N_r
      DOUBLE PRECISION:: N0_r, N0_g, lam_exp, lamg, lamr
      DOUBLE PRECISION:: massg, massr, dvg, dvr, t1, t2, z1, z2, y1, y2
!+---+
      do n2 = 1, nbr
!        vr(n2) = av_r*Dr(n2)**bv_r * DEXP(-fv_r*Dr(n2))
         vr(n2) = -0.1021 + 4.932E3*Dr(n2) - 0.9551E6*Dr(n2)*Dr(n2)     &
              + 0.07934E9*Dr(n2)*Dr(n2)*Dr(n2)                          &
              - 0.002362E12*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
      enddo
      do n = 1, nbg
         vg(n) = av_g*Dg(n)**bv_g
      enddo
!..Note values returned from wrf_dm_decomp1d are zero-based, add 1 for
!.. fortran indices.  J. Michalakes, 2009Oct30.
      km_s = 0
      km_e = ntb_r*ntb_r1 - 1
      do km = km_s, km_e
         m = km / ntb_r1 + 1
         k = mod( km , ntb_r1 ) + 1
         lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
         lamr = lam_exp * (crg(3)*org2*org1)**obmr
         N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
         do n2 = 1, nbr
            N_r(n2) = N0_r*Dr(n2)**mu_r *DEXP(-lamr*Dr(n2))*dtr(n2)
         enddo
         do j = 1, ntb_g
         do i = 1, ntb_g1
            lam_exp = (N0g_exp(i)*am_g*cgg(1)/r_g(j))**oge1
            lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
            N0_g = N0g_exp(i)/(cgg(2)*lam_exp) * lamg**cge(2)
            do n = 1, nbg
               N_g(n) = N0_g*Dg(n)**mu_g * DEXP(-lamg*Dg(n))*dtg(n)
            enddo
            t1 = 0.0d0
            t2 = 0.0d0
            z1 = 0.0d0
            z2 = 0.0d0
            y1 = 0.0d0
            y2 = 0.0d0
            do n2 = 1, nbr
               massr = am_r * Dr(n2)**bm_r
               do n = 1, nbg
                  massg = am_g * Dg(n)**bm_g
                  dvg = 0.5d0*((vr(n2) - vg(n)) + DABS(vr(n2)-vg(n)))
                  dvr = 0.5d0*((vg(n) - vr(n2)) + DABS(vg(n)-vr(n2)))
                  t1 = t1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                      *dvg*massg * N_g(n)* N_r(n2)
                  z1 = z1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                      *dvg*massr * N_g(n)* N_r(n2)
                  y1 = y1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                      *dvg       * N_g(n)* N_r(n2)
                  t2 = t2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                      *dvr*massr * N_g(n)* N_r(n2)
                  y2 = y2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                      *dvr       * N_g(n)* N_r(n2)
                  z2 = z2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                      *dvr*massg * N_g(n)* N_r(n2)
               enddo
 97            continue
            enddo
            tcg_racg(i,j,k,m) = t1
            tmr_racg(i,j,k,m) = DMIN1(z1, r_r(m)*1.0d0)
            tcr_gacr(i,j,k,m) = t2
            tmg_gacr(i,j,k,m) = z2
            tnr_racg(i,j,k,m) = y1
            tnr_gacr(i,j,k,m) = y2
         enddo
         enddo
      enddo
!..Note wrf_dm_gatherv expects zero-based km_s, km_e (J. Michalakes, 2009Oct30).
      end subroutine qr_acr_qg
!+---+-----------------------------------------------------------------+
!
!+---+-----------------------------------------------------------------+
!..Rain collecting snow (and inverse).  Explicit CE integration.
!+---+-----------------------------------------------------------------+
      subroutine qr_acr_qs
      implicit none
!..Local variables
      INTEGER:: i, j, k, m, n, n2
      INTEGER:: km, km_s, km_e
      DOUBLE PRECISION, DIMENSION(nbr):: vr, D1, N_r
      DOUBLE PRECISION, DIMENSION(nbs):: vs, N_s
      DOUBLE PRECISION:: loga_, a_, b_, second, M0, M2, M3, Mrat, oM3
      DOUBLE PRECISION:: N0_r, lam_exp, lamr, slam1, slam2
      DOUBLE PRECISION:: dvs, dvr, masss, massr
      DOUBLE PRECISION:: t1, t2, t3, t4, z1, z2, z3, z4
      DOUBLE PRECISION:: y1, y2, y3, y4
!+---+
      do n2 = 1, nbr
!        vr(n2) = av_r*Dr(n2)**bv_r * DEXP(-fv_r*Dr(n2))
         vr(n2) = -0.1021 + 4.932E3*Dr(n2) - 0.9551E6*Dr(n2)*Dr(n2)     &
              + 0.07934E9*Dr(n2)*Dr(n2)*Dr(n2)                          &
              - 0.002362E12*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
         D1(n2) = (vr(n2)/av_s)**(1./bv_s)
      enddo
      do n = 1, nbs
         vs(n) = 1.5*av_s*Ds(n)**bv_s * DEXP(-fv_s*Ds(n))
      enddo
!..Note values returned from wrf_dm_decomp1d are zero-based, add 1 for
!.. fortran indices.  J. Michalakes, 2009Oct30.
      km_s = 0
      km_e = ntb_r*ntb_r1 - 1
      do km = km_s, km_e
         m = km / ntb_r1 + 1
         k = mod( km , ntb_r1 ) + 1
         lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
         lamr = lam_exp * (crg(3)*org2*org1)**obmr
         N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
         do n2 = 1, nbr
            N_r(n2) = N0_r*Dr(n2)**mu_r * DEXP(-lamr*Dr(n2))*dtr(n2)
         enddo
         do j = 1, ntb_t
            do i = 1, ntb_s
!..From the bm_s moment, compute plus one moment.  If we are not
!.. using bm_s=2, then we must transform to the pure 2nd moment
!.. (variable called "second") and then to the bm_s+1 moment.
               M2 = r_s(i)*oams *1.0d0
               if (bm_s.gt.2.0-1.E-3 .and. bm_s.lt.2.0+1.E-3) then
                  loga_ = sa(1) + sa(2)*Tc(j) + sa(3)*bm_s &
                     + sa(4)*Tc(j)*bm_s + sa(5)*Tc(j)*Tc(j) &
                     + sa(6)*bm_s*bm_s + sa(7)*Tc(j)*Tc(j)*bm_s &
                     + sa(8)*Tc(j)*bm_s*bm_s + sa(9)*Tc(j)*Tc(j)*Tc(j) &
                     + sa(10)*bm_s*bm_s*bm_s
                  a_ = 10.0**loga_
                  b_ = sb(1) + sb(2)*Tc(j) + sb(3)*bm_s &
                     + sb(4)*Tc(j)*bm_s + sb(5)*Tc(j)*Tc(j) &
                     + sb(6)*bm_s*bm_s + sb(7)*Tc(j)*Tc(j)*bm_s &
                     + sb(8)*Tc(j)*bm_s*bm_s + sb(9)*Tc(j)*Tc(j)*Tc(j) &
                     + sb(10)*bm_s*bm_s*bm_s
                  second = (M2/a_)**(1./b_)
               else
                  second = M2
               endif
               loga_ = sa(1) + sa(2)*Tc(j) + sa(3)*cse(1) &
                  + sa(4)*Tc(j)*cse(1) + sa(5)*Tc(j)*Tc(j) &
                  + sa(6)*cse(1)*cse(1) + sa(7)*Tc(j)*Tc(j)*cse(1) &
                  + sa(8)*Tc(j)*cse(1)*cse(1) + sa(9)*Tc(j)*Tc(j)*Tc(j) &
                  + sa(10)*cse(1)*cse(1)*cse(1)
               a_ = 10.0**loga_
               b_ = sb(1)+sb(2)*Tc(j)+sb(3)*cse(1) + sb(4)*Tc(j)*cse(1) &
                  + sb(5)*Tc(j)*Tc(j) + sb(6)*cse(1)*cse(1) &
                  + sb(7)*Tc(j)*Tc(j)*cse(1) + sb(8)*Tc(j)*cse(1)*cse(1) &
                  + sb(9)*Tc(j)*Tc(j)*Tc(j)+sb(10)*cse(1)*cse(1)*cse(1)
               M3 = a_ * second**b_
               oM3 = 1./M3
               Mrat = M2*(M2*oM3)*(M2*oM3)*(M2*oM3)
               M0   = (M2*oM3)**mu_s
               slam1 = M2 * oM3 * Lam0
               slam2 = M2 * oM3 * Lam1
               do n = 1, nbs
                  N_s(n) = Mrat*(Kap0*DEXP(-slam1*Ds(n)) &
                      + Kap1*M0*Ds(n)**mu_s * DEXP(-slam2*Ds(n)))*dts(n)
               enddo
               t1 = 0.0d0
               t2 = 0.0d0
               t3 = 0.0d0
               t4 = 0.0d0
               z1 = 0.0d0
               z2 = 0.0d0
               z3 = 0.0d0
               z4 = 0.0d0
               y1 = 0.0d0
               y2 = 0.0d0
               y3 = 0.0d0
               y4 = 0.0d0
               do n2 = 1, nbr
                  massr = am_r * Dr(n2)**bm_r
                  do n = 1, nbs
                     masss = am_s * Ds(n)**bm_s
      
                     dvs = 0.5d0*((vr(n2) - vs(n)) + DABS(vr(n2)-vs(n)))
                     dvr = 0.5d0*((vs(n) - vr(n2)) + DABS(vs(n)-vr(n2)))
                     if (massr .gt. 1.5*masss) then
                     t1 = t1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvs*masss * N_s(n)* N_r(n2)
                     z1 = z1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvs*massr * N_s(n)* N_r(n2)
                     y1 = y1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvs       * N_s(n)* N_r(n2)
                     else
                     t3 = t3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvs*masss * N_s(n)* N_r(n2)
                     z3 = z3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvs*massr * N_s(n)* N_r(n2)
                     y3 = y3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvs       * N_s(n)* N_r(n2)
                     endif
                     if (massr .gt. 1.5*masss) then
                     t2 = t2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvr*massr * N_s(n)* N_r(n2)
                     y2 = y2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvr       * N_s(n)* N_r(n2)
                     z2 = z2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvr*masss * N_s(n)* N_r(n2)
                     else
                     t4 = t4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvr*massr * N_s(n)* N_r(n2)
                     y4 = y4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvr       * N_s(n)* N_r(n2)
                     z4 = z4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                         *dvr*masss * N_s(n)* N_r(n2)
                     endif
                  enddo
               enddo
               tcs_racs1(i,j,k,m) = t1
               tmr_racs1(i,j,k,m) = DMIN1(z1, r_r(m)*1.0d0)
               tcs_racs2(i,j,k,m) = t3
               tmr_racs2(i,j,k,m) = z3
               tcr_sacr1(i,j,k,m) = t2
               tms_sacr1(i,j,k,m) = z2
               tcr_sacr2(i,j,k,m) = t4
               tms_sacr2(i,j,k,m) = z4
               tnr_racs1(i,j,k,m) = y1
               tnr_racs2(i,j,k,m) = y3
               tnr_sacr1(i,j,k,m) = y2
               tnr_sacr2(i,j,k,m) = y4
            enddo
         enddo
      enddo
!..Note wrf_dm_gatherv expects zero-based km_s, km_e (J. Michalakes, 2009Oct30).
      end subroutine qr_acr_qs
!+---+-----------------------------------------------------------------+
!
!+---+-----------------------------------------------------------------+
!..This is a literal adaptation of Bigg (1954) probability of drops of
!..a particular volume freezing.  Given this probability, simply freeze
!..the proportion of drops summing their masses.
!+---+-----------------------------------------------------------------+
      subroutine freezeH2O
      implicit none
!..Local variables
      INTEGER:: i, j, k, n, n2
      DOUBLE PRECISION, DIMENSION(nbr):: N_r, massr
      DOUBLE PRECISION, DIMENSION(nbc):: N_c, massc
      DOUBLE PRECISION:: sum1, sum2, sumn1, sumn2, &
                         prob, vol, Texp, orho_w, &
                         lam_exp, lamr, N0_r, lamc, N0_c, y
!+---+
      orho_w = 1./rho_w
      do n2 = 1, nbr
         massr(n2) = am_r*Dr(n2)**bm_r
      enddo
      do n = 1, nbc
         massc(n) = am_r*Dc(n)**bm_r
      enddo
!..Freeze water (smallest drops become cloud ice, otherwise graupel).
      do k = 1, 45
!         print*, ' Freezing water for temp = ', -k
         Texp = DEXP( DFLOAT(k) ) - 1.0D0
         do j = 1, ntb_r1
            do i = 1, ntb_r
               lam_exp = (N0r_exp(j)*am_r*crg(1)/r_r(i))**ore1
               lamr = lam_exp * (crg(3)*org2*org1)**obmr
               N0_r = N0r_exp(j)/(crg(2)*lam_exp) * lamr**cre(2)
               sum1 = 0.0d0
               sum2 = 0.0d0
               sumn1 = 0.0d0
               sumn2 = 0.0d0
               do n2 = nbr, 1, -1
                  N_r(n2) = N0_r*Dr(n2)**mu_r*DEXP(-lamr*Dr(n2))*dtr(n2)
                  vol = massr(n2)*orho_w
                  prob = 1.0D0 - DEXP(-120.0D0*vol*5.2D-4 * Texp)
                  if (massr(n2) .lt. xm0g) then
                     sumn1 = sumn1 + prob*N_r(n2)
                     sum1 = sum1 + prob*N_r(n2)*massr(n2)
                  else
                     sumn2 = sumn2 + prob*N_r(n2)
                     sum2 = sum2 + prob*N_r(n2)*massr(n2)
                  endif
                  if ((sum1+sum2) .ge. r_r(i)) EXIT
               enddo
               tpi_qrfz(i,j,k) = sum1
               tni_qrfz(i,j,k) = sumn1
               tpg_qrfz(i,j,k) = sum2
               tnr_qrfz(i,j,k) = sumn2
            enddo
         enddo
         do i = 1, ntb_c
            lamc = 1.0D-6 * (Nt_c*am_r* ccg(2) * ocg1 / r_c(i))**obmr
            N0_c = 1.0D-18 * Nt_c*ocg1 * lamc**cce(1)
            sum1 = 0.0d0
            sumn2 = 0.0d0
            do n = nbc, 1, -1
               y = Dc(n)*1.0D6
               vol = massc(n)*orho_w
               prob = 1.0D0 - DEXP(-120.0D0*vol*5.2D-4 * Texp)
               N_c(n) = N0_c* y**mu_c * EXP(-lamc*y)*dtc(n)
               N_c(n) = 1.0D24 * N_c(n)
               sumn2 = sumn2 + prob*N_c(n)
               sum1 = sum1 + prob*N_c(n)*massc(n)
               if (sum1 .ge. r_c(i)) EXIT
            enddo
            tpi_qcfz(i,k) = sum1
            tni_qcfz(i,k) = sumn2
         enddo
      enddo
      end subroutine freezeH2O
!+---+-----------------------------------------------------------------+
!
!+---+-----------------------------------------------------------------+
!..Cloud ice converting to snow since portion greater than min snow
!.. size.  Given cloud ice content (kg/m**3), number concentration
!.. (#/m**3) and gamma shape parameter, mu_i, break the distrib into
!.. bins and figure out the mass/number of ice with sizes larger than
!.. D0s.  Also, compute incomplete gamma function for the integration
!.. of ice depositional growth from diameter=0 to D0s.  Amount of
!.. ice depositional growth is this portion of distrib while larger
!.. diameters contribute to snow growth (as in Harrington et al. 1995).
!+---+-----------------------------------------------------------------+
      subroutine qi_aut_qs
      implicit none
!..Local variables
      INTEGER:: i, j, n2
      DOUBLE PRECISION, DIMENSION(nbi):: N_i
      DOUBLE PRECISION:: N0_i, lami, Di_mean, t1, t2
      REAL:: xlimit_intg
!+---+
      do j = 1, ntb_i1
         do i = 1, ntb_i
            lami = (am_i*cig(2)*oig1*Nt_i(j)/r_i(i))**obmi
            Di_mean = (bm_i + mu_i + 1.) / lami
            N0_i = Nt_i(j)*oig1 * lami**cie(1)
            t1 = 0.0d0
            t2 = 0.0d0
            if (SNGL(Di_mean) .gt. 5.*D0s) then
             t1 = r_i(i)
             t2 = Nt_i(j)
             tpi_ide(i,j) = 0.0D0
            elseif (SNGL(Di_mean) .lt. D0i) then
             t1 = 0.0D0
             t2 = 0.0D0
             tpi_ide(i,j) = 1.0D0
            else
             xlimit_intg = lami*D0s
             tpi_ide(i,j) = GAMMP(mu_i+2.0, xlimit_intg) * 1.0D0
             do n2 = 1, nbi
               N_i(n2) = N0_i*Di(n2)**mu_i * DEXP(-lami*Di(n2))*dti(n2)
               if (Di(n2).ge.D0s) then
                  t1 = t1 + N_i(n2) * am_i*Di(n2)**bm_i
                  t2 = t2 + N_i(n2)
               endif
             enddo
            endif
            tps_iaus(i,j) = t1
            tni_iaus(i,j) = t2
         enddo
      enddo
      end subroutine qi_aut_qs
!
!+---+-----------------------------------------------------------------+
!..Variable collision efficiency for rain collecting cloud water using
!.. method of Beard and Grover, 1974 if a/A less than 0.25; otherwise
!.. uses polynomials to get close match of Pruppacher & Klett Fig 14-9.
!+---+-----------------------------------------------------------------+
      subroutine table_Efrw
      implicit none
!..Local variables
      DOUBLE PRECISION:: vtr, stokes, reynolds, Ef_rw
      DOUBLE PRECISION:: p, yc0, F, G, H, z, K0, X
      INTEGER:: i, j
      do j = 1, nbc
      do i = 1, nbr
         Ef_rw = 0.0
         p = Dc(j)/Dr(i)
         if (Dr(i).lt.50.E-6 .or. Dc(j).lt.3.E-6) then
          t_Efrw(i,j) = 0.0
         elseif (p.gt.0.25) then
          X = Dc(j)*1.D6
          if (Dr(i) .lt. 75.e-6) then
             Ef_rw = 0.026794*X - 0.20604
          elseif (Dr(i) .lt. 125.e-6) then
             Ef_rw = -0.00066842*X*X + 0.061542*X - 0.37089
          elseif (Dr(i) .lt. 175.e-6) then
             Ef_rw = 4.091e-06*X*X*X*X - 0.00030908*X*X*X               &
                   + 0.0066237*X*X - 0.0013687*X - 0.073022
          elseif (Dr(i) .lt. 250.e-6) then
             Ef_rw = 9.6719e-5*X*X*X - 0.0068901*X*X + 0.17305*X        &
                   - 0.65988
          elseif (Dr(i) .lt. 350.e-6) then
             Ef_rw = 9.0488e-5*X*X*X - 0.006585*X*X + 0.16606*X         &
                   - 0.56125
          else
             Ef_rw = 0.00010721*X*X*X - 0.0072962*X*X + 0.1704*X        &
                   - 0.46929
          endif
         else
          vtr = -0.1021 + 4.932E3*Dr(i) - 0.9551E6*Dr(i)*Dr(i) &
              + 0.07934E9*Dr(i)*Dr(i)*Dr(i) &
              - 0.002362E12*Dr(i)*Dr(i)*Dr(i)*Dr(i)
          stokes = Dc(j)*Dc(j)*vtr*rho_w/(9.*1.718E-5*Dr(i))
          reynolds = 9.*stokes/(p*p*rho_w)
          F = DLOG(reynolds)
          G = -0.1007D0 - 0.358D0*F + 0.0261D0*F*F
          K0 = DEXP(G)
          z = DLOG(stokes/(K0+1.D-15))
          H = 0.1465D0 + 1.302D0*z - 0.607D0*z*z + 0.293D0*z*z*z
          yc0 = 2.0D0/PI * ATAN(H)
          Ef_rw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))
         endif
         t_Efrw(i,j) = MAX(0.0, MIN(SNGL(Ef_rw), 0.95))
      enddo
      enddo
      end subroutine table_Efrw
!
!+---+-----------------------------------------------------------------+
!..Variable collision efficiency for snow collecting cloud water using
!.. method of Wang and Ji, 2000 except equate melted snow diameter to
!.. their "effective collision cross-section."
!+---+-----------------------------------------------------------------+
      subroutine table_Efsw
      implicit none
!..Local variables
      DOUBLE PRECISION:: Ds_m, vts, vtc, stokes, reynolds, Ef_sw
      DOUBLE PRECISION:: p, yc0, F, G, H, z, K0
      INTEGER:: i, j
      do j = 1, nbc
      vtc = 1.19D4 * (1.0D4*Dc(j)*Dc(j)*0.25D0)
      do i = 1, nbs
         vts = av_s*Ds(i)**bv_s * DEXP(-fv_s*Ds(i)) - vtc
         Ds_m = (am_s*Ds(i)**bm_s / am_r)**obmr
         p = Dc(j)/Ds_m
         if (p.gt.0.25 .or. Ds(i).lt.D0s .or. Dc(j).lt.6.E-6 &
               .or. vts.lt.1.E-3) then
          t_Efsw(i,j) = 0.0
         else
          stokes = Dc(j)*Dc(j)*vts*rho_w/(9.*1.718E-5*Ds_m)
          reynolds = 9.*stokes/(p*p*rho_w)
          F = DLOG(reynolds)
          G = -0.1007D0 - 0.358D0*F + 0.0261D0*F*F
          K0 = DEXP(G)
          z = DLOG(stokes/(K0+1.D-15))
          H = 0.1465D0 + 1.302D0*z - 0.607D0*z*z + 0.293D0*z*z*z
          yc0 = 2.0D0/PI * ATAN(H)
          Ef_sw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))
          t_Efsw(i,j) = MAX(0.0, MIN(SNGL(Ef_sw), 0.95))
         endif
      enddo
      enddo
      end subroutine table_Efsw
!
!+---+-----------------------------------------------------------------+
!..Integrate rain size distribution from zero to D-star to compute the
!.. number of drops smaller than D-star that evaporate in a single
!.. timestep.  Drops larger than D-star dont evaporate entirely so do
!.. not affect number concentration.
!+---+-----------------------------------------------------------------+
      subroutine table_dropEvap
      implicit none
!..Local variables
      DOUBLE PRECISION:: Nt_r, N0, lam_exp, lam
      REAL:: xlimit_intg
      INTEGER:: i, j, k
      do k = 1, ntb_r
      do j = 1, ntb_r1
         lam_exp = (N0r_exp(j)*am_r*crg(1)/r_r(k))**ore1
         lam = lam_exp * (crg(3)*org2*org1)**obmr
         N0 = N0r_exp(j)/(crg(2)*lam_exp) * lam**cre(2)
         Nt_r = N0 * crg(2) / lam**cre(2)
         do i = 1, nbr
            xlimit_intg = lam*Dr(i)
            tnr_rev(i,j,k) = GAMMP(mu_r+1.0, xlimit_intg) * Nt_r
         enddo
      enddo
      enddo
      end subroutine table_dropEvap
! TO APPLY TABLE ABOVE
!..Rain lookup table indexes.
!         Dr_star = DSQRT(-2.D0*DT * t1_evap/(2.*PI) &
!                 * 0.78*4.*diffu(k)*xsat*rvs/rho_w)
!         idx_d = NINT(1.0 + FLOAT(nbr) * DLOG(Dr_star/D0r)             &
!               / DLOG(Dr(nbr)/D0r))
!         idx_d = MAX(1, MIN(idx_d, nbr))
!
!         nir = NINT(ALOG10(rr(k)))
!         do nn = nir-1, nir+1
!            n = nn
!            if ( (rr(k)/10.**nn).ge.1.0 .and. &
!                 (rr(k)/10.**nn).lt.10.0) goto 154
!         enddo
!154      continue
!         idx_r = INT(rr(k)/10.**n) + 10*(n-nir2) - (n-nir2)
!         idx_r = MAX(1, MIN(idx_r, ntb_r))
!
!         lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
!         lam_exp = lamr * (crg(3)*org2*org1)**bm_r
!         N0_exp = org1*rr(k)/am_r * lam_exp**cre(1)
!         nir = NINT(DLOG10(N0_exp))
!         do nn = nir-1, nir+1
!            n = nn
!            if ( (N0_exp/10.**nn).ge.1.0 .and. &
!                 (N0_exp/10.**nn).lt.10.0) goto 155
!         enddo
!155      continue
!         idx_r1 = INT(N0_exp/10.**n) + 10*(n-nir3) - (n-nir3)
!         idx_r1 = MAX(1, MIN(idx_r1, ntb_r1))
!
!         pnr_rev(k) = MIN(nr(k)*odts, SNGL(tnr_rev(idx_d,idx_r1,idx_r) &   ! RAIN2M
!                    * odts))
!
!
!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
!     --- RETURNS THE INCOMPLETE GAMMA FUNCTION Q(A,X) EVALUATED BY ITS
!     --- CONTINUED FRACTION REPRESENTATION AS GAMMCF.  ALSO RETURNS
!     --- LN(GAMMA(A)) AS GLN.  THE CONTINUED FRACTION IS EVALUATED BY
!     --- A MODIFIED LENTZ METHOD.
!     --- USES GAMMLN
      IMPLICIT NONE
      INTEGER, PARAMETER:: ITMAX=100
      REAL, PARAMETER:: gEPS=3.E-7
      REAL, PARAMETER:: FPMIN=1.E-30
      REAL, INTENT(IN):: A, X
      REAL:: GAMMCF,GLN
      INTEGER:: I
      REAL:: AN,B,C,D,DEL,H
      GLN=GAMMLN(A)
      B=X+1.-A
      C=1./FPMIN
      D=1./B
      H=D
      DO 11 I=1,ITMAX
        AN=-I*(I-A)
        B=B+2.
        D=AN*D+B
        IF(ABS(D).LT.FPMIN)D=FPMIN
        C=B+AN/C
        IF(ABS(C).LT.FPMIN)C=FPMIN
        D=1./D
        DEL=D*C
        H=H*DEL
        IF(ABS(DEL-1.).LT.gEPS)GOTO 1
 11   CONTINUE
      PRINT *, 'A TOO LARGE, ITMAX TOO SMALL IN GCF'
 1    GAMMCF=EXP(-X+A*LOG(X)-GLN)*H
      END SUBROUTINE GCF
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!+---+-----------------------------------------------------------------+
      SUBROUTINE GSER(GAMSER,A,X,GLN)
!     --- RETURNS THE INCOMPLETE GAMMA FUNCTION P(A,X) EVALUATED BY ITS
!     --- ITS SERIES REPRESENTATION AS GAMSER.  ALSO RETURNS LN(GAMMA(A)) 
!     --- AS GLN.
!     --- USES GAMMLN
      IMPLICIT NONE
      INTEGER, PARAMETER:: ITMAX=100
      REAL, PARAMETER:: gEPS=3.E-7
      REAL, INTENT(IN):: A, X
      REAL:: GAMSER,GLN
      INTEGER:: N
      REAL:: AP,DEL,SUM
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.) PRINT *, 'X < 0 IN GSER'
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*gEPS)GOTO 1
 11   CONTINUE
      PRINT *,'A TOO LARGE, ITMAX TOO SMALL IN GSER'
 1    GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      END SUBROUTINE GSER
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!+---+-----------------------------------------------------------------+
      REAL FUNCTION GAMMLN(XX)
!     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.
      IMPLICIT NONE
      REAL, INTENT(IN):: XX
      DOUBLE PRECISION, PARAMETER:: STP = 2.5066282746310005D0
      DOUBLE PRECISION, DIMENSION(6), PARAMETER:: &
               COF = (/76.18009172947146D0, -86.50532032941677D0, &
                       24.01409824083091D0, -1.231739572450155D0, &
                      .1208650973866179D-2, -.5395239384953D-5/)
      DOUBLE PRECISION:: SER,TMP,X,Y
      INTEGER:: J
      X=XX
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*LOG(TMP)-TMP
      SER=1.000000000190015D0
      DO 11 J=1,6
        Y=Y+1.D0
        SER=SER+COF(J)/Y
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER/X)
      END FUNCTION GAMMLN
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!+---+-----------------------------------------------------------------+
      REAL FUNCTION GAMMP(A,X)
!     --- COMPUTES THE INCOMPLETE GAMMA FUNCTION P(A,X)
!     --- SEE ABRAMOWITZ AND STEGUN 6.5.1
!     --- USES GCF,GSER
      IMPLICIT NONE
      REAL, INTENT(IN):: A,X
      REAL:: GAMMCF,GAMSER,GLN
      GAMMP = 0.
      IF((X.LT.0.) .OR. (A.LE.0.)) THEN
        PRINT *, 'BAD ARGUMENTS IN GAMMP'
        RETURN
      ELSEIF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.-GAMMCF
      ENDIF
      END FUNCTION GAMMP
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!+---+-----------------------------------------------------------------+
      REAL FUNCTION WGAMMA(y)
      IMPLICIT NONE
      REAL, INTENT(IN):: y
      WGAMMA = EXP(GAMMLN(y))
      END FUNCTION WGAMMA
!+---+-----------------------------------------------------------------+
! THIS FUNCTION CALCULATES THE LIQUID SATURATION VAPOR MIXING RATIO AS
! A FUNCTION OF TEMPERATURE AND PRESSURE
!
      REAL FUNCTION RSLF(P,T)
      IMPLICIT NONE
      REAL, INTENT(IN):: P, T
      REAL:: ESL,X
      REAL, PARAMETER:: C0= .611583699E03
      REAL, PARAMETER:: C1= .444606896E02
      REAL, PARAMETER:: C2= .143177157E01
      REAL, PARAMETER:: C3= .264224321E-1
      REAL, PARAMETER:: C4= .299291081E-3
      REAL, PARAMETER:: C5= .203154182E-5
      REAL, PARAMETER:: C6= .702620698E-8
      REAL, PARAMETER:: C7= .379534310E-11
      REAL, PARAMETER:: C8=-.321582393E-13
      X=MAX(-80.,T-273.16)
!      ESL=612.2*EXP(17.67*X/(T-29.65))
      ESL=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
      RSLF=.622*ESL/(P-ESL)
!    ALTERNATIVE
!  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
!             supercooled water for atmospheric applications, Q. J. R.
!             Meteorol. Soc (2005), 131, pp. 1539-1565.
!    ESL = EXP(54.842763 - 6763.22 / T - 4.210 * ALOG(T) + 0.000367 * T
!        + TANH(0.0415 * (T - 218.8)) * (53.878 - 1331.22
!        / T - 9.44523 * ALOG(T) + 0.014025 * T))
      END FUNCTION RSLF
!+---+-----------------------------------------------------------------+
! THIS FUNCTION CALCULATES THE ICE SATURATION VAPOR MIXING RATIO AS A
! FUNCTION OF TEMPERATURE AND PRESSURE
!
      REAL FUNCTION RSIF(P,T)
      IMPLICIT NONE
      REAL, INTENT(IN):: P, T
      REAL:: ESI,X
      REAL, PARAMETER:: C0= .609868993E03
      REAL, PARAMETER:: C1= .499320233E02
      REAL, PARAMETER:: C2= .184672631E01
      REAL, PARAMETER:: C3= .402737184E-1
      REAL, PARAMETER:: C4= .565392987E-3
      REAL, PARAMETER:: C5= .521693933E-5
      REAL, PARAMETER:: C6= .307839583E-7
      REAL, PARAMETER:: C7= .105785160E-9
      REAL, PARAMETER:: C8= .161444444E-12
      X=MAX(-80.,T-273.16)
      ESI=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
      RSIF=.622*ESI/(P-ESI)
!    ALTERNATIVE
!  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
!             supercooled water for atmospheric applications, Q. J. R.
!             Meteorol. Soc (2005), 131, pp. 1539-1565.
!     ESI = EXP(9.550426 - 5723.265/T + 3.53068*ALOG(T) - 0.00728332*T)
      END FUNCTION RSIF
!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
!..Compute _radiation_ effective radii of cloud water, ice, and snow.
!.. These are entirely consistent with microphysics assumptions, not
!.. constant or otherwise ad hoc as is internal to most radiation
!.. schemes.  Since only the smallest snowflakes should impact
!.. radiation, compute from first portion of complicated Field number
!.. distribution, not the second part, which is the larger sizes.
!+---+-----------------------------------------------------------------+
      subroutine calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d,   &
     &                re_qc1d, re_qi1d, re_qs1d, kts, kte)
      IMPLICIT NONE
!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte
      REAL, DIMENSION(kts:kte), INTENT(IN):: &
     &                    t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: re_qc1d, re_qi1d, re_qs1d
!..Local variables
      INTEGER:: k
      REAL, DIMENSION(kts:kte):: rho, rc, nc, ri, ni, rs
      REAL:: smo2, smob, smoc
      REAL:: tc0, loga_, a_, b_
      DOUBLE PRECISION:: lamc, lami
      LOGICAL:: has_qc, has_qi, has_qs
      INTEGER:: inu_c
      real, dimension(15), parameter:: g_ratio = (/24,60,120,210,336,   &
     &                504,720,990,1320,1716,2184,2730,3360,4080,4896/)
      has_qc = .false.
      has_qi = .false.
      has_qs = .false.
      do k = kts, kte
         rho(k) = 0.622*p1d(k)/(R*t1d(k)*(qv1d(k)+0.622))
         rc(k) = MAX(R1, qc1d(k)*rho(k))
         nc(k) = MAX(R2, nc1d(k)*rho(k))
         if (rc(k).gt.R1 .and. nc(k).gt.R2) has_qc = .true.
         ri(k) = MAX(R1, qi1d(k)*rho(k))
         ni(k) = MAX(R2, ni1d(k)*rho(k))
         if (ri(k).gt.R1 .and. ni(k).gt.R2) has_qi = .true.
         rs(k) = MAX(R1, qs1d(k)*rho(k))
         if (rs(k).gt.R1) has_qs = .true.
      enddo
      if (has_qc) then
      do k = kts, kte
         if (rc(k).le.R1 .or. nc(k).le.R2) CYCLE
         inu_c = MIN(15, NINT(1000.E6/nc(k)) + 2)
         lamc = (nc(k)*am_r*g_ratio(inu_c)/rc(k))**obmr
         re_qc1d(k) = MAX(2.51E-6, MIN(SNGL(0.5D0 * DBLE(3.+inu_c)/lamc), 50.E-6))
      enddo
      endif
      if (has_qi) then
      do k = kts, kte
         if (ri(k).le.R1 .or. ni(k).le.R2) CYCLE
         lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
         re_qi1d(k) = MAX(5.01E-6, MIN(SNGL(0.5D0 * DBLE(3.+mu_i)/lami), 125.E-6))
      enddo
      endif
      if (has_qs) then
      do k = kts, kte
         if (rs(k).le.R1) CYCLE
         tc0 = MIN(-0.1, t1d(k)-273.15)
         smob = rs(k)*oams
!..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
!.. then we must compute actual 2nd moment and use as reference.
         if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
            smo2 = smob
         else
            loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
     &         + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
     &         + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
     &         + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*bm_s*bm_s*bm_s
            a_ = 10.0**loga_
            b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
     &         + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
     &         + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
     &         + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
     &         + sb(10)*bm_s*bm_s*bm_s
            smo2 = (smob/a_)**(1./b_)
         endif
!..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
     &         + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
     &         + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
     &         + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*cse(1)*cse(1)*cse(1)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
     &        + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
     &        + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
     &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
         smoc = a_ * smo2**b_
         re_qs1d(k) = MAX(25.E-6, MIN(0.5*smoc/smob, 999.E-6))
      enddo
      endif
      end subroutine calc_effectRad
!+---+-----------------------------------------------------------------+
!..Compute radar reflectivity assuming 10 cm wavelength radar and using
!.. Rayleigh approximation.  Only complication is melted snow/graupel
!.. which we treat as water-coated ice spheres and use Uli Blahak's
!.. library of routines.  The meltwater fraction is simply the amount
!.. of frozen species remaining from what initially existed at the
!.. melting level interface.
!+---+-----------------------------------------------------------------+
      subroutine calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d,     &
                          t1d, p1d, dBZ, kts, kte, ii, jj)
      IMPLICIT NONE
!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte, ii, jj
      REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
                          qv1d, qc1d, qr1d, nr1d, qs1d, qg1d, t1d, p1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ
!     REAL, DIMENSION(kts:kte), INTENT(INOUT):: vt_dBZ
!..Local variables
      REAL, DIMENSION(kts:kte):: temp, pres, qv, rho, rhof
      REAL, DIMENSION(kts:kte):: rc, rr, nr, rs, rg
      DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilamg, N0_r, N0_g
      REAL, DIMENSION(kts:kte):: mvd_r
      REAL, DIMENSION(kts:kte):: smob, smo2, smoc, smoz
      REAL:: oM3, M0, Mrat, slam1, slam2, xDs
      REAL:: ils1, ils2, t1_vts, t2_vts, t3_vts, t4_vts
      REAL:: vtr_dbz_wt, vts_dbz_wt, vtg_dbz_wt
      REAL, DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel
      DOUBLE PRECISION:: N0_exp, N0_min, lam_exp, lamr, lamg
      REAL:: a_, b_, loga_, tc0
      DOUBLE PRECISION:: fmelt_s, fmelt_g
      INTEGER:: i, k, k_0, kbot, n
      LOGICAL:: melti
      LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg
      DOUBLE PRECISION:: cback, x, eta, f_d
      REAL:: xslw1, ygra1, zans1
!+---+
      do k = kts, kte
         dBZ(k) = -35.0
      enddo
!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         rhof(k) = SQRT(RHO_NOT/rho(k))
         rc(k) = MAX(R1, qc1d(k)*rho(k))
         if (qr1d(k) .gt. R1) then
            rr(k) = qr1d(k)*rho(k)
            nr(k) = MAX(R2, nr1d(k)*rho(k))
            lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
            ilamr(k) = 1./lamr
            N0_r(k) = nr(k)*org2*lamr**cre(2)
            mvd_r(k) = (3.0 + mu_r + 0.672) * ilamr(k)
            L_qr(k) = .true.
         else
            rr(k) = R1
            nr(k) = R1
            mvd_r(k) = 50.E-6
            L_qr(k) = .false.
         endif
         if (qs1d(k) .gt. R2) then
            rs(k) = qs1d(k)*rho(k)
            L_qs(k) = .true.
         else
            rs(k) = R1
            L_qs(k) = .false.
         endif
         if (qg1d(k) .gt. R2) then
            rg(k) = qg1d(k)*rho(k)
            L_qg(k) = .true.
         else
            rg(k) = R1
            L_qg(k) = .false.
         endif
      enddo
!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope, and useful moments for snow.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         tc0 = MIN(-0.1, temp(k)-273.15)
         smob(k) = rs(k)*oams
!..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
!.. then we must compute actual 2nd moment and use as reference.
         if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
            smo2(k) = smob(k)
         else
            loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
               + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
               + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
               + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
               + sa(10)*bm_s*bm_s*bm_s
            a_ = 10.0**loga_
            b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
               + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
               + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
               + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
               + sb(10)*bm_s*bm_s*bm_s
            smo2(k) = (smob(k)/a_)**(1./b_)
         endif
!..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
               + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
               + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
               + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(1)*cse(1)*cse(1)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
              + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
              + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
         smoc(k) = a_ * smo2(k)**b_
!..Calculate bm_s*2 (th) moment.  Useful for reflectivity.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(3) &
               + sa(4)*tc0*cse(3) + sa(5)*tc0*tc0 &
               + sa(6)*cse(3)*cse(3) + sa(7)*tc0*tc0*cse(3) &
               + sa(8)*tc0*cse(3)*cse(3) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(3)*cse(3)*cse(3)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(3) + sb(4)*tc0*cse(3) &
              + sb(5)*tc0*tc0 + sb(6)*cse(3)*cse(3) &
              + sb(7)*tc0*tc0*cse(3) + sb(8)*tc0*cse(3)*cse(3) &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(3)*cse(3)*cse(3)
         smoz(k) = a_ * smo2(k)**b_
      enddo
!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for graupel.
!+---+-----------------------------------------------------------------+
      N0_min = gonv_max
      do k = kte, kts, -1
         if (temp(k).lt.270.65 .and. L_qr(k) .and. mvd_r(k).gt.100.E-6) then
            xslw1 = 4.01 + alog10(mvd_r(k))
         else
            xslw1 = 0.01
         endif
         ygra1 = 4.31 + alog10(max(5.E-5, rg(k)))
         zans1 = 3.1 + (100./(300.*xslw1*ygra1/(10./xslw1+1.+0.25*ygra1)+30.+10.*ygra1))
         N0_exp = 10.**(zans1)
         N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
         N0_min = MIN(N0_exp, N0_min)
         N0_exp = N0_min
         lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
         lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
         ilamg(k) = 1./lamg
         N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
      enddo
!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
      melti = .false.
      k_0 = kts
      do k = kte-1, kts, -1
         if ( (temp(k).gt.273.15) .and. L_qr(k)                         &
                                  .and. (L_qs(k+1).or.L_qg(k+1)) ) then
            k_0 = MAX(k+1, k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue
!+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*crg(4)*ilamr(k)**cre(4)
         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
                                 * (am_s/900.0)*(am_s/900.0)*smoz(k)
         if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
                                    * (am_g/900.0)*(am_g/900.0)         &
                                    * N0_g(k)*cgg(4)*ilamg(k)**cge(4)
      enddo
!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+
      if (.not. iiwarm .and. melti .and. k_0.ge.2) then
       do k = k_0-1, kts, -1
!..Reflectivity contributed by melting snow
          if (L_qs(k) .and. L_qs(k_0) ) then
           fmelt_s = MAX(0.05d0, MIN(1.0d0-rs(k)/rs(k_0), 0.99d0))
           eta = 0.d0
           oM3 = 1./smoc(k)
           M0 = (smob(k)*oM3)
           Mrat = smob(k)*M0*M0*M0
           slam1 = M0 * Lam0
           slam2 = M0 * Lam1
           do n = 1, nrbins
              x = am_s * xxDs(n)**bm_s
              call rayleigh_soak_wetgraupel (x, DBLE(ocms), DBLE(obms), &
                    fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_s, matrixstring_s,          &
                    inclusionstring_s, hoststring_s,                    &
                    hostmatrixstring_s, hostinclusionstring_s)
              f_d = Mrat*(Kap0*DEXP(-slam1*xxDs(n))                     &
                    + Kap1*(M0*xxDs(n))**mu_s * DEXP(-slam2*xxDs(n)))
              eta = eta + f_d * CBACK * simpson(n) * xdts(n)
           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif
!..Reflectivity contributed by melting graupel
          if (L_qg(k) .and. L_qg(k_0) ) then
           fmelt_g = MAX(0.05d0, MIN(1.0d0-rg(k)/rg(k_0), 0.99d0))
           eta = 0.d0
           lamg = 1./ilamg(k)
           do n = 1, nrbins
              x = am_g * xxDg(n)**bm_g
              call rayleigh_soak_wetgraupel (x, DBLE(ocmg), DBLE(obmg), &
                    fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_g, matrixstring_g,          &
                    inclusionstring_g, hoststring_g,                    &
                    hostmatrixstring_g, hostinclusionstring_g)
              f_d = N0_g(k)*xxDg(n)**mu_g * DEXP(-lamg*xxDg(n))
              eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
           enddo
           ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif
       enddo
      endif
      do k = kte, kts, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo
!..Reflectivity-weighted terminal velocity (snow, rain, graupel, mix).
!     do k = kte, kts, -1
!        vt_dBZ(k) = 1.E-3
!        if (rs(k).gt.R2) then
!         Mrat = smob(k) / smoc(k)
!         ils1 = 1./(Mrat*Lam0 + fv_s)
!         ils2 = 1./(Mrat*Lam1 + fv_s)
!         t1_vts = Kap0*csg(5)*ils1**cse(5)
!         t2_vts = Kap1*Mrat**mu_s*csg(11)*ils2**cse(11)
!         ils1 = 1./(Mrat*Lam0)
!         ils2 = 1./(Mrat*Lam1)
!         t3_vts = Kap0*csg(6)*ils1**cse(6)
!         t4_vts = Kap1*Mrat**mu_s*csg(12)*ils2**cse(12)
!         vts_dbz_wt = rhof(k)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)
!         if (temp(k).ge.273.15 .and. temp(k).lt.275.15) then
!            vts_dbz_wt = vts_dbz_wt*1.5
!         elseif (temp(k).ge.275.15) then
!            vts_dbz_wt = vts_dbz_wt*2.0
!         endif
!        else
!         vts_dbz_wt = 1.E-3
!        endif
!        if (rr(k).gt.R1) then
!         lamr = 1./ilamr(k)
!         vtr_dbz_wt = rhof(k)*av_r*crg(13)*(lamr+fv_r)**(-cre(13))      &
!                    / (crg(4)*lamr**(-cre(4)))
!        else
!         vtr_dbz_wt = 1.E-3
!        endif
!        if (rg(k).gt.R2) then
!         lamg = 1./ilamg(k)
!         vtg_dbz_wt = rhof(k)*av_g*cgg(5)*lamg**(-cge(5))               &
!                    / (cgg(4)*lamg**(-cge(4)))
!        else
!         vtg_dbz_wt = 1.E-3
!        endif
!        vt_dBZ(k) = (vts_dbz_wt*ze_snow(k) + vtr_dbz_wt*ze_rain(k)      &
!                     + vtg_dbz_wt*ze_graupel(k))                        &
!                     / (ze_rain(k)+ze_snow(k)+ze_graupel(k))
!     enddo
      end subroutine calc_refl10cm
!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
END MODULE module_mp_thompiso
!+---+-----------------------------------------------------------------+
