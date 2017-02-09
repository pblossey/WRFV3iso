MODULE module_mp_water_isotopologues

  !bloss: Constants for computing isotope fractionation
  ! NOTE iso_index=18, H2O16
  ! NOTE iso_index=19, HDO16
  ! NOTE iso_index=20, H2O18
  ! NOTE iso_index=17, H2O17 (!!NOT IMPLEMENTED YET!!)
  !
  !                 DEUTERIUM         O18
  REAL, DIMENSION(18:20), PARAMETER, PRIVATE :: &
       calI1 = (/  0.,   16288.,      0.       /), & ! coefs in \alpha_I computation
       calI2 = (/  0.,   0.,          11.839    /), &
       calI3=  (/  0.,  -0.0934,     -0.028224 /), &
       calW1 = (/  0.,   24884.,      1137.       /), & ! coefs in \alpha_W computation
       calW2=  (/  0.,  -76.248,     -0.4156   /), &
       calW3=  (/  0.,   0.052612,   -0.0020667/)
  REAL, DIMENSION(18:20), PARAMETER :: & ! Merlivat (1978), ratio of molecular diffusivities
       Drat_light_over_heavy =  (/     1.0000,  1./0.9755,     1./0.9723 /), &  ! D_v / D_v'
       Mrat_light_over_heavy =  (/     18./18.,    18./19.,      18./20. /)    ! M/M' where ' is the heavy isotope

  INTEGER, DIMENSION(1:3), PARAMETER :: &
       iso_index_ref = (/  18,  19,  20 /)         ! indices needed for alpha_* computations

  CHARACTER(LEN=3), DIMENSION(18:20), PARAMETER :: &
       iso_string = (/ 'H2O', 'HDO', 'O18' /)

  !iso_factor is a small number to make the contribution of isotopes to
  !  the mass negligible It is an inconvenience due to putting isotopes in
  !  microphysics.
  REAL     , PARAMETER ::            iso_factor=1.e-10
  !
CONTAINS 

  !****************************************************************
  !****************************************************************
  !bloss 061209
  !Add function that only computes equilibrium alpha for liquid.
  !  Coefficients are from Majoube (1971, Journal
  !  de Chimie Physique et Physico-Chimie Biologique, vol. 68,
  !  no. 10, pp. 1423-1436).
  REAL FUNCTION alfaW_equilibrium(Tk,iso_index)
    !----------------------------------------------------------------
    IMPLICIT NONE
    !----------------------------------------------------------------
    REAL, INTENT(IN   ) :: Tk
    INTEGER, INTENT(IN) :: iso_index !==19 for HDO, ==20 for H2O18

    alfaW_equilibrium = &
         exp( (calW1(iso_index) + Tk*(calW2(iso_index) + Tk*calW3(iso_index))) &
              /(Tk*Tk) )

  END FUNCTION alfaW_equilibrium

  !****************************************************************
  !****************************************************************
  !bloss 061209
  !Add function that only computes equilibrium alpha for ice.
  !  Coefficients for O18 are from Majoube (1970, Nature, vol. 226,
  !  p. 1242, 27 June 1970, doi: 10.1038/2261242a0).  Those for
  !  deuterium are from Merlivat & Nief (1967, Tellus, doi:
  !  10.1111/j.2153-3490.1967.tb01465.x) as presented in Jouzel (1986,
  !  Handbook of Environmental Isotope Geochemistry) with a couple
  !  of differences in the last decimal places.
  !   Merlivat & Nief: alfaI = exp(-0.0945 + 16289/T^2)
  !   Jouzel 1986: alfaI = exp(-0.0934 + 16288/T^2)

  REAL FUNCTION alfaI_equilibrium(Tk,iso_index)
    !----------------------------------------------------------------
    IMPLICIT NONE
    !----------------------------------------------------------------
    REAL, INTENT(IN   ) :: Tk
    INTEGER, INTENT(IN) :: iso_index !==19 for HDO, ==20 for H2O18

    alfaI_equilibrium = &
         exp( (calI1(iso_index) + Tk*(calI2(iso_index) + Tk*calI3(iso_index))) &
              /(Tk*Tk) )

  END FUNCTION alfaI_equilibrium

  !****************************************************************
  !****************************************************************
  REAL FUNCTION alfaK(alpha_equil, S_tilde, Vrat_light_over_heavy, iso_index)
    !bloss (061509, modified 2014-03-31): Note that alpha_k formula 
    !  is identical for liquid and ice.  Use a single function for both.
    !  The difference is in the equilibrium fractionation coefficient
    !  and the value of the saturation ratio S_tilde that are input.
    !  These will reflect whether the exchange occurs between vapor 
    !  and ice or vapor and liquid.
    !  Modification (in 2014): change formulation to remove dependence
    !  on material properties through the constant b in Gedzelman
    !  & Arnold (1994, doi: 10.1029/93JD03518) and the appendix of 
    !  Blossey et al (2010, doi: 10.1029/2010JD014554, eqn B26).  
    !  The current form is that of equation 11 in Jouzel & Merlivat 
    !  (1984, doi: 10.1029/JD089iD07p11749) where our S_tilde 
    !  corresponds to their Si, which is defined in their equation 13.
    !  In addition, we have included the ratio of the ventilation
    !  factors.
    !----------------------------------------------------------------
    IMPLICIT NONE
    !----------------------------------------------------------------
    REAL, INTENT(IN   ) :: alpha_equil  ! equilibrium fractionation
    REAL, INTENT(IN   ) :: S_tilde  ! saturation ratio = rv*^drop/rv_ambient
                                    !  equivalent to Si (Jouzel & Merlivat, 1984, eqn 13).
    REAL, INTENT(IN   ) :: Vrat_light_over_heavy ! = V/V' where V and V' are
                                ! ventilation coefficients for the 
                                ! standard and heavy isotopes, respectively. 
    INTEGER, INTENT(IN) :: iso_index !==19 for HDO, ==20 for H2O18

    alfaK = S_tilde / (1. + alpha_equil*(S_tilde - 1.) &
                             *Drat_light_over_heavy(iso_index) &
                             *Vrat_light_over_heavy)

  END FUNCTION alfaK

END MODULE module_mp_water_isotopologues
