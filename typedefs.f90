MODULE typedefs
  ! Declare type for a node of the binary tree.
  TYPE :: node
    DOUBLE PRECISION          :: tau
    DOUBLE PRECISION          :: omega
    INTEGER                   :: primary_species
    INTEGER                   :: secondary_species
    INTEGER                   :: ix,iy,iz !integer coords
    INTEGER                   :: hop_dir ! The direction of the next hop
    INTEGER                   :: f_i ! Number of hops from site in direction of collecting electrode
    INTEGER                   :: b_i ! Number of hops from site in direction of injecting electrode
    DOUBLE PRECISION          :: erand ! E_rand + E_applied
    DOUBLE PRECISION          :: eapplied
    DOUBLE PRECISION          :: eself
    DOUBLE PRECISION          :: esr
    DOUBLE PRECISION          :: edisk
    DOUBLE PRECISION          :: elayer
    DOUBLE PRECISION          :: etotal ! E_fixed + E_variable
    DOUBLE PRECISION          :: J_i ! The current density of the site
    TYPE (node)     , POINTER :: before
    TYPE (node)     , POINTER :: after
    TYPE (node)     , POINTER :: parent
    INTEGER                   :: leftRight
    ! The variables below are for spin and MD calculations
    INTEGER                   :: ic ! Unique identifyer
    DOUBLE PRECISION          :: x,y,z !x=coord1,y=coord2,z=coord3
    DOUBLE PRECISION          :: vx,vy,vz !velocities
    DOUBLE PRECISION          :: fx,fy,fz !forces
    DOUBLE PRECISION          :: sx,sy,sz !spins
    DOUBLE PRECISION          :: bx,by,bz !fields
    DOUBLE PRECISION          :: socx,socy,socz !spin-orbit-coupling factors
    DOUBLE PRECISION          :: wsx,wsy,wsz
    DOUBLE PRECISION          :: wkx,wky,wkz
END TYPE node

CONTAINS
  SUBROUTINE init_node(temp)
    IMPLICIT NONE
    TYPE(node), POINTER :: temp
    temp%ic                = 0
    temp%tau               = 0.0
    temp%omega             = 0.0
    temp%primary_species   = 0
    temp%secondary_species = 0
    temp%ix                = 0
    temp%iy                = 0
    temp%iz                = 0
    temp%x                 = 0.0
    temp%y                 = 0.0
    temp%z                 = 0.0
    temp%hop_dir           = 0
    temp%erand             = 0.0
    temp%eapplied          = 0.0
    temp%eself             = 0.0
    temp%esr               = 0.0
    temp%edisk             = 0.0
    temp%elayer            = 0.0
    temp%etotal           = 0.0
    temp%J_i               = 0.0
    temp%f_i               = 0
    temp%b_i               = 0
    temp%leftRight         = -1
    temp%vx                = 0.0
    temp%vy                = 0.0
    temp%vz                = 0.0
    temp%fx                = 0.0
    temp%fy                = 0.0
    temp%fz                = 0.0
    temp%bx                = 0.0
    temp%by                = 0.0
    temp%bz                = 0.0
    temp%sx                = 0.0
    temp%sy                = 0.0
    temp%sz                = 0.0
    temp%socx              = 0.0
    temp%socy              = 0.0
    temp%socz              = 0.0
    temp%wsx               = 0.0
    temp%wsy               = 0.0
    temp%wsz               = 0.0
    temp%wkx               = 0.0
    temp%wky               = 0.0
    temp%wkz               = 0.0
    NULLIFY(temp%before, temp%after, temp%parent)
  END SUBROUTINE init_node

  SUBROUTINE wipe_node(temp)
    IMPLICIT NONE
    TYPE(node), POINTER :: temp
    temp%tau               = 0.0
    temp%omega             = 0.0
    temp%secondary_species = 0
    temp%hop_dir           = 0
    temp%eself             = 0.0
    temp%etotal            = temp%erand + temp%eapplied + temp%esr + temp%elayer - temp%edisk
    temp%sx                = 0.0
    temp%sy                = 0.0
    temp%sz                = 0.0
    temp%socx              = 0.0
    temp%socy              = 0.0
    temp%socz              = 0.0
    temp%wsx               = 0.0
    temp%wsy               = 0.0
    temp%wsz               = 0.0
    temp%wkx               = 0.0
    temp%wky               = 0.0
    temp%wkz               = 0.0
  END SUBROUTINE wipe_node
END MODULE typedefs
