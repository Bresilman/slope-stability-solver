module data_structure
        implicit none
        
        integer, parameter :: DP = kind(1.0D0)
        
        ! --- UNIFIED REINFORCEMENT STRUCTURE ---
        type :: reinforcement_t
            integer  :: type_id       ! 1 = Geogrid/Geotextile, 2 = Soil Nail/Anchor
            integer  :: id            ! Identification number
            
            ! --- Geometry ---
            real(DP) :: x_start       ! X coordinate of the "Head" (Slope face)
            real(DP) :: y_start       ! Y coordinate of the "Head"
            real(DP) :: x_end         ! X coordinate of the "Tail" (Anchorage)
            real(DP) :: y_end         ! Y coordinate of the "Tail"
            
            ! --- Strength Properties ---
            real(DP) :: tensile_str   ! Tensile Strength (kN/m or kN/nail)
            real(DP) :: shear_str     ! Shear Strength (kN/nail)
            real(DP) :: spacing       ! Horizontal spacing (m)
            
            ! --- Pullout Properties ---
            real(DP) :: bond_str      ! Cohesive Bond Strength (kN/m). If > 0, overrides friction.
            real(DP) :: friction_coeff! Soil-Reinforcement Friction Coefficient (tan(delta)). Default ~0.8*tan(phi).
            
            ! --- Derived/Helper ---
            real(DP) :: orientation_rad ! Inclination angle from horizontal
        end type reinforcement_t

        ! --- SURCHARGE STRUCTURE ---
        type :: surcharge_t
            real(DP) :: x_start
            real(DP) :: x_end
            real(DP) :: magnitude
        end type surcharge_t

        type :: slice_t
            real(DP) :: b, W, base_ang_rad, base_length
            real(DP) :: pore_pressure, cohesion, friction_ang_rad
            real(DP) :: x_mid, y_base, height
            integer  :: material_id
        end type slice_t    
        
        type :: soil_t
            real(DP) :: gamma, cohesion, friction_ang_rad 
        end type soil_t
        
    end module data_structure
