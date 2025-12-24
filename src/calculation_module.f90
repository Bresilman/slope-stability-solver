module calculation_module
    use data_structure
    implicit none
    
contains
    
    function calculate_fos_fellenius(slices, reinforcements, xc, yc, r) result(fos)
        type(slice_t), intent(in), dimension(:) :: slices
        type(reinforcement_t), intent(in), dimension(:), optional :: reinforcements
        real(DP), intent(in) :: xc, yc, r
        
        real(DP) :: fos
        real(DP) :: sum_numerator, sum_denominator      
        real(DP) :: driving_force, resisting_force, effective_normal_force
        integer :: i
        real(DP) :: term_reinforcement
        
        sum_numerator = 0.0_DP
        sum_denominator = 0.0_DP
        term_reinforcement = 0.0_DP
        
        ! 1. Soil Forces
        do i = 1, size(slices)
            driving_force = slices(i)%W * sin(slices(i)%base_ang_rad)
            effective_normal_force = (slices(i)%W * cos(slices(i)%base_ang_rad)) - &
                                     (slices(i)%pore_pressure * slices(i)%base_length)
            if (effective_normal_force < 0.0_DP) effective_normal_force = 0.0_DP
            
            resisting_force = slices(i)%cohesion * slices(i)%base_length + &
                              (effective_normal_force * tan(slices(i)%friction_ang_rad))
            
            sum_denominator = sum_denominator + driving_force
            sum_numerator = sum_numerator + resisting_force
        end do
        
        ! 2. Reinforcement Forces
        if (present(reinforcements)) then
            call calculate_reinforcement_contribution(slices, reinforcements, xc, yc, r, term_reinforcement)
            sum_numerator = sum_numerator + term_reinforcement
        end if
        
        if (abs(sum_denominator) < 1.0E-9_DP) then
            fos = 999.0_DP 
        else
            fos = sum_numerator / abs(sum_denominator)
        end if
    end function calculate_fos_fellenius

    ! --------------------------------------------------------------------------
    ! Bishop Simplified Method (Iterative) - ROBUST VERSION
    ! --------------------------------------------------------------------------
    function calculate_fos_bishop(slices, reinforcements, xc, yc, r, initial_guess) result(fos)
        type(slice_t), intent(in), dimension(:) :: slices
        type(reinforcement_t), intent(in), dimension(:), optional :: reinforcements
        real(DP), intent(in) :: xc, yc, r
        real(DP), intent(in), optional :: initial_guess
        
        real(DP) :: fos, fos_trial, fos_new, error_diff
        real(DP) :: sum_driving, sum_resisting
        real(DP) :: m_alpha, term_soil_strength, term_reinforcement
        real(DP) :: sin_alpha, cos_alpha, tan_alpha, tan_phi
        integer :: i, iter
        integer, parameter :: MAX_ITER = 50
        real(DP), parameter :: TOLERANCE = 0.005_DP
        
        ! Initial Guess
        if (present(initial_guess)) then
            fos_trial = initial_guess
        else
            fos_trial = 1.0_DP
        end if
        
        ! Ensure initial guess is sane
        if (fos_trial < 0.1_DP) fos_trial = 1.0_DP
        
        ! Pre-calculate Reinforcement (Constant in Bishop Simplified)
        term_reinforcement = 0.0_DP
        if (present(reinforcements)) then
            call calculate_reinforcement_contribution(slices, reinforcements, xc, yc, r, term_reinforcement)
        end if
        
        ! Iteration Loop
        fos = 999.0_DP ! Default fail state
        
        do iter = 1, MAX_ITER
            sum_driving = 0.0_DP
            sum_resisting = 0.0_DP
            
            do i = 1, size(slices)
                sin_alpha = sin(slices(i)%base_ang_rad)
                cos_alpha = cos(slices(i)%base_ang_rad)
                tan_alpha = tan(slices(i)%base_ang_rad) ! Calculated here
                tan_phi   = tan(slices(i)%friction_ang_rad)
                
                ! Denominator: Sum of W * sin(alpha)
                sum_driving = sum_driving + (slices(i)%W * sin_alpha)
                
                ! Numerator Term: m_alpha Calculation
                ! m_alpha = cos(alpha) * (1 + tan(alpha)*tan(phi)/F)
                
                if (fos_trial < 0.001_DP) fos_trial = 0.001_DP
                
                ! Use tan_alpha variable here
                m_alpha = cos_alpha * (1.0_DP + (tan_alpha * tan_phi) / fos_trial)
                
                ! Stability: Prevent m_alpha from being too small
                if (m_alpha < 0.2_DP) m_alpha = 0.2_DP 
                
                term_soil_strength = (slices(i)%cohesion * slices(i)%b) + &
                                     ((slices(i)%W - slices(i)%pore_pressure * slices(i)%b) * tan_phi)
                
                sum_resisting = sum_resisting + (term_soil_strength / m_alpha)
            end do
            
            ! Add reinforcement to resisting moment
            sum_resisting = sum_resisting + term_reinforcement
            
            if (abs(sum_driving) < 1.0E-9_DP) then
                fos = 999.0_DP 
                return
            end if
            
            fos_new = sum_resisting / sum_driving
            
            if (fos_new < 0.0_DP) then
                 fos = 999.0_DP
                 return 
            end if
            
            error_diff = abs(fos_new - fos_trial)
            if (error_diff < TOLERANCE) then
                fos = fos_new
                exit
            end if
            
            fos_trial = 0.5_DP * fos_trial + 0.5_DP * fos_new
        end do
        
        fos = fos_trial 
        
    end function calculate_fos_bishop

    subroutine calculate_reinforcement_contribution(slices, reinforcements, xc, yc, r, total_force)
        type(slice_t), intent(in), dimension(:) :: slices
        type(reinforcement_t), intent(in), dimension(:) :: reinforcements
        real(DP), intent(in) :: xc, yc, r
        real(DP), intent(inout) :: total_force
        
        integer :: j
        real(DP) :: reinf_y, dy, dx, t1, t2, x_int, y_int
        integer :: num_intersections
        
        do j = 1, size(reinforcements)
            if (reinforcements(j)%type_id == 1) then
                reinf_y = reinforcements(j)%y_start
                if (reinf_y > (yc - r) .and. reinf_y < (yc + r)) then
                    dy = reinf_y - yc
                    dx = sqrt(r**2 - dy**2)
                    call apply_geogrid_force(xc + dx, reinf_y, reinforcements(j), slices, xc, yc, total_force)
                    call apply_geogrid_force(xc - dx, reinf_y, reinforcements(j), slices, xc, yc, total_force)
                end if
            else if (reinforcements(j)%type_id == 2) then
                call intersect_line_circle(reinforcements(j)%x_start, reinforcements(j)%y_start, &
                                           reinforcements(j)%x_end, reinforcements(j)%y_end, &
                                           xc, yc, r, t1, t2, num_intersections)
                if (num_intersections >= 1) then
                    x_int = reinforcements(j)%x_start + t1 * (reinforcements(j)%x_end - reinforcements(j)%x_start)
                    y_int = reinforcements(j)%y_start + t1 * (reinforcements(j)%y_end - reinforcements(j)%y_start)
                    call apply_nail_force(x_int, y_int, reinforcements(j), slices, xc, yc, total_force)
                end if
                if (num_intersections == 2) then
                    x_int = reinforcements(j)%x_start + t2 * (reinforcements(j)%x_end - reinforcements(j)%x_start)
                    y_int = reinforcements(j)%y_start + t2 * (reinforcements(j)%y_end - reinforcements(j)%y_start)
                    call apply_nail_force(x_int, y_int, reinforcements(j), slices, xc, yc, total_force)
                end if
            end if
        end do
    end subroutine calculate_reinforcement_contribution

    subroutine apply_geogrid_force(x_int, y_int, reinf, slices, xc, yc, sum_resisting)
        real(DP), intent(in) :: x_int, y_int
        type(reinforcement_t), intent(in) :: reinf
        type(slice_t), intent(in), dimension(:) :: slices
        real(DP), intent(in) :: xc, yc
        real(DP), intent(inout) :: sum_resisting
        
        real(DP) :: alpha_base, force_tangential
        real(DP) :: l_eff, sigma_v, t_pullout, t_available, fric_coeff
        real(DP) :: y_surface, depth_z
        integer :: k
        
        if (x_int >= reinf%x_start .and. x_int <= reinf%x_end) then
            if (abs(y_int - yc) > 1.0E-9_DP) then
                alpha_base = atan( -1.0_DP * (x_int - xc) / (y_int - yc) )
            else
                alpha_base = 1.570796_DP 
            end if
            
            if (reinf%x_end > reinf%x_start) then
                l_eff = reinf%x_end - x_int
            else
                l_eff = reinf%x_start - x_int
            end if
            if (l_eff < 0.0_DP) l_eff = 0.0_DP
            
            y_surface = y_int 
            do k = 1, size(slices)
                if (x_int >= (slices(k)%x_mid - slices(k)%b/2.0_DP) .and. &
                    x_int <= (slices(k)%x_mid + slices(k)%b/2.0_DP)) then
                    y_surface = slices(k)%y_base + slices(k)%height
                    exit
                end if
            end do
            
            depth_z = y_surface - y_int
            if (depth_z < 0.0_DP) depth_z = 0.0_DP
            sigma_v = 18.0_DP * depth_z
            
            if (reinf%bond_str > 0.0_DP) then
                t_pullout = l_eff * reinf%bond_str
            else
                if (reinf%friction_coeff > 0.0_DP) then
                    fric_coeff = reinf%friction_coeff
                else
                    fric_coeff = 0.46_DP 
                end if
                t_pullout = 2.0_DP * l_eff * sigma_v * fric_coeff
            end if
            
            t_available = min(reinf%tensile_str, t_pullout)
            force_tangential = t_available * cos(alpha_base)
            if (force_tangential > 0.0_DP) sum_resisting = sum_resisting + force_tangential
        end if
    end subroutine apply_geogrid_force

    subroutine apply_nail_force(x_int, y_int, reinf, slices, xc, yc, sum_resisting)
        real(DP), intent(in) :: x_int, y_int
        type(reinforcement_t), intent(in) :: reinf
        type(slice_t), intent(in), dimension(:) :: slices
        real(DP), intent(in) :: xc, yc
        real(DP), intent(inout) :: sum_resisting
        
        real(DP) :: alpha_base, beta_angle
        real(DP) :: t_normalized, s_normalized, t_available, t_pullout, fric_coeff
        real(DP) :: force_contribution, l_eff_stable, y_surface, depth_z, sigma_v
        integer :: k
        
        if (abs(y_int - yc) > 1.0E-9_DP) then
            alpha_base = atan( -1.0_DP * (x_int - xc) / (y_int - yc) )
        else
            alpha_base = 1.570796_DP 
        end if
        
        l_eff_stable = sqrt((reinf%x_end - x_int)**2 + (reinf%y_end - y_int)**2)
        
        if (reinf%bond_str > 0.0_DP) then
            t_pullout = l_eff_stable * reinf%bond_str
        else
            y_surface = y_int 
            do k = 1, size(slices)
                if (x_int >= (slices(k)%x_mid - slices(k)%b/2.0_DP) .and. &
                    x_int <= (slices(k)%x_mid + slices(k)%b/2.0_DP)) then
                    y_surface = slices(k)%y_base + slices(k)%height
                    exit
                end if
            end do
            
            depth_z = y_surface - y_int
            if (depth_z < 0.0_DP) depth_z = 0.0_DP
            sigma_v = 18.0_DP * depth_z
            
            if (reinf%friction_coeff > 0.0_DP) then
                fric_coeff = reinf%friction_coeff
            else
                fric_coeff = 0.577_DP 
            end if
            
            t_pullout = l_eff_stable * (3.14159_DP * 0.15_DP) * (sigma_v * fric_coeff) 
        end if
        
        t_available = min(reinf%tensile_str, t_pullout)
        
        if (reinf%spacing > 0.0_DP) then
            t_normalized = t_available / reinf%spacing
            s_normalized = reinf%shear_str / reinf%spacing
        else
            t_normalized = 0.0_DP; s_normalized = 0.0_DP
        end if
        
        beta_angle = alpha_base - reinf%orientation_rad
        force_contribution = t_normalized * cos(beta_angle) + s_normalized * sin(beta_angle)
        
        if (force_contribution > 0.0_DP) then
            sum_resisting = sum_resisting + force_contribution
        end if
        
    end subroutine apply_nail_force

    subroutine intersect_line_circle(x1, y1, x2, y2, xc, yc, r, t1, t2, n)
        real(DP), intent(in) :: x1, y1, x2, y2, xc, yc, r
        real(DP), intent(out) :: t1, t2
        integer, intent(out) :: n 
        real(DP) :: dx, dy, fx, fy, a, b, c, delta
        
        n = 0; t1 = -1.0_DP; t2 = -1.0_DP
        dx = x2 - x1; dy = y2 - y1
        fx = x1 - xc; fy = y1 - yc
        a = dx*dx + dy*dy
        b = 2.0_DP * (fx*dx + fy*dy)
        c = (fx*fx + fy*fy) - r*r
        delta = b*b - 4.0_DP*a*c
        
        if (delta >= 0.0_DP) then
            delta = sqrt(delta)
            t1 = (-b - delta) / (2.0_DP * a)
            t2 = (-b + delta) / (2.0_DP * a)
            if (t1 >= 0.0_DP .and. t1 <= 1.0_DP) n = n + 1
            if (t2 >= 0.0_DP .and. t2 <= 1.0_DP) n = n + 1
        end if
    end subroutine intersect_line_circle

end module calculation_module
