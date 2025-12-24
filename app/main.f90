program slope_solver
    use data_structure
    use io_module
    use calculation_module
    use slicing_module
    implicit none
    
    type(input_parameters_t) :: params
    type(slice_t), allocatable, dimension(:) :: slices
    logical :: is_valid
    
    real(DP) :: xc, yc, r, fos, fos_fellenius
    real(DP) :: min_fos
    real(DP) :: crit_xc, crit_yc, crit_r
    real(DP) :: gamma_water 
    
    real(DP) :: dx, dy, dr
    integer :: ix, iy, ir
    integer :: valid_circles_count
    integer :: log_unit, err_code
    character(len=20) :: method_str ! Temp string for comparison
    
    print *, "================================================="
    print *, "          SLOPE STABILITY SOLVER V1.7            "
    print *, "          (Fellenius & Bishop / Surcharges)      "
    print *, "================================================="

    call read_input_file("input.txt", params)
    
    ! --- DETERMINE UNIT WEIGHT OF WATER ---
    if (trim(adjustl(params%units)) == 'IMPERIAL') then
        gamma_water = 62.4_DP
        print *, "Units: IMPERIAL (Gamma_w = 62.4 pcf)"
    else
        gamma_water = 9.81_DP
        print *, "Units: METRIC (Gamma_w = 9.81 kN/m3)"
    end if
    
    min_fos = 999.0_DP
    crit_xc = 0.0_DP; crit_yc = 0.0_DP; crit_r = 0.0_DP
    valid_circles_count = 0
    
    ! Calculate steps
    if (params%grid_steps_x > 0) dx = (params%grid_bottom_right(1) - params%grid_top_left(1)) / real(params%grid_steps_x, DP)
    if (params%grid_steps_y > 0) dy = (params%grid_bottom_right(2) - params%grid_top_left(2)) / real(params%grid_steps_y, DP)
    if (params%radius_steps > 0) dr = (params%radius_max - params%radius_min) / real(params%radius_steps, DP)

    call open_log_file("log_results.txt", log_unit)

    ! Clean the method string for comparison
    method_str = trim(adjustl(params%solver_method))
    print *, "Starting Grid Search..."
    print *, "Method: ", method_str

    ! --- SEARCH LOOP ---
    do ix = 0, params%grid_steps_x
        xc = params%grid_top_left(1) + real(ix, DP) * dx
        
        do iy = 0, params%grid_steps_y
            yc = params%grid_top_left(2) + real(iy, DP) * dy
            
            do ir = 0, params%radius_steps
                r = params%radius_min + real(ir, DP) * dr
                
                call create_slices(xc, yc, r, params%num_slices, &
                                   params%slope_surface_coords, &
                                   params%soils, &
                                   params%all_soil_boundaries, &
                                   params%water_table_coords, &
                                   gamma_water, & 
                                   params%surcharges, & 
                                   slices, is_valid, err_code)
                
                if (is_valid) then
                    ! Always calc Fellenius first (used as initial guess for Bishop)
                    fos_fellenius = calculate_fos_fellenius(slices, params%reinforcements, xc, yc, r)
                    
                    if (trim(method_str) == 'BISHOP') then
                        ! DEBUG CHECK: Uncomment if unsure
                        ! print *, "Running Bishop..." 
                        fos = calculate_fos_bishop(slices, params%reinforcements, xc, yc, r, fos_fellenius)
                    else
                        fos = fos_fellenius
                    end if
                    
                    if (fos < 100.0_DP) then
                        call write_to_log(log_unit, xc, yc, r, fos, "OK")
                    end if
                    
                    valid_circles_count = valid_circles_count + 1
                    
                    if (fos < min_fos) then
                        min_fos = fos
                        crit_xc = xc
                        crit_yc = yc
                        crit_r  = r
                        print '(A, F8.4, A)', " New Min FoS found: ", min_fos, " ..."
                    end if
                else
                    if (err_code == 2) call write_to_log(log_unit, xc, yc, r, 999.0_DP, "INVALID_SOIL_ID")
                end if
                
            end do
        end do
    end do
    
    close(log_unit)

    print *, ""
    print *, "================================================="
    print *, "               SEARCH COMPLETED                  "
    print *, "================================================="
    print *, "Total Valid Circles: ", valid_circles_count
    print *, "Minimum FoS:         ", min_fos
    print *, "Critical Center Xc:  ", crit_xc
    print *, "Critical Center Yc:  ", crit_yc
    print *, "Critical Radius R:   ", crit_r
    print *, "================================================="

    call write_results_summary("results.txt", min_fos, crit_xc, crit_yc, crit_r)
    
    if (allocated(slices)) deallocate(slices)

end program slope_solver
