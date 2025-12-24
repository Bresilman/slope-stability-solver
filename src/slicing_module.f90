module slicing_module
    use data_structure
    implicit none
    
contains

    function linear_interpolation(polyline, x)
        implicit none
        real(DP), intent(in), dimension(:,:) :: polyline
        real(DP), intent(in) :: x
        real(DP) :: linear_interpolation
        real(DP) :: y
        integer :: i, num_points
        real(DP) :: x1, y1, x2, y2
        
        linear_interpolation = 0.0_DP 
        num_points = size(polyline, dim=1)

        if (x <= polyline(1,1)) then
            y = polyline(1,2)
            linear_interpolation = y
            return
        else if (x >= polyline(num_points,1)) then
            y = polyline(num_points,2)
            linear_interpolation = y
            return
        end if
        
        do i=1 , num_points - 1
            if (x >= polyline(i,1) .and. x <= polyline(i+1,1)) then
                x1 = polyline(i,1)
                y1 = polyline(i,2)
                x2 = polyline(i+1,1)
                y2 = polyline(i+1,2)

                if (abs(x2 - x1) < 1.0E-9_DP) then
                    y = y1
                    linear_interpolation = y
                    return
                end if
                
                y = y1 + (x - x1) * ((y2 - y1) / (x2 - x1))
                linear_interpolation = y
                return
            end if
        end do
        y = polyline(num_points,2)   
        linear_interpolation = y
    end function linear_interpolation

    ! FIXED: Split the long subroutine definition line using '&' to prevent truncation errors
    subroutine create_slices(xc, yc, r, num_slices, slope_coords, all_soils, &
                             all_soil_boundaries, water_coords, gamma_water, &
                             surcharges, slices, is_valid, error_code)
                             
        real(DP), intent(in) :: xc, yc, r
        integer, intent(in) :: num_slices
        real(DP), intent(in), dimension(:,:) :: slope_coords, water_coords
        type(soil_t), intent(in), dimension(:) :: all_soils
        real(DP), intent(in), dimension(:,:,:) :: all_soil_boundaries
        real(DP), intent(in) :: gamma_water
        type(surcharge_t), intent(in), dimension(:), optional :: surcharges
        type(slice_t), allocatable, dimension(:), intent(out) :: slices
        logical, intent(out) :: is_valid
        integer, intent(out), optional :: error_code

        real(DP) :: x_entry, x_exit, total_width, slice_width
        real(DP) :: x_right, x_left, x_mid
        real(DP) :: y_top_left, y_top_right, y_bottom_left, y_bottom_right, y_mid_base
        real(DP) :: h_avg, area, y_water, h_water
        integer :: i, k, soil_id
        real(DP) :: overlap_left, overlap_right, surcharge_weight
        
        is_valid = .false.
        if (present(error_code)) error_code = 0 

        call find_circle_ground_intersections(xc, yc, r, slope_coords, x_entry, x_exit, is_valid)
        if (.not. is_valid) then
            if (present(error_code)) error_code = 1 
            return
        end if

        total_width = x_exit - x_entry  
        slice_width = total_width / real(num_slices, DP)

        if (allocated(slices)) deallocate(slices)
        allocate(slices(num_slices))

        do i = 1, num_slices
            x_left = x_entry + real(i - 1, DP) * slice_width
            x_right = x_left + slice_width
            x_mid = (x_left + x_right) / 2.0_DP

            if ((r**2 - (x_left - xc)**2) < 0.0_DP .or. &
                (r**2 - (x_right - xc)**2) < 0.0_DP .or. &
                (r**2 - (x_mid - xc)**2) < 0.0_DP) then
                is_valid = .false.
                if (present(error_code)) error_code = 1
                return
            end if
            
            y_top_left = linear_interpolation(slope_coords, x_left)
            y_top_right = linear_interpolation(slope_coords, x_right)
            y_bottom_left = yc - sqrt(r**2 - (x_left - xc)**2)
            y_bottom_right = yc - sqrt(r**2 - (x_right - xc)**2)
            y_mid_base = yc - sqrt(r**2 - (x_mid - xc)**2)

            slices(i)%x_mid = x_mid
            slices(i)%y_base = y_mid_base

            if (y_bottom_left > y_top_left + 1.0E-3_DP .or. y_bottom_right > y_top_right + 1.0E-3_DP) then
                 slices(i)%height = 0.0_DP; slices(i)%W = 0.0_DP; slices(i)%cohesion = 0.0_DP
                 slices(i)%friction_ang_rad = 0.0_DP; slices(i)%base_ang_rad = 0.0_DP
                 slices(i)%base_length = 0.0_DP; slices(i)%material_id = 0
                 cycle 
            end if

            h_avg = ((y_top_left - y_bottom_left) + (y_top_right - y_bottom_right)) / 2.0_DP
            if (h_avg < 0.0_DP) h_avg = 0.0_DP
            slices(i)%height = h_avg
            area = slice_width * h_avg

            slices(i)%b = slice_width
            slices(i)%base_length = sqrt((x_right - x_left)**2 + (y_bottom_right - y_bottom_left)**2)
            slices(i)%base_ang_rad = atan2((y_bottom_right - y_bottom_left), (x_right - x_left))

            call get_soil_at_point(x_mid, y_mid_base, all_soils, all_soil_boundaries, soil_id)
            slices(i)%material_id = soil_id 

            if (soil_id == 0) then
                 if (present(error_code)) error_code = 2 
                 is_valid = .false.
                 return
            end if

            if (soil_id >= 1 .and. soil_id <= size(all_soils)) then
                slices(i)%cohesion = all_soils(soil_id)%cohesion
                slices(i)%friction_ang_rad = all_soils(soil_id)%friction_ang_rad
                slices(i)%W = all_soils(soil_id)%gamma * area
            else
                slices(i)%cohesion = 0.0_DP; slices(i)%friction_ang_rad = 0.0_DP; slices(i)%W = 0.0_DP
            end if

            y_water = linear_interpolation(water_coords, x_mid)
            h_water = y_water - y_mid_base
            if (h_water < 0.0_DP) h_water = 0.0_DP
            slices(i)%pore_pressure = h_water * gamma_water
            
            if (present(surcharges)) then
                do k = 1, size(surcharges)
                    overlap_left = max(x_left, surcharges(k)%x_start)
                    overlap_right = min(x_right, surcharges(k)%x_end)
                    if (overlap_right > overlap_left) then
                        surcharge_weight = surcharges(k)%magnitude * (overlap_right - overlap_left)
                        slices(i)%W = slices(i)%W + surcharge_weight
                    end if
                end do
            end if
            
        end do
        is_valid = .true.
    end subroutine create_slices

    subroutine find_circle_ground_intersections(xc, yc, r, polyline, x_min, x_max, found)
        real(DP), intent(in) :: xc, yc, r
        real(DP), intent(in), dimension(:,:) :: polyline
        real(DP), intent(out) :: x_max, x_min
        logical, intent(out) :: found
        integer, parameter :: N_STEPS = 2000 
        integer :: i, num_points
        real(DP) :: x_start, x_end, x_step, current_x, y_ground, y_circle
        logical :: is_below_ground
        
        num_points = size(polyline, dim=1)
        found = .false.
        x_min = huge(1.0_DP); x_max = -huge(1.0_DP) 
        x_start = polyline(1, 1); x_end = polyline(num_points, 1)
        x_step = (x_end - x_start) / real(N_STEPS, DP)
        if (x_step <= 0.0_DP) return

        do i = 0, N_STEPS
            current_x = x_start + (real(i, DP) * x_step)
            if (abs(current_x - xc) >= r) cycle  
            y_ground = linear_interpolation(polyline, current_x)
            if (r**2 - (current_x - xc)**2 < 0.0_DP) cycle
            y_circle = yc - sqrt(r**2 - (current_x - xc)**2)
            is_below_ground = (y_circle <= y_ground + 1.0E-3_DP) 
            if (is_below_ground) then
                found = .true.
                x_min = min(x_min, current_x)
                x_max = max(x_max, current_x)
            end if
        end do
        if (x_min > x_max) found = .false.
    end subroutine find_circle_ground_intersections

    subroutine get_soil_at_point(x, y, all_soils, all_soil_boundaries, soil_id)
        real(DP), intent(in) :: x, y
        type(soil_t), intent(in), dimension(:) :: all_soils
        real(DP), intent(in), dimension(:,:,:) :: all_soil_boundaries  
        integer, intent(out) :: soil_id
        real(DP) :: y_boundary
        integer :: i, num_boundaries, num_points
        real(DP), allocatable, dimension(:,:) :: temp_boundary

        num_boundaries = size(all_soil_boundaries, dim=3)
        num_points = size(all_soil_boundaries, dim=1)
        soil_id = 1 
        if (num_boundaries == 0) return 
        allocate(temp_boundary(num_points, 2))

        do i = 1, num_boundaries
            temp_boundary = all_soil_boundaries(:,:,i)
            y_boundary = linear_interpolation(temp_boundary, x)
            if (y <= y_boundary) then 
                soil_id = i + 1       
            else
                exit
            end if
        end do
        deallocate(temp_boundary)
        if (soil_id > size(all_soils)) soil_id = size(all_soils)
    end subroutine get_soil_at_point
    
end module slicing_module
