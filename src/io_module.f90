! ==============================================================================
! MÓDULO: IO_MODULE
! ==============================================================================
! Este módulo gerencia todas as operações de entrada e saída (I/O) do programa.
! ...
! ==============================================================================

module io_module
    use data_structure
    implicit none
    
    ! Estrutura que contém TODOS os dados lidos do arquivo de entrada.
    type :: input_parameters_t
        character(len=100) :: title          ! Título do projeto
        character(len=20)  :: units          ! Sistema de unidades (ex: METRIC)
        character(len=20)  :: solver_method  ! Método de solução (ex: FELLENIUS)
        character(len=20)  :: weight_method  ! Método de cálculo de peso (ex: SIMPLIFIED)
        integer            :: num_slices     ! Número de fatias a serem usadas
        integer            :: num_soils      ! Quantidade total de tipos de solo
        integer            :: num_reinforcements
        integer            :: num_surcharges
        
        type(soil_t), allocatable, dimension(:) :: soils ! Array com propriedades dos solos
        type(reinforcement_t), allocatable, dimension(:) :: reinforcements
        type(surcharge_t), allocatable, dimension(:) :: surcharges
        
        ! Polilinhas geométricas (Arrays de coordenadas X, Y)
        real(DP), allocatable, dimension(:,:)   :: slope_surface_coords
        real(DP), allocatable, dimension(:,:)   :: water_table_coords
        
        ! Array 3D para fronteiras: (ponto, coord, indice_fronteira)
        real(DP), allocatable, dimension(:,:,:) :: all_soil_boundaries  
        
        ! Parâmetros da grade de busca (Search Grid)
        real(DP) :: grid_top_left(2)      ! Canto superior esquerdo da busca (X, Y)
        real(DP) :: grid_bottom_right(2)  ! Canto inferior direito da busca (X, Y)
        integer  :: grid_steps_x          ! Número de divisões no eixo X
        integer  :: grid_steps_y          ! Número de divisões no eixo Y
        
        ! Parâmetros de busca do raio
        real(DP) :: radius_min            ! Raio mínimo
        real(DP) :: radius_max            ! Raio máximo
        integer  :: radius_steps          ! Número de incrementos de raio
    end type input_parameters_t
    
contains

    ! --------------------------------------------------------------------------
    ! Função Auxiliar: Converte string para maiúsculas
    ! Usada para garantir que a leitura de palavras-chave (keywords) não seja
    ! sensível a maiúsculas/minúsculas (case-insensitive).
    ! --------------------------------------------------------------------------
    function uppercase(input_str) result(output_str)
        implicit none
        character(len=*), intent(in) :: input_str
        character(len=len(input_str)) :: output_str
        integer :: i
        character(len=1) :: char
        
        do i = 1, len(input_str)
            char = input_str(i:i)
            if (char >= 'a' .and. char <= 'z') then
                output_str(i:i) = achar(iachar(char) - 32)
            else
                output_str(i:i) = char
            end if
        end do
    end function uppercase

    ! --------------------------------------------------------------------------
    ! Sub-rotina: read_input_file
    ! Lê o arquivo de texto, identifica seções e valores, e preenche 'params'.
    ! --------------------------------------------------------------------------
    subroutine read_input_file(filename, params)
        character(len=*), intent(in) :: filename
        type(input_parameters_t), intent(out) :: params
        
        integer :: file_unit, i, status
        integer :: slope_coord_count, water_coord_count, soil_count, reinf_count, surch_count
        integer :: num_boundaries, points_per_boundary_count
        integer :: current_boundary_index, point_count_in_boundary
        character(len=200) :: line, keyword, value
        logical :: reading_slope_surface, reading_water_table, reading_soil_boundary
        logical :: is_first_boundary
        integer :: dummy_id, dummy_read_count ! Variável temporária para ler o ID e contagem de leitura
        
        ! Temp vars
        integer :: r_type
        real(DP) :: r_x1, r_y1, r_x2, r_y2, r_t_str, r_s_str, r_space, s_x1, s_x2, s_mag
        
        open(newunit=file_unit, file=filename, status='old', action='read', iostat=status)
        
        if (status /= 0) then
            print *, 'Erro critico! Nao foi possivel abrir o arquivo: ', trim(filename)
            stop
        end if
        
        ! --- PASSO 1: Contagem para Alocação ---
        slope_coord_count = 0
        water_coord_count = 0
        num_boundaries = 0
        points_per_boundary_count = 0
        reading_soil_boundary = .false.
        reading_slope_surface = .false.
        reading_water_table = .false.
        is_first_boundary = .false.
        
        do
            read(file_unit, '(A)', end=100, iostat=dummy_read_count) line
            if (dummy_read_count /= 0) exit
            if (trim(line) == "" .or. line(1:1) == "#" .or. line(1:1) == "!") cycle 
            
            select case (trim(uppercase(line)))
                case ('SLOPE_SURFACE_START'); reading_slope_surface = .true.; cycle
                case ('SLOPE_SURFACE_END');   reading_slope_surface = .false.; cycle
                case ('WATER_TABLE_START');   reading_water_table = .true.; cycle
                case ('WATER_TABLE_END');     reading_water_table = .false.; cycle  
                case ('SOIL_BOUNDARY_START')
                    reading_soil_boundary = .true.
                    num_boundaries = num_boundaries + 1
                    if (num_boundaries == 1) is_first_boundary = .true.
                    cycle
                case ('SOIL_BOUNDARY_END')
                    reading_soil_boundary = .false.
                    is_first_boundary = .false.
                    cycle
            end select  
            
            if (reading_slope_surface) slope_coord_count = slope_coord_count + 1
            if (reading_water_table) water_coord_count = water_coord_count + 1
            if (reading_soil_boundary .and. is_first_boundary) then
                points_per_boundary_count = points_per_boundary_count + 1
            end if
            
        end do
100     continue
        
        ! Alocação dinâmica dos arrays com base nas contagens
        allocate(params%slope_surface_coords(slope_coord_count, 2))
        allocate(params%water_table_coords(water_coord_count, 2))
        if (num_boundaries > 0) then
            allocate(params%all_soil_boundaries(points_per_boundary_count, 2, num_boundaries))
        end if
        
        ! --- PASSO 2: Leitura e Armazenamento dos Dados ---
        rewind(file_unit) ! Volta ao início do arquivo
        
        slope_coord_count = 0
        water_coord_count = 0
        soil_count = 0
        reinf_count = 0
        surch_count = 0
        current_boundary_index = 0
        point_count_in_boundary = 0
        reading_slope_surface = .false.
        reading_water_table = .false.
        reading_soil_boundary = .false.
        
        do
            read(file_unit, '(A)', end=200, iostat=dummy_read_count) line
            if (dummy_read_count /= 0) exit
            if (trim(line) == "" .or. line(1:1) == "#" .or. line(1:1) == "!") cycle
            
            ! Verifica se a linha é um marcador de seção
            select case (trim(uppercase(line)))
                case ('SLOPE_SURFACE_START'); reading_slope_surface = .true.; cycle
                case ('SLOPE_SURFACE_END');   reading_slope_surface = .false.; cycle
                case ('WATER_TABLE_START');   reading_water_table = .true.; cycle
                case ('WATER_TABLE_END');     reading_water_table = .false.; cycle  
                case ('SOIL_BOUNDARY_START')
                    reading_soil_boundary = .true.
                    current_boundary_index = current_boundary_index + 1
                    point_count_in_boundary = 0
                    cycle
                case ('SOIL_BOUNDARY_END');   reading_soil_boundary = .false.; cycle    
            end select  
            
            ! Leitura de coordenadas geométricas (dentro das seções)
            if (reading_slope_surface) then
                slope_coord_count = slope_coord_count + 1
                read(line, *) params%slope_surface_coords(slope_coord_count, :)
                cycle
            end if
            if (reading_water_table) then
                water_coord_count = water_coord_count + 1
                read(line, *) params%water_table_coords(water_coord_count, :)
                cycle
            end if
            if (reading_soil_boundary) then
                point_count_in_boundary = point_count_in_boundary + 1
                if (point_count_in_boundary <= points_per_boundary_count .and. current_boundary_index <= num_boundaries) then
                    read(line, *) params%all_soil_boundaries(point_count_in_boundary, :, current_boundary_index)
                end if
                cycle
            end if
            
            ! Leitura de pares Chave: Valor (KEYWORD: VALUE)
            i = index(line, ':')

            if (i > 0) then
                keyword = trim(uppercase(adjustl(line(1:i-1))))
                value = trim(adjustl(line(i+1:)))
                
                select case (keyword)
                    case ('TITLE'); read(value, '(A)') params%title
                    
                    ! FIX: Ensure UNITS is Uppercase
                    case ('UNITS'); 
                        read(value, '(A)') params%units
                        params%units = uppercase(trim(params%units))
                        
                    ! FIX: Ensure SOLVER is Uppercase ("Bishop" -> "BISHOP")
                    case ('SOLVER'); 
                        read(value, '(A)') params%solver_method
                        params%solver_method = uppercase(trim(params%solver_method))
                        
                    case ('WEIGHT_METHOD'); read(value, '(A)') params%weight_method
                    case ('NUM_SLICES'); read(value, *) params%num_slices
                    case ('NUM_SOILS')
                        read(value, *) params%num_soils
                        if (allocated(params%soils)) deallocate(params%soils)
                        allocate(params%soils(params%num_soils))
                    case ('SOIL')
                        soil_count = soil_count + 1
                        if (soil_count <= params%num_soils) then
                            ! Lê ID (descartado), Coesão, Atrito e Peso
                            read(value, *) dummy_id, &
                                             params%soils(soil_count)%cohesion, &
                                             params%soils(soil_count)%friction_ang_rad, &
                                             params%soils(soil_count)%gamma
                            ! Converte graus para radianos imediatamente - SPLIT LINE FOR FORTRAN LENGTH LIMIT
                            params%soils(soil_count)%friction_ang_rad = &
                                params%soils(soil_count)%friction_ang_rad * (acos(-1.0_DP) / 180.0_DP)
                        else
                            print *, 'Aviso: Mais solos definidos do que o valor de NUM_SOILS.'
                        end if
                    
                    case ('NUM_REINFORCEMENTS')
                        read(value, *) params%num_reinforcements
                        if (allocated(params%reinforcements)) deallocate(params%reinforcements)
                        if (params%num_reinforcements > 0) allocate(params%reinforcements(params%num_reinforcements))
                    case ('REINFORCEMENT')
                        reinf_count = reinf_count + 1
                        if (reinf_count <= params%num_reinforcements) then
                            ! Defaults
                            params%reinforcements(reinf_count)%bond_str = 0.0_DP
                            params%reinforcements(reinf_count)%friction_coeff = 0.0_DP 
                            
                            ! Try reading 11 values
                            read(value, *, iostat=status) dummy_id, r_type, r_x1, r_y1, r_x2, r_y2, r_t_str, r_s_str, r_space, &
                                                          params%reinforcements(reinf_count)%bond_str, &
                                                          params%reinforcements(reinf_count)%friction_coeff
                            
                            if (status /= 0) then
                                ! Fallback: Read 10 values
                                read(value, *, iostat=status) dummy_id, r_type, r_x1, r_y1, r_x2, r_y2, r_t_str, r_s_str, r_space, &
                                                              params%reinforcements(reinf_count)%bond_str
                                if (status /= 0) then
                                    ! Fallback: Read 9 values
                                    read(value, *) dummy_id, r_type, r_x1, r_y1, r_x2, r_y2, r_t_str, r_s_str, r_space
                                    params%reinforcements(reinf_count)%bond_str = 0.0_DP
                                end if
                                params%reinforcements(reinf_count)%friction_coeff = 0.0_DP
                            end if
                            
                            params%reinforcements(reinf_count)%id = dummy_id
                            params%reinforcements(reinf_count)%type_id = r_type
                            params%reinforcements(reinf_count)%x_start = r_x1
                            params%reinforcements(reinf_count)%y_start = r_y1
                            params%reinforcements(reinf_count)%x_end = r_x2
                            params%reinforcements(reinf_count)%y_end = r_y2
                            params%reinforcements(reinf_count)%tensile_str = r_t_str
                            params%reinforcements(reinf_count)%shear_str = r_s_str
                            params%reinforcements(reinf_count)%spacing = r_space
                            
                            if (abs(r_x2 - r_x1) > 1.0E-9_DP) then
                                params%reinforcements(reinf_count)%orientation_rad = atan2(r_y2 - r_y1, r_x2 - r_x1)
                            else
                                params%reinforcements(reinf_count)%orientation_rad = 1.57079632679_DP
                            end if
                        end if

                    case ('NUM_SURCHARGES')
                        read(value, *) params%num_surcharges
                        if (allocated(params%surcharges)) deallocate(params%surcharges)
                        if (params%num_surcharges > 0) allocate(params%surcharges(params%num_surcharges))
                    case ('SURCHARGE')
                        surch_count = surch_count + 1
                        if (surch_count <= params%num_surcharges) then
                            read(value, *) s_x1, s_x2, s_mag
                            params%surcharges(surch_count)%x_start = s_x1
                            params%surcharges(surch_count)%x_end = s_x2
                            params%surcharges(surch_count)%magnitude = s_mag
                        end if

                    case ('GRID_TOP_LEFT'); read(value, *) params%grid_top_left
                    case ('GRID_BOTTOM_RIGHT'); read(value, *) params%grid_bottom_right
                    case ('GRID_STEPS_X'); read(value, *) params%grid_steps_x
                    case ('GRID_STEPS_Y'); read(value, *) params%grid_steps_y
                    case ('RADIUS_MIN'); read(value, *) params%radius_min
                    case ('RADIUS_MAX'); read(value, *) params%radius_max
                    case ('RADIUS_STEPS'); read(value, *) params%radius_steps
                    case default
                        print *, 'Aviso: Palavra-chave desconhecida no input: ', trim(keyword)
                end select
            end if
        end do
200     continue
        close(file_unit)
        
    end subroutine read_input_file

    ! --------------------------------------------------------------------------
    ! Sub-rotina: write_results_summary
    ! Escreve o relatório final com o Fator de Segurança Mínimo encontrado.
    ! --------------------------------------------------------------------------
    subroutine write_results_summary(filename, min_fos, crit_xc, crit_yc, crit_r)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: min_fos, crit_xc, crit_yc, crit_r
        integer :: file_unit, status

        open(newunit=file_unit, file=filename, status='replace', action='write', iostat=status)

        if (status /= 0) then
            print *, 'Erro: Nao foi possivel criar arquivo de saida: ', filename
            stop
        end if

        write(file_unit, '(A)') '======================================'
        write(file_unit, '(A)') ' Relatorio de Estabilidade de Taludes '
        write(file_unit, '(A)') '======================================'
        write(file_unit, *) ""
        write(file_unit, '(A, F8.4)')   'Fator de Seguranca Minimo (FoS): ', min_fos
        write(file_unit, *) ""
        write(file_unit, '(A)')         'Parametros do Circulo Critico:'
        write(file_unit, '(A, F12.4)')  '  Centro Xc: ', crit_xc
        write(file_unit, '(A, F12.4)')  '  Centro Yc: ', crit_yc
        write(file_unit, '(A, F12.4)')  '  Raio R:    ', crit_r
        write(file_unit, '(A)') '======================================'

        close(file_unit)
        print *, 'Resumo dos resultados salvo em: ', trim(filename)
    end subroutine write_results_summary

    ! --------------------------------------------------------------------------
    ! Sub-rotina: open_log_file
    ! Inicializa o arquivo de log CSV.
    ! --------------------------------------------------------------------------
    subroutine open_log_file(filename, file_unit)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: file_unit
        integer :: status

        open(newunit=file_unit, file=filename, status='replace', action='write', iostat=status)

        if (status /= 0) then
            print *, 'Erro: Nao foi possivel criar arquivo de log: ', filename
            stop
        end if
        ! Escreve o cabeçalho CSV
        write(file_unit, '(A)') 'Xc,Yc,Radius,FoS,Status'
    end subroutine open_log_file
            
    ! --------------------------------------------------------------------------
    ! Sub-rotina: write_to_log
    ! Escreve uma linha de dados no arquivo de log CSV aberto.
    ! --------------------------------------------------------------------------
    subroutine write_to_log(file_unit, xc, yc, r, fos, status_msg)
        integer, intent(in) :: file_unit
        real(DP), intent(in) :: xc, yc, r, fos
        character(len=*), intent(in) :: status_msg

        ! Formato CSV: Floats separados por vírgula
        write(file_unit, '(F12.4, ",", F12.4, ",", F12.4, ",", F12.4, ",", A)') xc, yc, r, fos, trim(status_msg)
    end subroutine write_to_log

end module io_module
