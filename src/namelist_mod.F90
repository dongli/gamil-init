module namelist_mod

  implicit none

  character(300) :: model_grid_file = "N/A"
  character(300) :: uvtq_data_file  = "N/A"
  character(300) :: ozone_data_file = "N/A"
  character(300) :: aero_data_file = "N/A"
  character(300) :: topo_data_file  = "N/A"
  character(300) :: sst_data_file = "N/A"
  character(300) :: ice_data_file = "N/A"
  logical :: correct_sst_sice = .true.
  logical :: model_ic_allow_extrap = .false.
  integer num_model_lon, num_model_lat
  real(8) :: latmesh_B = 2
  integer :: uvtq_data_type = -1

  integer, parameter :: equal_interval_grid = 1
  integer, parameter :: even_area_grid = 2
  integer :: model_grid_type = even_area_grid

  namelist /gamil_init_config/  &
    num_model_lon             , &
    num_model_lat             , &
    model_grid_type           , &
    uvtq_data_type            , &
    model_grid_file           , &
    uvtq_data_file            , &
    ozone_data_file           , &
    topo_data_file            , &
    aero_data_file            , &
    sst_data_file             , &
    ice_data_file             , &
    model_ic_allow_extrap     , &
    latmesh_B                 , &
    correct_sst_sice

end module namelist_mod
