program test_interp_aero

    use source_data
    use model_grids
    use model_fc

    call calc_model_grids(128, 60)

    call source_data_start_aerosol("../data/aero_1.9x2.5_L26_1850-2009.nc")
    call model_fc_interp_aerosol
    call source_data_end_aerosol

end program test_interp_aero
