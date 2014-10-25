program main

    use model_grids

    call model_grids_read("../grids/grid.gamil.128x60.nc")
    call model_grids_write("test_model_grids.nc")

end program main
