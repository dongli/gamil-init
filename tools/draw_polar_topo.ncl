load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    file_names = (/"ic.gamil.360x154x26.0011-01-01-00000.nc", \
                   "ic.gamil.720x308x26.0011-01-01-00000.nc", \
                   "ic.gamil.1440x616x26.0011-01-01-00000.nc"/)
    f = addfiles(file_names, "r")

    wks = gsn_open_wks("pdf", "polar_topo")

    gsn_define_colormap(wks, "OceanLakeLandSnow")

    res = True
    res@gsnFrame = False
    res@gsnDraw = False
    res@gsnPolar = "NH"
    res@mpMinLatF = 70.0
    res@cnFillMode = "RasterFill"
    res@cnFillOn = True
    res@cnLinesOn = False
    res@gsnSpreadColors = True
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/10.0,100.0,500.0,1000.0,1500.0,2000.0,2500.0,3000.0,3500.0,4000.0/)

    plot = new(dimsizes(file_names), graphic)
    do j = 0, dimsizes(file_names)-1
        PHIS = f[j]->PHIS/10.0
        copy_VarCoords(f[j]->PHIS, PHIS)
        plot(j) = gsn_csm_contour_map(wks, PHIS(0,:,:), res)
        delete(PHIS)
    end do

    gsn_panel(wks, plot, (/2,2/), False)

end
