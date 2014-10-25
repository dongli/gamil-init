load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    system("echo Input the file name:")
    file_name = systemfunc("read -e -p '> ' ans && echo $ans")

    system("echo Input variable name:")
    var_name = systemfunc("read -e -p '> ' ans && echo $ans")

    f = addfile(file_name, "r")

    wks = gsn_open_wks("pdf", var_name)

    gsn_define_colormap(wks, "rainbow")

    res = True
    res@gsnMaximize = True
    res@gsnSpreadColors = True
    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnFillMode = "CellFill"
    res@lbLabelAngleF = 45.0

    num_dim = dimsizes(dimsizes(f->$var_name$))
    if (num_dim .eq. 3) then
        plot = gsn_csm_contour_map_ce(wks, f->$var_name$(0,:,:), res)
    end if
    if (num_dim .eq. 4) then
        system("echo Input the level to plot:")
        k = stringtointeger(systemfunc("read -p '> ' ans && echo $ans"))
        plot = gsn_csm_contour_map_ce(wks, f->$var_name$(0,:,k,:), res)
    end if

end
