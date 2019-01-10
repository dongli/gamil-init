load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    file_names = systemfunc("ls ic.gamil.*")
    f = addfiles(file_names, "r")

    wks = gsn_open_wks("pdf", "lat_grids")

    res = True
    res@gsnFrame = False
    res@gsnDraw = False
    res@gsnMaximize = True
    res@vpWidthF = 0.8
    res@vpHeightF = 0.4
    res@tmYLMode = "Manual"
    res@tmYLTickStartF = 0.0
    res@tmYLTickEndF = 3.0
    res@tmYLTickSpacingF = 0.2
    res@trYMinF = res@tmYLTickStartF
    res@xyLineThicknessF = 3.0
    colors = (/"black","green","red","blue"/)

    plot = new(dimsizes(file_names), graphic)
    do j = 0, dimsizes(file_names)-1
        lat = new(dimsizes(f[j]->lat)-1, typeof(f[j]->lat))
        dlat = new(dimsizes(lat), typeof(f[j]->lat))
        do i = 0, dimsizes(f[j]->lat)-2
            lat(i) = (f[j]->lat(i)+f[j]->lat(i+1))/2.0
            dlat(i) = f[j]->lat(i+1)-f[j]->lat(i)
        end do
        res@xyLineColor = colors(j)
        plot(j) = gsn_csm_xy(wks, lat, dlat, res)
        delete(lat)
        delete(dlat)
    end do

    do j = 1, dimsizes(plot)-1
        overlay(plot(0), plot(j))
    end do
    draw(plot(0))
    frame(wks)

end
