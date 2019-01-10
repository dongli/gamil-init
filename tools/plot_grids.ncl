if (.not. isdefined("grid_file")) then
  system("echo '[Error]: Argument grid_file is needed!'")
  exit
end if

f = addfile(grid_file, "r")

wks = gsn_open_wks("pdf", systemfunc("basename " + grid_file + " .nc"))

min_lat = -90.0

res_map = True
res_map@gsnDraw = True
res_map@gsnFrame = False
;res_map@gsnPolar = "NH"
;res_map@mpMinLatF = min_lat
res_map@mpProjection = "Satellite"
res_map@mpCenterLatF = 45.0

plot = gsn_csm_map(wks, res_map)

center_lon = f->grid_center_lon
center_lat = f->grid_center_lat
corner_lon = f->grid_corner_lon
corner_lat = f->grid_corner_lat
imask = f->grid_imask

res_grid = True
res_grid@gsMarkerSizeF = 0.00001
res_grid@gsMarkerColor = "red"

gsn_polymarker(wks, plot, center_lon, center_lat, res_grid)

res_cell = True
res_cell@gsLineThicknessF = 0.0001

dims = dimsizes(corner_lon)
do i = 0, dims(0) - 1
  if (corner_lat(i,0) .lt. min_lat) then
    continue
  end if
  if (imask(i) .eq. 0) then
    gsn_polygon(wks, plot, corner_lon(i,:), corner_lat(i,:), res_cell)
  else
    gsn_polyline(wks, plot, corner_lon(i,:), corner_lat(i,:), res_cell)
  end if
end do

frame(wks)
