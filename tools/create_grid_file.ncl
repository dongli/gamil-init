load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  ic_path = "./ic.gamil.128x60x26.0011-01-01-00000.nc"

  f = addfile(ic_path, "r")

  lat = f->lat
  lat_bnds = f->lat_bnds
  lon = f->lon
  lon_bnds = f->lon_bnds

  num_lat = dimsizes(lat)
  num_lon = dimsizes(lon)
  grid_size = num_lon*num_lat

  grid_dims = (/num_lon,num_lat/)
  grid_imask = new(grid_size, integer)
  grid_imask@units = "unitless"
  grid_center_lat = new(grid_size, double)
  grid_center_lat@units = "degrees"
  grid_center_lon = new(grid_size, double)
  grid_center_lon@units = "degrees"
  grid_corner_lat = new((/grid_size,4/), double)
  grid_corner_lat@units = "degrees"
  grid_corner_lon = new((/grid_size,4/), double)
  grid_corner_lon@units = "degrees"

  grid_imask = 1
  k = 0
  do j = 0, num_lat-1
    do i = 0, num_lon-1
      grid_center_lat(k) = lat(j)
      grid_center_lon(k) = lon(i)
      if (j .ne. 0) then
        grid_corner_lat(k,:) = (/lat_bnds(j-1),lat_bnds(j-1),lat_bnds(j),lat_bnds(j)/)
      else
        grid_corner_lat(k,:) = (/-90.0d0,-90.0d0,lat_bnds(j),lat_bnds(j)/)
      end if
      if (i .ne. 0) then
        grid_corner_lon(k,:) = (/lon_bnds(i-1),lon_bnds(i),lon_bnds(i),lon_bnds(i-1)/)
      else
        grid_corner_lon(k,:) = (/lon_bnds(num_lon-1),lon_bnds(i),lon_bnds(i),lon_bnds(num_lon-1)/)
      end if
      k = k+1
    end do
  end do

  grid_path = "grid.gamil." + sprinti("%d", num_lon) + "x" + sprinti("%d", num_lat) + ".nc"

  system("rm -rf " + grid_path)
  f = addfile(grid_path, "c")

  setfileoption(f, "DefineMode", True)

  f_att = True
  f_att@history = "Created by GAMIL-INIT (" + systemfunc("date") + ")"
  f_att@author = "Li Dong <dongli@lasg.iap.ac.cn>"
  fileattdef(f, f_att)

  dim_names = (/"grid_rank","grid_size","grid_corners"/)
  dim_sizes = (/2,grid_size,4/)
  dim_unlim = (/False,False,False/)
  filedimdef(f, dim_names, dim_sizes, dim_unlim)

  filevardef(f, "grid_dims", "integer", "grid_rank")
  filevardef(f, "grid_imask", "integer", "grid_size")
  filevardef(f, "grid_center_lat", "double", "grid_size")
  filevardef(f, "grid_center_lon", "double", "grid_size")
  filevardef(f, "grid_corner_lat", "double", (/"grid_size","grid_corners"/))
  filevardef(f, "grid_corner_lon", "double", (/"grid_size","grid_corners"/))

  filevarattdef(f, "grid_imask", grid_imask)
  filevarattdef(f, "grid_center_lat", grid_center_lat)
  filevarattdef(f, "grid_center_lon", grid_center_lon)
  filevarattdef(f, "grid_corner_lat", grid_corner_lat)
  filevarattdef(f, "grid_corner_lon", grid_corner_lon)

  setfileoption(f, "DefineMode", False)

  f->grid_dims = (/grid_dims/)
  f->grid_imask = (/grid_imask/)
  f->grid_center_lat = (/grid_center_lat/)
  f->grid_center_lon = (/grid_center_lon/)
  f->grid_corner_lat = (/grid_corner_lat/)
  f->grid_corner_lon = (/grid_corner_lon/)

end
