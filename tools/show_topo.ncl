load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    f = addfile("../data/topo.2min.nc", "r")

    num_lon = dimsizes(f->lon)

    do i = 0, num_lon-2
        system("echo "+(f->lon(i+1)-f->lon(i))+" "+f->lon(i))
    end do

end
