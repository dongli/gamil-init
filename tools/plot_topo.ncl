load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  if (.not. isdefined("ic_file")) then
    system("echo '[Error]: Argument ic_file is needed!'")
    exit
  end if

	f = addfile(ic_file, "r")

	wks = gsn_open_wks("pdf", "topo")

	res = True
	res@cnLinesOn = False
	res@cnFillOn = True

	plot = gsn_csm_contour_map(wks, f->PHIS(0,:,:), res)

  plot = gsn_csm_contour_map(wks, f->SGH(0,:,:), res)

  plot = gsn_csm_contour_map(wks, f->LANDFRAC(0,:,:), res)

end
