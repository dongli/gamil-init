load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    fice = addfile("HadISST_ice.nc", "w")

    date = cd_calendar(fice->time, -2)
    copy_VarCoords(fice->time, date)
    fice->date = date

    fice->sic = lonFlip(fice->sic)
    fice->sic = fice->sic(:,::-1,:)
    fice->sic@_FillValue = fice->sic@missing_value

    fsst = addfile("HadISST_sst.nc", "w")

    fsst->time = fice->time
    fsst->date = date

    fsst->sst = lonFlip(fsst->sst)
    fsst->sst = fsst->sst(:,::-1,:)
    fsst->sst@_FillValue = fsst->sst@missing_value

end
