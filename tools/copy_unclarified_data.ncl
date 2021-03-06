load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    if (.not. isdefined("file_name1")) then
        system("echo Input the base initial condition:")
        file_name1 = systemfunc("read -e -p '> ' ans && echo $ans")
    end if

    if (.not. isdefined("file_name2")) then
        system("echo Input the generated initial condition:")
        file_name2 = systemfunc("read -e -p '> ' ans && echo $ans")
    end if

    f1 = addfile(file_name1, "r")
    f2 = addfile(file_name2, "rw")

    lon1 = f1->lon
    lat1 = f1->lat
    lon2 = f2->lon
    lat2 = f2->lat

    setfileoption(f2, "DefineMode", True)

    do i = 1, 4
        var_name = "TS"+i
        filevardef(f2, var_name, typeof(f1->$var_name$), (/"time","lat","lon"/))
        f2->$var_name$ = linint2(lon1, lat1, f1->$var_name$, True, lon2, lat2, 0)
        f2->$var_name$@long_name = f1->$var_name$@long_name
        f2->$var_name$@units = f1->$var_name$@units
        f2->$var_name$@comments = "This variable is copied from old dataset!"
        system("echo [Notice]: Variable "+var_name+" has been copied.")
    end do

    filevardef(f2, "SNOWHICE", typeof(f1->SNOWHICE), (/"time","lat","lon"/))
    f2->SNOWHICE = linint2(lon1, lat1, f1->SNOWHICE, True, lon2, lat2, 0)
    f2->SNOWHICE@long_name = f1->SNOWHICE@long_name
    f2->SNOWHICE@units = f1->SNOWHICE@units
    f2->SNOWHICE@comments = "This variable is copied from old dataset!"
    system("echo [Notice]: Variable SNOWHICE has been copied.")

    filevardef(f2, "TSICE", typeof(f1->TSICE), (/"time","lat","lon"/))
    f2->TSICE = linint2(lon1, lat1, f1->TSICE, True, lon2, lat2, 0)
    f2->TSICE@long_name = f1->TSICE@long_name
    f2->TSICE@units = f1->TSICE@units
    f2->TSICE@comments = "This variable is copied from old dataset!"
    system("echo [Notice]: Variable TSICE has been copied.")

    filevardef(f2, "CWAT", typeof(f1->CWAT), (/"time","lat","lev","lon"/))
    do k = 0, dimsizes(f1->lev)-1
        f2->CWAT(:,:,k,:) = linint2(lon1, lat1, f1->CWAT(:,:,k,:), True, lon2, lat2, 0)
    end do
    f2->CWAT@long_name = f1->CWAT@long_name
    f2->CWAT@units = f1->CWAT@units
    f2->CWAT@comments = "This variable is copied from old dataset!"
    system("echo [Notice]: Variable CWAT has been copied.")

end
