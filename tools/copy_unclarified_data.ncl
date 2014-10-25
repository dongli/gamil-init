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

    setfileoption(f2, "DefineMode", True)

    do i = 1, 4
        var_name = "TS"+i
        filevardef(f2, var_name, typeof(f1->$var_name$), (/"time","lat","lon"/))
        f2->$var_name$ = f1->$var_name$
        f2->$var_name$@comments = "This variable is copied from old dataset!"
        system("echo [Notice]: Variable "+var_name+" has been copied.")
    end do

    filevardef(f2, "SNOWHICE", typeof(f1->SNOWHICE), (/"time","lat","lon"/))
    f2->SNOWHICE = f1->SNOWHICE
    f2->SNOWHICE@comments = "This variable is copied from old dataset!"
    system("echo [Notice]: Variable SNOWHICE has been copied.")

    filevardef(f2, "TSICE", typeof(f1->TSICE), (/"time","lat","lon"/))
    f2->TSICE = f1->TSICE
    f2->TSICE@comments = "This variable is copied from old dataset!"
    system("echo [Notice]: Variable TSICE has been copied.")

    filevardef(f2, "CWAT", typeof(f1->CWAT), (/"time","lat","lev","lon"/))
    f2->CWAT = f1->CWAT
    f2->CWAT@comments = "This variable is copied from old dataset!"
    system("echo [Notice]: Variable CWAT has been copied.")

end
