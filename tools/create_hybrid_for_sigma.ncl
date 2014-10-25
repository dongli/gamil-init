load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    sigm = (/0.00133607736846712, 0.00513899915929914, 0.0116468036573691, \
             0.021517142227496,   0.0346602323216384,  0.0503740850440076, \
             0.067136789924738,   0.0824488125980373,  0.0975739588097264, \
             0.115367910909135,   0.136301572054213,   0.160928952769942, \
             0.189901763775095,   0.22398679286417,    0.264086091221295, \
             0.311260796214845,   0.366759384193747,   0.432050593290324, \
             0.508862340965566,   0.599227418892952,   0.695816036901329, \
             0.787020691248154,   0.866740639339174,   0.929434205785941, \
             0.970474653327663,   0.992539942855718/)

    sigi = (/0, \
             0.00267215473693425, 0.00760584358166403, 0.0156877637330741, \
             0.0273465207219179,  0.0419739439213588,  0.0587742261666565, \
             0.0754993536828195,  0.0893982715132551,  0.105749646106198, \
             0.124986175712072,   0.147616968396354,   0.174240937143529, \
             0.205562590406661,   0.242410995321679,   0.285761187120911, \
             0.33676040530878,    0.396758363078715,   0.467342823501934, \
             0.550381858429198,   0.648072979356706,   0.743559094445952, \
             0.830482288050356,   0.902998990627992,   0.95586942094389, \
             0.985079885711436,   1/)

    pt = 2.194

    hyam = sigm
    hybm = 1-sigm

    hyai = sigi
    hybi = 1-sigi

    P0 = pt

    f = addfile("hybpara.26.nc", "c")

    setfileoption(f, "DefineMode", True)

    f_atts = True
    f_atts@comments = "This hybrid sigma-pressure coordinate is equivalent to sigma coordinate"
    fileattdef(f, f_atts)

    dim_names = (/"lev","ilev"/)
    dim_sizes = (/26,27/)
    dim_unlim = (/False,False/)
    filedimdef(f, dim_names, dim_sizes, dim_unlim)

    filevardef(f, "P0", "float", "ncl_scalar")
    filevardef(f, "hyam", "float", "lev")
    filevardef(f, "hybm", "float", "lev")
    filevardef(f, "hyai", "float", "ilev")
    filevardef(f, "hybi", "float", "ilev")

    f->P0 = (/P0/)
    f->hyam = (/hyam/)
    f->hybm = (/hybm/)
    f->hyai = (/hyai/)
    f->hybi = (/hybi/)

end
