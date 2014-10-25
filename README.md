Description
===========

This package is used to generate the initial condition of GAMIL (or includes
forcing data in future). There are two data sources:

* reanalysis data, e.g. ERA-interim
* model initial condition of different resolution

The second source can be used to change the original settings, such as the
topography, resolution, etc.

Usage
=====

* Compile the codes

Use `CMake` to generate `Makefile`. Note use `gfortran` since the codes contain
some Fortran 2003 syntax, and other Fortran compilers are not happy with it.

```
$ cd build
$ FC=gfortran cmake ..
$ make
```

* Edit namelist

When reanalysis data is used, the original data should be organized into one
file, the utility script will be added in future. Then, edit the namelist
(there is a template can be copied and edited, called "namelist-era").

When existed model initial condition is used, just edit the namelist (the
template is called "namelist-model").

Then run:

```
$ ./gamil-init <namelist>
```

* Copy unclaried data (for the time being)

There are several variables that I do not know how to set, then copy them from
old initial condition. Run NCL script:

```
$ ncl <GAMIL-INIT ROOT>/tools/copy_unclarified_data.ncl
```

Authors
=======

* Li Dong <dongli@lasg.iap.ac.cn>
* Wenyu Huang
* Fabo Zhang
