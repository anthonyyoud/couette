Time-series plots
=================
The data for time-series is output in such a way so as to be readable by any
plotting program which can read columnar ascii text.  As of writing this, the
time-series data is output to the file `u_growth.dat`.  It contains 12 columns
which are:

1. time
2. radial velocity at mid-height, mid-gap
3. radial velocity at z=0, mid-gap
4. growth rate of the radial velocity
5. growth rate of the axial velocity
6. axial velocity at mid-height, mid-gap
7. stream function at mid-height, mid-gap
8. total azimuthal velocity at mid-height, mid-gap
9. azimuthal vorticity at quarter-height, mid-gap
10. azimuthal magnetic field at quarter-height, mid-gap
11. azimuthal current at mid-height, mid-gap
12. Reynolds number.

Any one or more of these columns can be plotted using
[gnuplot](http://www.gnuplot.info/), for example, with a command of the form:

    p "u_growth.dat" u 1:n w l

which plots column `n` versus column `1` using lines as the style.

Other common tasks include plotting two columns on the same graph which can be
done with:

    p "u_growth.dat" u 1:n w l, "" u 1:m w l

which plots column `n` versus column `1` and column `m` versus column `1` on
the same graph with the same axis range.  A more sophisticated command which
does the same plot but uses the left-hand y axis for the column `n` data and
the right-hand y axis for the column `m` data is as follows:

    p "u_growth.dat" u 1:n ax x1y1 w l, "" u 1:m ax x1y2 w l

Contour plots
=============
Data for contour plots is stored in the files containing a sequence of numbers.
The sequence of numbers represents the time at which a particular snapshot was
taken.  The actual (dimensionless) diffusion time at which the snapshots are
taken can be found by multiplying by the time step.

All files except `xsect1234567.dat` contain data which can be plotted using
gnuplot.  The prefix of the file determines what data is contained within, eg.
`j` for current, `vr` for radial velocity, `p` for the stream function psi,
etc.

To plot the data use the command

    sp "<filename>" u 1:2:3 w l

which by default shows a surface plot with no contours.  Contours can be
enabled with

    set contour

and the number of contour levels can be set with

    set cntrparam level n

where `n` is the number of levels.  In gnuplot v4.0 and above the surface can
be shaded with the command

    set pm3d

or the contours themselves can be shaded using

    set pm3d map

which then disables the surface and shows a 2D (shaded) contour plot.

There are many other options to control the display of surface/contour plots.
gnuplot's built-in help (accessed with the `help` command) has good
descriptions of all the possibilities.  Start with `help splot` for a
description of surface plots.

Contour plots of binary data
----------------------------
Additional files named `xsect*.dat` are also output which allows the
possibility of animations of contour plots using IDL (Interactive Data
Language).  These files are read in to IDL using the `vcsect.pro`, `jpeg.pro`,
and/or `jpeg_vcsect.pro` IDL procedures (not currently part of the code).
These procedures can be used to output jpeg images of the contour plots which
can then be converted into an animation using either `mencoder` or `convert`.

To convert a series of jpeg images into an avi animation using mencoder (mainly
so that animations can easily be played in Windows) use the command

    mencoder "mf://*.jpg" -mf fps=25 -o out.avi -ovc lavc -lavcopts
        vcodec=msmpeg4v2:vhq (pass=n) -nosound -noaspect

`Pass = n` allows two passes of the input files to enable (supposedly) better
quality.  First call mencoder with `pass = 1` and no output option; then call
it again with `pass = 2` and the output option.

A similar task can be performed if the animation has already been created in
mpeg format.  Then the command

    mencoder -ovc lavc -lavcopts vcodec=msmpeg4v2:vhq (pass=n) -nosound
        -ofps 25 -noaspect in.mpg -o out.avi

will convert in.mpg to out.avi (again with the option of two passes).

Another option of converting the jpeg files into an animation is to use the
`convert` command which is part of the
[ImageMagick](http://www.imagemagick.org) set of programs.  The simple command
is

    convert -quality n *.jpg out.mpg

which converts all jpeg files into `out.mpg` with a quality option where `n` is
an integer between 0 and 100, with 100 being best quality.

The mencoder option generally gives smaller file sizes but lower quality
output.
