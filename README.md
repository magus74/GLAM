# GLAM
Glycogen Lactate Absorption Model

source code tested on Linux Gentoo (gcc 4.8.4) and on Windows 8 (Microsoft Visual Studio 15)

Usage:
./glam.exe --help gives the following output:

Allowed options:
  -h [ --help ]                    produce help message
  -i [ --input ] arg               set input file (ply)
  -s [ --source ] arg              set source file (ply)
  -b [ --laplacian ]               compute laplacian field (default false)
  -w [ --weight ] arg              set weight max
  -r [ --radius ] arg              set influence radius (default 0.5)
  -l [ --lambertian ]              set lambertian light model (default false =
                                   gaussian)
  -c [ --colormap ] arg            set colormap ascii file (t_norm R_byte
                                   G_byte B_byte )
  -t [ --threshold ] arg           set normalized threshold for clustering
  -p [ --peaks ] arg               output cluster peak ply file
  -d [ --cluster-data ] arg        output cluster data file
  -f [ --object-data ] arg         output object data filename
  -x [ --colorized ] arg           output cluster peak ply filename
  -m [ --colormap-per-object ] arg colormap per object ascii file (t_norm
                                   R_byte G_byte B_byte )
  -e [ --expected-absorption ] arg set maximum expected absorption
  -o [ --output ] arg              output file ply


Dependecies:
* Boost 
* SpaceLand (https://sourceforge.net/projects/spacelib)
* PLY and PLY I/O

Instructions:
* compile ply and ply_io by using cmake (CMake gui for generating Visual Studio solution, cmake for generating Makefile for Unix systems)
* set up SL, ply, and ply_io libraries and update CMakeLists.txt accordingly 
