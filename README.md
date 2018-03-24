# Stickkit

Copyright 2007,08,13-15,17 Mark J. Stock mstock@umich.edu

## Introduction

Stickkit is a command-line tool for manipulating segmented vector geometry in 2D or 3D space. It is a single ANSI C file and requires no special libraries to build on most machines. It fills a gap that I found in command-line geometry processing software: a tool designed specifically for segmented geometries (graphs and networks).

Note that similar software is available for three-dimesional triangle meshes: first my own Rocktools, and also the very powerful tool Meshlab---a GUI triangle mesh manipulation tool that can read and write a number of formats, plus perform a large variety of smoothing and simplification operations.
Files

# Get and build

The current version of Stickkit is 0.5. On a Linux or OSX system, build it with:

    sudo yum install wxGTK-devel wxsvg-devel
    git clone https://github.com/markstock/stickkit.git
    cd stickkit
    make

or, to cross-compile for Windows on Linux, use:

    /usr/local/bin/i386-mingw32-gcc -O2 -o stickkit.exe stickkit.c -lm

# Features

Stickkit can do a few interesting things:

* Reads Wavefront-.obj-like files (which I call ".seg" to differentiate from the tri-mesh standard) with tangents and radii and reads Radiance lists of cylinders, cones, and spheres.
* All node locations, tangent vectors, and radii are dynamically stored in a binary tree and duplicate entries are quickly found and removed, so meshes using 10-20 million segments load, process, and write relatively fast.
* Basic operations include multi-axis scaling and translation.
* Coarsen operation removes nodes such that all segments have a minimum length. New nodes are placed according to spline interpolation of the original paths.
* Likewise, the refine operation adds new nodes along spline-interpolated paths such that no segment exceeds a maximum length.
* Roughen operation displaces each node according to a Gaussian (diffusion-like) distribution.
* Split operation splits input geometry smoothly along its longest axis and writes two separate files.
* Most operations can be stacked, for example: -coarsen 0.01 -roughen 0.3 -refine .005 -roughen 0.3 -refine .002 -roughen 0.3 and will be performed in order.
* VTK files can be opened in ParaView, an open source scientific visualization package from Kitware/ASC/Sandia/LANL/ARL. Most of the images below were made using ParaView. 

# To do

* Allow writing of accurate png files
* Use the command-line options from rockxray for aiming
* Use the voxel stuff from bob.c for pixel values
* Use png-writing code from rockxray for the output
* Get stickkit to dump .rad segments and nodes as 10 different "def6"-type things
* That way, I can turn def0 into "light" and the others into "glow"

