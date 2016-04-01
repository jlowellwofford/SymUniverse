README
======

SymUniverse is a flexible n-body simulation code.

This project lives on GitHub: [SymUniverse](https://github.com/jlowellwofford/SymUniverse).

** This software is licensed under the GNU Lesser General Public License v 3.0 (see LICENSE file). **

Overview
--------

SymUniverse includes a symulator (`sym`), various universe building tools (`u*`), and analysis tools (`a*`).

The core of SymUniverse is `sym`.  `sym` doesn't actually do any physics.  Rather, it is a pluggable interface
that iterates through "Slices" (i.e. time slices) of "Universes", and then passes the Slices to transformation
modules.  The modules can do essentially anything to the slices, but the point is for them to perform physically
relevant transformations, such as calculating accelerations due to forces, moving the particles based on their
velocity, or detecting and resolving collisions.

Universes are stored in Universe Data Files.  This is a proprietry binary file format that is described in the
file "universe.h".  Universe Slices do not have a fixed particle matrix size; particles can be added or removed
by transformation modules.

Finally, particles can have user-definable flags.  Certain modules can filter based on these flags.  This allows
for different fundamental physics based on file flags.  For instance, a user could use one flag for photons and one
for electrons, then evolve the system using different modules for photons and electrons.

Downloading & Compiling
-----------------------

### You need the following tools:

* git (version control)
* compiler toolchain (e.g. Xcode for Mac OS X)
* CMake (see [CMake](http://cmake.org))

### Downloading & compiling instructions (Mac OS X):

1. `git clone git@github.com:jlowellwofford/SymUniverse.git`
2. `cd SymUniverse`
3. `mkdir xcode ; cd xcode`
4. `/Applications/CMake.app/Contents/bin/cmake -G "Xcode" ..`
5. `xcodebuild`

Binaries are placed in `xcode/src/{bin,modules,lib}/Debug`.

** SymUniverse has not yet been ported to Linux (but will be soon). **

Using SymUniverse
-----------------

### Creating a universe

You can use the following universe creation tools.  See -h option for specifics for each.

* `ubuild`      - build a new universe based on a set of parameters (e.g. particle/velocity distribution, etc.)
                  ubuild's options are still rapidly evolving.
* `umerge`      - merge two universes. An example use case is to setup two types of particles in seperate `ubuilds`,
                    then merge them into one universe. (umerge is presently incomplete)
* `ucat`        - concattenate two universes together. (ucat is presently incomplete)
* `udivide`     - slice up a universe. (udivide is presently incomplete)
* `uextract`    - extract specified slices from a universe and put them in a new universe. (uextract is presently incomplete).

### Running a sym

** Make sure you specify the path to your modules with -M <path> (e.g. -M ../../modules/Debug) **

To see all of the available options, run:

`sym -M <mod_path> -h`

This will give options for both sym and all available modules.

The general syntax is:

`sym -M <mod_path> -i <in_file> -o <out_file> -t <num_steps> -m mod1[op1=val1,op2,...] -m mod2 -mod3[opt1=val1,...]`

Notes:

* Modules are added to the pipeline in command line order.
* If input and output files are the same, sym resumes after the last Slice in the file.
* By default, num_steps = -1, meaning infinite.  `sym` can exit safely, finishing the current step, using Ctrl^C.
  A second Ctrl^C causes an immediate quit.

### Analysing universes

Analysis tools haven't been created yet.  Coming soon!

### An example simulation

Say we want to evolve 1000 bodies (we're not going to worry about physical constants for the example so we'll leave them
as defaults).  We want gravitational interactions and hard sphere collisions and 100 timesteps.  We could do the following:

1. `ubuild -n 1000 -T 300 -r 1E-6 -b 0,1 -B 0,1 -o my.univ`
2. `sym -M ../../modules/Debug -t 100 -i my.univ -o my.univ -m fgrav[cleara=1,plummer=0] -m integrate[boundary=none] -m hscollide -m boundary[boundary=elastic]`

The result will be stored in my.univ.

