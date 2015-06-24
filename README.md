# talfer
Fluid simulation game based on the Lattice Boltzmann Method (LBM) where the fluid interacts with rigid bodies and vice versa. The rigid bodies can be steered by hand movements recognized by the Kinect.


For physically correct Poiseuille flow:
./build/fa_2013_release -p8 -i data/poiseuille_small.png -a 1 -x 0.0

# compiling
There are several compiling instructions available by a simple Makefile:

cd path/to/talfer

make <kinect> <debug> <release>

Check the "run" script to know how to start the game.

