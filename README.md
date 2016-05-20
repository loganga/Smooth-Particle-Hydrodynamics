# Smooth-Particle-Hydrodynamics
Project for CSSE 451 by Gene Logan and Jack Porter

#Overview
This is a real-time implementation of smooth particle hydrodynamics (SPH). This is essentially a CPU simulation of the SPH method that solves a system of differential equations and renders particles with OpenGL.

#Running the Simulation
To run the code, download the repo and open the project in Visual Studio. Build in release mode to produce an executable or run directly in the IDE if you prefer, but know that tweaking physical constants in the equation must be done directly in the code.

#The Method
In a nutshell, our implementation is a relatively simple process that, for each time step, determines acceleration as a function of all interacting forces on a given particle. In this case, forces are affected by physical qualities at each particle such as the radius, viscosity, bulk modulus, and density. With the calculated acceleration, displacement is determined and uploaded to the GPU. For detailed information on the equations used in this project, refer to the "Smooth Particle Hydrodynamics" section of the following paper: http://matthias-mueller-fischer.ch/publications/sca03.pdf

#Performance
Our method does not make use of trees or spatial partitioning as a means of acceleration.  It does, however, gain a slight performance boost by symetrically updating particles as it iterates through all particles' interactions on a given particle. As a result, 100 particles render around 45 fps, 150 around 14 fps, and 200 around 6 fps.
