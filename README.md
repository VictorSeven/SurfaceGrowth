# SurfaceGrowth


Simulations about various surface growth models done for the subject Cooperative and Critical Phenomena of the Master in Physics of Complex Systems of the Institute for Cross-Disciplinary Physics and Complex Systems.

The C++ code is able to simulate various growth models: diffusion and ballistic, as well as a numerical solution to the KPZ equation via the Milstein algorithm. The functions can be called to do simply the hard computations, or also to write a variable to visualize the surface. All the data is written into a file, which can be then analyzed to get critical exponents. The Python code used for this task is also incluided.

For example, next graph shows results for correlation lenght vs time, for different system sizes. You can see how it saturates when the correlation gets the size of the system. For the largest size critical exponent has been computed:

![Ballistic](https://github.com/VictorSeven/SurfaceGrowth/blob/master/source/corbal.png "Ballistic")


It is possible also to do the scaling of this graphs:

![Scaling](https://github.com/VictorSeven/SurfaceGrowth/blob/master/source/corbalsca.png "Ballistic Scaling")


You can use this code as an implementation example. 
