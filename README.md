# sourcePy

sourcePy is a python package developed to streamline [HYSPLIT](https://www.arl.noaa.gov/hysplit/) back-trajectory analyses for identifying the source regions of pollutants and other atmospheric components.

In general, back-trajectory source analysis experiments involve taking measurements of atmospheric concentrations and running back-trajectories from that location during the times the observations were taken. Once you have pairings of trajectories and associated observations, an analysis method like the ones included in this package can be performed. sourcePy contains all the tools necessary to run one of these experiments and plot the results assuming a user has NOAA ARL meteorology data and a local installation of HYSPLIT. If a user has already run HYSPLIT trajectories, then just the desired analysis can be done.


### Available Analysis Methods

Potential Source Contribution Function ([PSCF](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2001JD900024))
Concentration-Weighted Trajectories ([CWT](https://www.sciencedirect.com/science/article/pii/S0169809511002948))
Simplified Quantitative Transport Bias Analysis ([SQTBA](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JD030380))
Redistributed Concentration Field ([RCF](https://www.sciencedirect.com/science/article/abs/pii/1352231095003142))


# Modules

#### hysplit.py
The hysplit module is intentended to allow users to quickly and easily generate trajectories using HYSPLIT. It contains the met class and a variety of methods to this end, including the met class, a class that contains descriptive information about available NOAA ARL files. Met objects can be paired with trajectory request objects to generate HYSPLIT's CONTROL and SETUP files, and then HYSPLIT can be run directly within the python kernel.

NOTE: This module currently only works on Windows machines. Non-windows users can still use the other modules.

Coming soon: full Unix support, HYSPLIT concentrations

#### source.py
The primary calculations and handling of trajectories and concentration observations are held in the source module. There are three primary classes in source: The measurement class, trajectory class, and trajectory request class. 

**Measurements** can be defined manually or read in from data files, and contain information about an observation of some atmospheric aerosol (location, timing, value, etc). Measurements are necessary for the source analyses implemented here and can also be useful for determining the location and timing of trajectories.

**Trajectory requests** are objects that contain information necessary for use in the hysplit module, like the start and end times and the start location of the trajectory.

**Trajectories** are objects that contain information from a HYSPLIT trajectory file along with associated measurement objects.

If the user has read in trajectories, paired them with measurements, and defined a few other parameters like critical values for PSCF and the analysis grid, then any of the available analysis methods can be used to generate 2D source analysis fields.

Coming soon: Faster/alternative RCF calculation methods, concentration analyses

#### plotSource.py
plotSource uses matplotlib and cartopy to plot data about arrays of trajectories and measurements, generate basic geographic maps, and plot the 2D source analysis fields. It is not a comprehensive plotting module, and is only intended to get new users off on the right foot when they look to start plotting their data.
