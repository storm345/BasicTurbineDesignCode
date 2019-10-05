# BasicWindTurbineDesign
MATLAB code for doing basic aero optimisation and structural predictions for turbine blades 
and then exporting the cross section geometry as a series of .pts files. 

Requires MATLAB parallel toolbox to work properly. Aero predictions come from XFoil so are 2D. 
Centrifugal and bending stresses are computed separately along the span and it is up to the designer 
to make sure the sum of these stresses is within allowable limits. 

This code was originally written for the L2 applications project at Imperial College London but is free for anybody to use 
in design of any turbine component, but be aware of the simplistic nature of the analysis and seek to further 
validate the generated design. 

First code run will be slow as it characterizes, using XFoil, the lift/drag of each blended airfoil section at 
various Re and M but subsequent runs will be significantly faster as this computation will already be done.

The designer must specify the airfoils (in form of files containing the coordinates) and their variation along the span and
must configure the parameters that define the chord distribution through modification of the variables within the file. The chord profile
used is a simple k/r kind.

