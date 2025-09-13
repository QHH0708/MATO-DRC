# MATO-DRC
Multi-Agent Trajectory Optimization with Dynamic Reachability Constraints

**Environment**

Matlab R2023(a) on a Windows operating system.  
We don't use any special features of this version of Matlab; thus, we believe that any version of Matlab will work.

**Dependencies**

1. Multi-Parametric Toolbox 3 (MPT3). [MPT3](https://www.mpt3.org/). For the convenient visualization of reachable sets.
2. OPTI Toolbox. [OPTI](https://github.com/jonathancurrie/OPTI). For solving the nonlinear optimization problem.
3. YALMIP. [YALMIP](https://github.com/yalmip/YALMIP). For solving the linear matrix inequality problem of the ellipsoid convex hull operation.

To install these three dependencies, please visit their official websites, and download their installation .m files. Then, run these .m files in MATLAB, and the installation will be completed automatically.

**File Structure**

All the Matlab code and the results of the numerical experiments are included in the folder "\MATO-DRC".  
All the figures (.fig files) are included in the folder "\MATO-DRC\Figures".  
All data from the numerical experiments are saved as .mat files and included in the folder "\MATO-DRC\Data".

**Numerical Examples**
1. Trajectory Optimization for a Single Quadrotor

   In these numerical experiments, the reachable set of the quadrotor is considered as a fixed ball with a radius of 0.5.
   
   Run "initial_parameter.m". Change "N" for different trajectory segments.
   
   Run "optimization_smoothing_MPCC.m". Directly copy the first line of "optimization_smoothing_MPCC.m" to the Command Window.
   
	 Run "plot_main.m" (for visualization). Directly copy the first line of "plot\_main.m" to the Command Window.  Change the input "sol_index" to a column index of "Solutions" that you want to plot.

	 To plot multiple solutions, set "sol_index" to a row vector where each element represents a solution index.
