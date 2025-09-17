# MATO-DRC
Multi-Agent Trajectory Optimization with Dynamic Reachability Constraints

**Environment**

Matlab R2023(a) on a Windows operating system.  
We don't use any special features of this version of Matlab; thus, we believe that any version of Matlab will work.

**Dependencies**

1. Multi-Parametric Toolbox 3 (MPT3). [MPT3](https://www.mpt3.org/). For the convenient visualization of reachable sets.
2. OPTI Toolbox. [OPTI](https://github.com/jonathancurrie/OPTI). For solving the nonlinear optimization problem.
3. YALMIP. [YALMIP](https://github.com/yalmip/YALMIP). For solving the linear matrix inequality problem of the ellipsoid convex hull operation.

To install these three dependencies, please visit their official websites.

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

   <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A1_multi_N.png" alt="Trajectory Optimization for 1 agents" width="500" align=center />

2. Trajectory Optimization for Multi-Quadrotor with Ellipsoid Reachable Sets

	Run "initial_parameter_MultiAgent.m"
	
	Run "opti_smoothing_MPCC_MultiAgent.m". Directly copy the first line of "opti_smoothing_MPCC_MultiAgent.m" to the Command Window.
	
	Run "plot_main_MultiAgent.m" (for visualization). Directly copy the first line of "plot_main_MultiAgent.m" to the Command Window. Change the input "sol_index" to a column index of "Solutions" that you want to plot.
	
	Change the input "zono_or_ellipsoid" to the number "0" for ellipsoid reachable set case.
	
	To visualize the trajectories, activate the code at line 180 of "test_b_to_trajectory_MultiAgent.m".

<table>
  <tr>
    <td align="center">
      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A5_N4_with_Bezier.png" alt="Trajectory Optimization for 5 agents with ellipsoid reachable sets" width="500" align=center />
      <br>
      <sub><b>Trajectories of 5 agents and their Bezier polygons</b></sub>
    </td>
    <td align="center">
      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A5_N4_with_RS.png" alt="Trajectory Optimization for 5 agents with ellipsoid reachable sets" width="500" align=center />
      <br>
      <sub><b>Trajectories of 5 agents and the ellipsoid reachable sets</b></sub>
    </td>
  </tr>
</table>



 <table>
  <tr>
    <td align="center">
      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A5_N4_elipRS_r12.png" alt="Trajectory Optimization for 5 agents with ellipsoid reachable sets" width="500" align=center />
      <br>
      <sub><b>r1-r2 view</b></sub>
    </td>
    <td align="center">
      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A5_N4_elipRS_r13.png" alt="Trajectory Optimization for 5 agents with ellipsoid reachable sets" width="500" align=center />
      <br>
      <sub><b>r1-r3 view</b></sub>
    </td>
  </tr>
</table>

 3. Trajectory Optimization for Multi-Quadrotor with Zonotope Reachable Sets

    Run "initial_parameter_MultiAgent_zonoRS.m"
	
	Run "opti_smoothing_MPCC_MultiAgent_zonoRS.m". Directly copy the first line of "opti_smoothing_MPCC_MultiAgent.m" to the Command Window.
	
	Run "plot_main_MultiAgent.m" (for visualization). Directly copy the first line of "plot_main_MultiAgent.m" to the Command Window. Change the input "sol_index" to a column index of "Solutions" that you want to plot.
	
	Change the input "zono_or_ellipsoid" to the number "1" for zonotope reachable set case.

<table>
  <tr>
    <td align="center">
      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A4_N4_with_RS_zono.png" alt="Trajectory Optimization for 5 agents with zonotope reachable sets" width="500" align=center />
      <br>
      <sub><b>Trajectories of 5 agents and the zonotope reachable sets</b></sub>
    </td>
  </tr>
</table>

 <table>
  <tr>
    <td align="center">
      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A5_N4_zonoRS_r12.png" alt="Trajectory Optimization for 5 agents with zonotope reachable sets" width="500" align=center />
      <br>
      <sub><b>r1-r2 view</b></sub>
    </td>
    <td align="center">
      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A5_N4_zonoRS_r13.png" alt="Trajectory Optimization for 5 agents with ellipsoid reachable sets" width="500" align=center />
      <br>
      <sub><b>r1-r3 view</b></sub>
    </td>
  </tr>
</table>

4. Examples for Trajectory Optimization of Larger-Scale Agent systems (The figures can be find in folder \MATO-DRC\Figures, corresponding code has not been published)

   Example with 12 agents
	<table>
	  <tr>
	    <td align="center">
	      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A12_N4_with_Bezier.png" alt="Trajectory Optimization for 12 agents" width="500" align=center />
	      <br>
	      <sub><b>Trajectory optimization results of 12 agents</b></sub>
	    </td>
	    <td align="center">
	      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A12_N4_with_Bezier_xy.png" alt="Trajectory Optimization for 12 agents" width="500" align=center />
	      <br>
	      <sub><b>x-y view</b></sub>
	    </td>
	  </tr>
	</table>

 	<table>
	  <tr>
	    <td align="center">
	      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A12_N4_with_Bezier_xz.png" alt="Trajectory Optimization for 12 agents" width="500" align=center />
	      <br>
	      <sub><b>x-z view</b></sub>
	    </td>
	    <td align="center">
	      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A12_N4_with_RS.png" alt="Trajectory Optimization for 12 agents" width="500" align=center />
	      <br>
	      <sub><b>Reachable sets</b></sub>
	    </td>
	  </tr>
	</table>
 
   Example with 16 agents
	<table>
	  <tr>
	    <td align="center">
	      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A16_N4_with_Bezier.png" alt="Trajectory Optimization for 12 agents" width="500" align=center />
	      <br>
	      <sub><b>Trajectory optimization results of 12 agents</b></sub>
	    </td>
	    <td align="center">
	      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A16_N4_with_Bezier_xy.png" alt="Trajectory Optimization for 12 agents" width="500" align=center />
	      <br>
	      <sub><b>x-y view</b></sub>
	    </td>
	  </tr>
	</table>

 	<table>
	  <tr>
	    <td align="center">
	      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A16_N4_with_Bezier_xz.png" alt="Trajectory Optimization for 12 agents" width="500" align=center />
	      <br>
	      <sub><b>x-z view</b></sub>
	    </td>
	    <td align="center">
	      <img src="https://github.com/QHH0708/MATO-DRC/blob/main/MATO-DRC/Figures/A16_N4_with_RS.png" alt="Trajectory Optimization for 12 agents" width="500" align=center />
	      <br>
	      <sub><b>Reachable sets</b></sub>
	    </td>
	  </tr>
	</table>
