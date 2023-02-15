# MATLAB_modeling_two_phase_flow (UPDATING!)

In this project, I formulate an IMPES (Implicit Pressure, Explicit Saturation) method to simulate two-phase fluid displacement experiments in porous media. Consider an experimental setup consisting of a pump, a core sample and a pressure valve. The core is initially saturated with oil and connate water i.e. So_init = 1 - Swc. The pump is connected to the inlet of the core sample and starts injecting pure water into the core with a constant rate (measured at Standard Condition). At the outlet a pressure valve is connected to the core that maintains a constant outlet pressure during the experiment.


Assumptions:

· The core is horizontal and gravity forces are neglected.

· The core is homogenous and isotropic.

· Capillary pressure is neglected.

· At the beginning of the experiment the entire core is pressurized to an initial pressure.

· Fluids are slightly compressible.

· The rock fluid properties are given.


Here are the data input for the model and some key results.

**Data input deck**

1. Relative permeabilities and capillary pressure:

2. Fluid Densities and viscosities:

· Fluids are slightly compressible with constant compressibility:

· Fluids have constant viscosity:

3. Rock properties:

· Rock is homogenous, isotropic and incompressible:

4. Injection and production mode and initial condition:

· Constant injection at the inlet:

· Production with constant pressure at the outlet:

· Initial condition:

   S</sup>o_ini<sup> = 1 - S</sup>wc<sup> , P</sup>ini<sup> = 100 [bar]


**Modeling output and sensitivities of the mobility ratio and time step**

The mobility ratio is defined as M = (krw.mu_o)/(kro.mu_w). Increasing the mobility ratio is equivalent to decreasing the ratio of viscosities (mu_w/mu_o). The results suggest that increasing the mobility ratio: reduces the time to water break through and lowers the efficiency of production (by increasing the time to produce all viable oil). However, increasing the mobility ratio does reduce the injection pressure necessary to support a constant injection rate.

**1. mu_o = mu_w**

![testAnimated_muy_o=muy_w_stable2](https://user-images.githubusercontent.com/86640902/219100477-30e5523c-839d-4723-962c-1a01ee403478.gif)


![Output_mu_w=mu_o_stable](https://user-images.githubusercontent.com/86640902/219101223-1076c638-db15-4844-aaf9-38801aa84b69.jpg)


**2. mu_o = 10mu_w**

![testAnimated_muy_o=10muy_w_stable2](https://user-images.githubusercontent.com/86640902/219100542-a2bb8634-b695-44d9-959a-bcf49ab8fc25.gif)

![Output_10mu_w=mu_o_stable](https://user-images.githubusercontent.com/86640902/219100790-de0c854f-1f95-4b78-840b-fd5ccbead061.jpg)

**3. Unstable solution with mu_o = mu_w**

Although I used an implicit method (IMPES) to solve for pressure, stability is not guaranteed due to the embedded assumption of constant saturation over a timestep. Reducing the grid size and time step both tend to increase the accuracy of results, as long as stability
requirements are met. However, reducing the time step or grid size increase the time to run the simulations. In this sensityvity, I increased the timestep size from Nt = 400 (above stable solutions) to Nt = 100 to illustrate how computational overhead affect the solution stability.   

![testAnimated_muy_o=muy_w](https://user-images.githubusercontent.com/86640902/219053413-d4b16f24-548b-4912-b58d-d7b184ec41a2.gif)

![Output_mu_w=mu_o_unstable](https://user-images.githubusercontent.com/86640902/219101361-1d59c26b-1282-47b2-ba3e-f49efd3ceb5b.jpg)


