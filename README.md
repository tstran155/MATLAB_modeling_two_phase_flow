# MATLAB_modeling_two_phase_flow (UPDATING!)

Here are the data input for the model and some key results.

**Data input deck**


**mu_o = mu_w**

![testAnimated_muy_o=muy_w_stable2](https://user-images.githubusercontent.com/86640902/219100477-30e5523c-839d-4723-962c-1a01ee403478.gif)


![Output_mu_w=mu_o_stable](https://user-images.githubusercontent.com/86640902/219101223-1076c638-db15-4844-aaf9-38801aa84b69.jpg)


**mu_o = 10mu_w**

![testAnimated_muy_o=10muy_w_stable2](https://user-images.githubusercontent.com/86640902/219100542-a2bb8634-b695-44d9-959a-bcf49ab8fc25.gif)

![Output_10mu_w=mu_o_stable](https://user-images.githubusercontent.com/86640902/219100790-de0c854f-1f95-4b78-840b-fd5ccbead061.jpg)

**Unstable solution with mu_o = mu_w**

Although I used an implicit method to solve for pressure, stability is not guaranteed due to the embedded assumption of constant saturation over a timestep. Reducing the grid size and time step both tend to increase the accuracy of results, as long as stability
requirements are met. However, reducing the time step or grid size increase the time to run the simulations.

![testAnimated_muy_o=muy_w](https://user-images.githubusercontent.com/86640902/219053413-d4b16f24-548b-4912-b58d-d7b184ec41a2.gif)

![Output_mu_w=mu_o_unstable](https://user-images.githubusercontent.com/86640902/219101361-1d59c26b-1282-47b2-ba3e-f49efd3ceb5b.jpg)


