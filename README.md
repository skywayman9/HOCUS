# HOCUS



HOCUS (High-Order Central Upwind Scheme), which combines the MP scheme and a linear-compact scheme using the BVD principle.

Novel hybrid nonlinear explicit-compact scheme for shock- capturing based on a boundary variation diminishing (BVD) reconstruction. In our approach, we combine a non-dissipative sixth-order central compact interpolation and a fifth-order monotonicity preserving scheme (MP5) through the BVD algorithm. For a smooth solution, the BVD reconstruction chooses the highest order possible interpolation, which is central, i.e. non-dissipative in the current approach and for the discontinuities, the algorithm selects the monotone scheme. This method provides an alternative to the existing adaptive upwind-central schemes in the literature. 

Example 4.9, 4.10 and 4.11 in the following paper, published in Journal of Computational Physcics, are solved in this code, github_hocus.f90.


https://doi.org/10.1016/j.jcp.2020.110067
https://arxiv.org/pdf/2012.09905


It doesn't exactly "reproduce" the result for all the cases. I haven't looked at this code in an year as I developed a new approach. 

https://arxiv.org/pdf/2205.01034
https://arxiv.org/pdf/2106.01738

This algorithm can also reuse gradients for viscous flux computations (du/dx =  u_i+1/2 - u_i-1/2, where u_i+1/2 is the central interpolation scheme, Eqn 28 in the paper) as in my current work. 

As always ---- If you find them helpful, use them, but if they don't work for your case, please don't blame me.


Sainath
