# Contact-Aware-Controller-Design-for-Complementarity-Systems
ICRA2020

arXiv link: https://arxiv.org/abs/1909.11221

Video: https://www.youtube.com/watch?v=WS4nMXtCxcQ

## Dependencies
The linear complementarity problems (LCPs) generated by this library are solved using [http://pages.cs.wisc.edu/~ferris/path.html](PATH). 

The bilinear matrix inequalities are solved using http://www.penopt.com/penbmi.html (PenBMI). Please contact the creators in order to obtain the license.

Optimization problems are formulated using https://yalmip.github.io/ (YALMIP).

`pathlcp`, `PenBMI`, `yalmip` will need to be in the MATLAB path for the examples to run.

## Functionality
The library can be used to design contact-aware controllers for linear complementarity systems (requires YALMIP and PenBMI). The code can be used to design controllers for any linear complementarity system models as long as the P-matrix assumption holds. It is important to note that the set related to `\bar{\lambda}` needs to be generated specifically wrt the system at hand. Make sure the S-procedure terms related to that set is correct for your model before running the code to design a controller.

The designed controller can be tested on the linear complementarity system (recommended as a sanity check). The code can be used to evaluate the dynamics of any linear complementarity model as long as the P-matrix assumption holds (requires PATH).

The designed controller can be tested on the nonlinear complementarity system model (requires PATH).

## Examples
`acrobot`: Controller design and its implementation on an acrobot with soft joint limits

`cartpole`: Controller design and its implementation on a cartpole with soft walls

`partial_feedback`: Controller design and its implementation on a model with a cartpole and two carts, where the cart in the middle is not observed


