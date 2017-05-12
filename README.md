# CMiE-2017-Manuel-Coquet
Final Project Computational Methods in Economics

## Policy Appications of Power Systems Modelling

### Power Flow Algorithms:


- In this project, I start by building Power Flow model for a small network and comparing the computational efficiency of Newton's method and the optimizing package JuMP in Julia.


- Afterwards, I expand the Power Flow model into an Optimal Power Flow to optimize the operation and planning of a power grid based on minimizing generation costs, and to determine Locational Marginal Prices that generators and load consumers face.


### Market power regulation:


- I use an Optimal Power Flow model to examine an example of Market Power and the need for regulation due to the incompatibility of incentives between profit-maximizing for a firm and optimizing the welfare of society. 


- This is done by determining the social optimal (cost-minimization) grid power flow, and then showing that a generation can withhold generation to manipulate prices to their own benefit. 


- At the profit maximization for the generator, I show that there is a substantial deadweight loss to society.


### Transmission expansion planning:


- First, I replicate a Mixed Integer Linear (MILP) approach to solving the transmission expansion planning problem from the paper Transmission Expansion Planning: A Mixed-Integer LP Approach (2003) by Natalia Alguacil, Alexis L. Motto and Antonio J. Conejo.


- I apply a DC algorithm without considering losses and another DC LP algorithm linearizing losses to the Garver 6-bus system, and I show that my results are the same as the paper.

    
- I reach the same conclusions as the author; using a DC model without considering losses to project transmission expansion planning leads to underinvestment - however, computing technology now allows us to approach better solutions using more elaborate algorithms.


- Nonetheless, I also found that considering losses becomes computationally intensive, it takes 25 times more time to solve the model with losses. For larger systems, this may create concerns as computing time rises exponentially.
  
  
- Since the DC Algorithm does not guarantee AC Feasibility, I also built a relaxed AC Transmission Expansion Planning model based on the paper Transmission Expansion Planning Using an AC Model: Formulations and Possible Relaxations (2012) by Zhang et al. 


- For the AC TEP algorithm, I ran the algorithm in KNITRO in Julia locally - for some reason there is a bug in Jupyter notebook -, but I pasted the code and results (local code is also found in github repository).


- AC algorithms do not guarantee global optimum since the problem is non convex. Nonetheless, I show that it is possible to find AC feasible local solutions using a multistart approach at 500 different initial points. The best local solution that I found requires twice as much investment as the DC solution.
