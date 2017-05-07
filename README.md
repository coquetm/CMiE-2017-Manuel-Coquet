# CMiE-2017-Manuel-Coquet
Final Project Power Flow and Optimal Power Flow Algorithms

## Policy Appications of Power Systems Modelling

- In this project, I start by building Power Flow model for a small network and comparing the computational efficiency of Newton's method and the optimizing package JuMP in Julia.

- Afterwards, I expand the Power Flow model into an Optimal Power Flow to optimize the operation and planning of a power grid based on minimizing generation costs, and to determine Locational Marginal Prices that generators and load consumers face.

- In addition, I use the Optimal Power Flow model to examine an example of Market Power and the incompatibility of incentives between profit-maximizing for a firm and optimizing the welfare of society.

- Finally, I replicate a Mixed Integer Linear (MILP) approach to solving the transmission expansion planning problem from the paper Transmission Expansion Planning: A Mixed-Integer LP Approach (2003) by Natalia Alguacil, Alexis L. Motto and Antonio J. Conejo.

## 1. AC Power Flow analysis
An alternating current power-flow model is a model used in electrical engineering to analyze power grids. It provides a nonlinear system which describes the energy flow through each transmission line. The goal of a power-flow study is to obtain complete voltages angle and magnitude information for each bus in a power system for specified load and generator real power and voltage conditions.  

Once this information is known, real and reactive power flow on each branch as well as generator reactive power output can be analytically determined. Power-flow or load-flow studies are important for planning future expansion of power systems as well as in determining the best operation of existing systems.

## AC Power Flow Problem description
- Consider a network with 5 buses that forms a cycle (i.e., the lines are (1,2), (2,3), (3,4),(4,5) and (5,1)).
- Assume that Bus 1 is a slack, Bus 2 is PV (generator) and Buses 3-5 are PQ (loads).
- Assume that the resistance and reactance of each line are both equal to 1.

Repeat the following lines 1-3 several times (say 100 times):
- 1: Generate a random vector V such that all voltage magnitudes are somehow close to 1 and all voltage angles are close to 0, and that the phase at the slack bus is zero.
- 2: Based on the random vector of voltages V, compute the loads at the PQ buses, the P and |V| values at the PV buses, and |V| at the slack bus. (constraints/input data needed to run the model)
- 3: Use the measurement data to solve the power flow problem in two ways: (1) Newton's Mehtod, (2) JuMP 
- 4: Declare a success if the obtained solution matches the original random state V. Compute how many times each of the above two methods was successful for different values of voltage angles (sensitivity analysis).

#### NOTE: The problem does not always have a solution, especially at large angles

## 2. AC Optimal Power Flow
It is an expansion of power flow analysis in which power flow are optimized in order to minimize the cost of generation subject to the power flow constraints and other operational constraints, such as generator minimum output constraints, transmission stability and voltage constraints, and limits on switching mechanical equipment.

Equality constraints
- Power balance at each node - power flow equations

Inequality constraints
- Network operating limits (line flows, voltages)
- Limits on control variables

Solving an OPF is necessary to determine the optimal operation and planning of the grid. In this algorithm, I will simulate an optimal power flow model for a 6-bus system and determine the locational marginal prices (LMPs).

## 3. Market Power
- In this exercise, I will show that when a generator can exert market power (manipulate prices) due to cheaper costs -could also be due to size-, they will seek to maximize profit and produce at levels suboptimal for society. 

- When left unregulated, generators with market power will withhold capacity to maximize profit, causing them to produce away from society's optimal point. This will create Deadweight Loss (DWL) that will be paid by ratepayers, justifying regulating generators with market power in markets.

#### The steps to show the impact of market power will be the following:
- Determine the optimal social point by optimizing the power system without restrcitions (Determine optimal power flows, generator profit and systemwide costs)
- Run power flow simulations withholding generation capacity from the generator with market power
- Calculate the new optimal power flow, generator profit, systemwide costs and DWL
- Determine the optimal operation point for the generator based on profit maximizing
- Determine the associated DWL to society based on the generator profit-maximizing

## 4. Transmission expansion Planning (Mixed Integer Programming)
The transmission planning process must identify and support development of transmission infrastructure that is sufficiently robust and can enable competition among wholesale capacity and energy suppliers in energy markets. However, it is a very complex mathematical problem. In this project, I replicate algorithms that try to tackle this problem through linearizations and relaxations.

#### The TEP algorithms I will replicate come from the Paper:
- Transmission Expansion Planning: A Mixed-Integer LP Approach (2003) by Natalia Alguacil, Alexis L. Motto and Antonio J. Conejo

- This paper presents a mixed-integer LP approach to the solution of the long-term transmission expansion planning problem

- In general, this problem is large-scale, mixed-integer, nonlinear, and nonconvex. The authors derive a mixed-integer linear formulation that considers losses and guarantees convergence to optimality using existing optimization software

- The proposed model is applied to Garver’s 6-bus system, the IEEE Reliability Test System, and a realistic Brazilian system. However, I only apply the algorithm to Garver’s 6-bus system.
