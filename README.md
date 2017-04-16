# CMiE-2017-Manuel-Coquet
Final Project Power Flow Algorithms

## AC Power Flow analysis
An alternating current power-flow model is a model used in electrical engineering to analyze power grids. It provides a nonlinear system which describes the energy flow through each transmission line. The goal of a power-flow study is to obtain complete voltages angle and magnitude information for each bus in a power system for specified load and generator real power and voltage conditions.  

Once this information is known, real and reactive power flow on each branch as well as generator reactive power output can be analytically determined. Power-flow or load-flow studies are important for planning future expansion of power systems as well as in determining the best operation of existing systems.

## Problem description
- Consider a network with 5 buses that forms a cycle (i.e., the lines are (1,2), (2,3), (3,4),(4,5) and (5,1)).
- Assume that Bus 1 is a slack, Bus 2 is PV (generator) and Buses 3-5 are PQ (loads).
- Assume that the resistance and reactance of each line are both equal to 1.

Repeat the following lines 1-3 several times (say 100 times):
- 1: Generate a random vector V such that all voltage magnitudes are somehow close to 1 and all voltage angles are close to 0, and that the phase at the slack bus is zero.
- 2: Based on the random vector of voltages V, compute the loads at the PQ buses, the P and |V| values at the PV buses, and |V| at the slack bus. (constraints/input data needed to run the model)
- 3: Use the measurement data to solve the power flow problem in two ways: (1) Newton's Mehtod, (2) JuMP 
- 4: Declare a success if the obtained solution matches the original random state V. Compute how many times each of the above two methods was successful for different values of voltage angles (sensitivity analysis).

#### NOTE: The problem does not always have a solution, especially at large angles
