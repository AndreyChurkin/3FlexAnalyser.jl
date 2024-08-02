# 3FlexAnalyser.jl
**Julia tool for modelling and analysing flexibility in low voltage distribution networks, with nonlinear three-phase power flow and explicit constraints on voltage unbalance and phase coordination**

The tool enables quantifying the impacts of phase unbalance, voltage unbalance and distributed energy resources (DER) coordination on flexibility services in low voltage distribution networks.
The main idea of aggregated flexibility modelling in the unbalanced setting was originally formulated by Wangwei Kong (The University of Manchester ➔ National Grid). The tool was further developed, tested and published by Andrey Churkin (The University of Manchester ➔ Imperial College London).

At the core of **3FlexAnalyser.jl** lies a nonlinear three-phase optimal power flow (OPF) model adapted from **PowerModelsDistribution.jl** [1]. The mathematical formulation of the three-phase OPF problem is extended by including flexible units, voltage unbalance limits, and phase coordination constraints. To characterise feasible flexibility services available in distribution networks, the concept of aggregated P-Q flexibility areas is applied. That is, the limits of aggregated P-Q flexibility at a selected location are estimated by iteratively solving the three-phase OPF model. Therefore, the tool allows to directly translate voltage unbalance and phase coordination constraints as reductions in the aggregated P-Q flexibility. Please see more details and examples in the accompanying manuscript [2].

A high-level overview of the proposed framework is presented in the diagram below.
The inputs include network data, available flexible resources, load and flexibility unbalances, voltage unbalance factor (VUF) limits, and phase coordination assumptions.
Then, the limits of flexibility services from DER are estimated using the concept of P-Q flexibility areas. Next, voltage unbalance and phase coordination constraints are imposed and the P-Q flexibility areas are estimated again.
Finally, the aggregated flexibility is compared for cases with and without these constraints and the impacts of phase unbalance and DER coordination are quantified.

<img src="framework_flowchart.png" alt="Framework" width="1000">

An example of the tool's output is shown in the figure below, where the aggregated P-Q flexibility of a distribution network is estimated for phase A.
Each solution of the three-phase OPF problem is displayed by a black circle marker.
The convex (can be concave in some cases) hull of these solutions in the P-Q space represents the set of feasible operating points and thus serves as the estimation of aggregated flexibility of the network.

Note that some of the OPF problems may not converge, be infeasible, or have the status "Converged to a point of local infeasibility". Such solutions are scattered using red circles. These circles are not included in the convex hull and are displayed to inform users of potential issues. A few points with convergence issues may be acceptable for the flexibility analysis. However, if multiple points are red, this indicates problems in the case study or simulation parameters.

<img src="PQ_area_example.png" alt="PQ_area_example" width="500">

A single run of the tool produces one P-Q flexibility area, saved in `results/` together with the simulations data in JLD format. To quantify the impact of voltage unbalance limits, phase coordination constraints or other parameters, the tool has to be run multiple times and its outputs should be compared.
For example, in [2], various combinations of voltage unbalance and phase coordination constraints have been examined. The outputs of the tool from multiple runs have been combined and presented in a single figure to make the analysis easy for readers.

There are two cases included and tested in this tool, as shown below.

Illustrative 5-bus system with 3 flexible units:

<img src="5_bus_scheme.png" alt="5_bus_scheme" width="500">

Real 221-bus low voltage distribution network in the UK with 12 flexible units:

<img src="221_bus_UK_graph.png" alt="221_bus_UK_graph" width="500">

Extensive simulations performed for these cases demonstrate that flexibility services in distribution networks can be significantly constrained due to inherent load unbalances, voltage unbalance limits, and lack of coordination between DER connected to different phases.
Specifically, the aggregated P-Q flexibility areas reduce when phase coordination constraints and voltage limits are imposed.
It is found that the worst conditions for providing flexibility services include the lack of coordination between flexible units connected to different phases and tight voltage unbalance constraints.
In such cases, flexible units cannot effectively manage voltage unbalance across different phases and locations, which results in the infeasibility of services provision.
The figures below show the impact of phase coordination and voltage unbalance limits for (a) phase A, (b) phase B, and (c) phase C.

Example of the results and analysis for the illustrative 5-bus system:

<img src="results_from_manuscript_5_bus.png" alt="results_from_manuscript_5_bus" width="1000">

Example of the results and analysis for the 221-bus distribution network in the UK:

<img src="results_from_manuscript_221_bus_UK.png" alt="results_from_manuscript_221_bus_UK" width="1000">

### RUNNING THE TOOL:

To run the tool, execute `functions/3FlexAnalyser.jl`.
This is the main code for calculating P-Q flexibility areas, plotting them, and imposing voltage unbalance and phase coordination constraints.
The code reads data about the network and flexibility from `/cases` in the OpenDSS format (.dss).
Yet, many parameters of the simulations have to be specified in this code to produce correct simulations.
For example, it is necessary to select the case study (.dss file), P-Q limits of flexible units, phase (A, B, or C) for which the P-Q flexibility area will be maximised, location of flexibility aggregation, number of simulations (points used to estimate P-Q flexibility areas), voltage unbalance constraints, phase coordination constraints, some plotting parameters, etc.

This code is not perfect and may be further improved. Nonetheless, as of August 2024, the code works reliably and reasonably fast for the two included case studies.

The tool has been tested in Julia v1.10.4 (2024-06-04) with the following packages:
- CSV v0.10.14
- ConcaveHull v1.1.0
- DataFrames v1.6.1
- Ipopt v1.6.3
- JLD v0.13.5
- JuMP v1.22.2
- LazySets v2.14.1
- Plots v1.40.4
- Polyhedra v0.7.8
- PowerModelsDistribution v0.15.2
- StatsPlots v0.15.7
- Suppressor v0.2.7

### REFERENCES:
[1] Fobes, David M., Sander Claeys, Frederik Geth, and Carleton Coffrin, "PowerModelsDistribution.jl: An open-source framework for exploring distribution power flow formulations," Electric Power Systems Research, vol. 189, 2020.

[2] To be updated.......
