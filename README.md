# 3FlexAnalyser.jl
**Julia tool for modelling and analysing flexibility in low voltage distribution networks, with nonlinear three-phase power flow and explicit constraints on voltage unbalance and phase coordination**

The tool enables quantifying the impacts of phase unbalance, voltage unbalance and distributed energy resources (DER) coordination on flexibility services in low voltage distribution networks.
The main idea of aggregated flexibility modelling in the unbalanced setting was originally formulated by Wangwei Kong (The University of Manchester ➔ National Grid). The tool was further developed, tested and published by Andrey Churkin (The University of Manchester ➔ Imperial College London).

At the core of **3FlexAnalyser.jl** lies a nonlinear three-phase optimal power flow (OPF) model adapted from **PowerModelsDistribution.jl** [1]. The mathematical formulation of the three-phase OPF problem is extended by including flexible units, voltage unbalance limits, and phase coordination constraints. To characterise feasible flexibility services available in distribution networks, the concept of aggregated P-Q flexibility areas is applied. That is, the limits of aggregated P-Q flexibility at a selected location are estimated by iteratively solving the three-phase OPF model. Therefore, the tool allows to directly translate voltage unbalance and phase coordination constraints as reductions in the aggregated P-Q flexibility.

A high-level overview of the proposed framework is presented in the diagram below.
The inputs include network data, available flexible resources, load and flexibility unbalances, VUF limit, and phase coordination assumptions.
Then, the limits of flexibility services from DER are estimated using the concept of P-Q flexibility areas. Next, voltage unbalance and phase coordination constraints are imposed and the P-Q flexibility areas are estimated again.
Finally, the aggregated flexibility is compared for cases with and without these constraints and the impacts of phase unbalance and DER coordination are quantified.

<img src="framework_flowchart.png" alt="Framework" width="1000">

An example of the tool's output is shown in the figure below, where the aggregated P-Q flexibility of a distribution network is estimated for phase A.
Each solution of the three-phase OPF problem is displayed by a black circle marker.
The convex (can be concave in some cases) hull of these solutions in the P-Q space represents the set of feasible operating points and thus serves as the estimation of aggregated flexibility of the network.

Note that some of the OPF problems may not converge, be infeasible, or have the status "Converged to a point of local infeasibility". Such solutions are scattered using red circles. These circles are not included in the convex hull and are displayed to inform users of potential issues. A few points with convergence issues may be acceptable for the flexibility analysis. However, if multiple points are red, this indicates problems in the case study or simulation parameters.

<img src="PQ_area_example.png" alt="PQ_area_example" width="500">

<img src="5_bus_scheme.png" alt="5_bus_scheme" width="500">

<img src="221_bus_UK_graph.png" alt="221_bus_UK_graph" width="500">

<img src="results_from_manuscript_5_bus.png" alt="results_from_manuscript_5_bus" width="1000">

<img src="results_from_manuscript_221_bus_UK.png" alt="results_from_manuscript_221_bus_UK" width="1000">



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

[2] To be updated....
