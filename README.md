# 3FlexAnalyser.jl
**Julia tool for modelling and analysing flexibility in low voltage distribution networks, with nonlinear three-phase power flow and explicit constraints on voltage unbalance and phase coordination**

... don't forget to mention PowerModelsDistribution.jl ...

<img src="framework_flowchart.png" alt="Framework" width="1000">

<img src="5_bus_scheme.png" alt="Framework" width="500">

<img src="221_bus_UK_graph.png" alt="Framework" width="500">


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