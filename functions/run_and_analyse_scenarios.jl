"""
This script solves OPF problems for scenarios generated in generate_load_and_flexibility_scenarios.jl
It outputs the highest voltage unbalances for each scenario.
"""

using PowerModelsDistribution
using JuMP, Ipopt
using CSV
using DataFrames
using JLD2

cd(dirname(@__FILE__))
println(pwd())



load_path = joinpath(@__DIR__, "..", "results", "load_and_flexibility_scenarios_221_bus_v1.jld2")
@load load_path math_generated_scenarios