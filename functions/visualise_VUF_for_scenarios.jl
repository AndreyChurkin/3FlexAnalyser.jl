"""
This script can be used to visualise and compare VUFs for scenarios computed in run_and_analyse_scenarios.jl
"""

using DataFrames
using JLD2
using Plots, Plots.PlotMeasures

cd(dirname(@__FILE__))
println(pwd())



fz = 16
plt_max_vuf = plot(
    # title = "VUF duration curve",
    xlabel = "Number of scenarios",
    ylabel = "Maximum VUF, %",

    size = (1000,1000), # width and height of the whole plot (in px)
    # aspect_ratio = :equal,
    ylim = (0.0, 2.3),
    # xlim = (0,1000),

    xtickfontsize=fz, ytickfontsize=fz,
    fontfamily = "Courier", 
    titlefontsize = fz,
    xguidefontsize = fz,
    yguidefontsize = fz,
    legendfontsize = fz-2,

    # legend = false,
    
    framestyle = :box,
    margin = 10mm,

    minorgrid = :true,
)

load_path = joinpath(@__DIR__, "..", "results", "VUFs_for_1000_scenarios_mode_no_flexibility.jld2")
@load load_path VUF_per_scenario

plot!(plt_max_vuf,
    collect(1:length(VUF_per_scenario)),
    VUF_per_scenario*100, 
    label = "No flexibility",
    color = :grey,
    w = 5
)


load_path = joinpath(@__DIR__, "..", "results", "VUFs_for_1000_scenarios_mode_random_and_fixed.jld2")
@load load_path VUF_per_scenario

plot!(plt_max_vuf,
    collect(1:length(VUF_per_scenario)),
    VUF_per_scenario*100, 
    label = "Random fixed flexibility",
    color = :black,
    w = 5
)


load_path = joinpath(@__DIR__, "..", "results", "VUFs_for_1000_scenarios_mode_random_min_PQ_in_one_phase1.jld2")
@load load_path VUF_per_scenario

plot!(plt_max_vuf,
    collect(1:length(VUF_per_scenario)),
    VUF_per_scenario*100, 
    label = "Phase A random consumption",
    # color = :darkred,
    color = "#B22222",
    # color = :navy,
    w = 5,
    # linestyle=:dash
)



display(plt_max_vuf)

plt_name = "VUF_visualisation_and_comparison_for_scenarios"
savefig("../results/"*plt_name*".png")
savefig("../results/"*plt_name*".pdf")
savefig("../results/"*plt_name*".svg")