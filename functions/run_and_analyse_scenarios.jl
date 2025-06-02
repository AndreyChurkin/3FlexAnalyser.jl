"""
This script solves OPF problems for scenarios generated in generate_load_and_flexibility_scenarios.jl
The outputs are the distribution of loads across phases and scenarios, and the highest voltage unbalance for each scenario.
"""

using PowerModelsDistribution
using JuMP, Ipopt
using CSV
using DataFrames
using JLD2
using Distributions


cd(dirname(@__FILE__))
println(pwd())

include("../functions/calculate_VUF_a_posteriori.jl")



# # Select a case to analyse:

# eng = parse_file("../cases/5_bus_case_illustrative/LV_balanced_flex_balanced.dss")
# eng = parse_file("../cases/5_bus_case_illustrative/LV_unbalanced_flex_balanced.dss")
# eng = parse_file("../cases/5_bus_case_illustrative/LV_unbalanced_flex_unbalanced.dss")

eng = parse_file("../cases/221_bus_real_UK_case/Master_221_bus_UK.dss")

eng["settings"]["sbase_default"] = 1 # if = 1, pm.model will be in kW
# eng["settings"]["voltage_scale_factor"] = 1
eng["settings"]["power_scale_factor"] = 1000

solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>1,"tol"=>1e-9)

math = transform_data_model(eng)



# Load the generated scenarios:

load_path = joinpath(@__DIR__, "..", "results", "load_and_flexibility_1000_scenarios_221_bus.jld2")
@load load_path math_generated_scenarios



# # Select parameters for analysing the scenarios:

# Select "true" to run the OPF model for each of the generated scenarios:
solve_OPF_per_scenario = true

# Assume how flexible units will be operating:
# flexibility_mode = "no_flexibility"
# flexibility_mode = "random_and_fixed"
# flexibility_mode = "minimum_P_and_Q"
# flexibility_mode = "maximum_P_and_Q"
# flexibility_mode = "maximum_PQ_in_one_phase"
# flexibility_mode = "minimum_PQ_in_one_phase"
flexibility_mode = "random_min_PQ_in_one_phase"




gen_ε = 1e-4 # <-- gap

# # Set the P-Q limits for flexible generators (in kW and kVAr):

# # (5-bus case)
# gen_lim_Pmax = 8
# gen_lim_Pmin = -8
# gen_lim_Qmax = 8
# gen_lim_Qmin = -8

# # (221-bus UK case)
gen_lim_Pmax = 5.0
gen_lim_Pmin = -5.0
gen_lim_Qmax = 5.0
gen_lim_Qmin = -5.0



# Analyse the scenarios one by one:

global total_load_kW_perphase = Array{Float64}(undef, 0, 3)
global total_load_kVAr_perphase = Array{Float64}(undef, 0, 3)

global nonconverging_solutions = false

global max_VUF_per_scenario = Array{Float64}(undef, 0, 1)

@time begin
for sc = 1:length(math_generated_scenarios)
    total_load_kW_perphase_sc = [0.0, 0.0, 0.0]'
    total_load_kVAr_perphase_sc = [0.0, 0.0, 0.0]'

    # Replacing values in the math model with the values from the scenarios:
    math["load"] = math_generated_scenarios["scenario_$(sc)"]["load"]
    math["gen"] = math_generated_scenarios["scenario_$(sc)"]["gen"]

    for load_i = 1:length(math["load"])
        n_phases = length(math["load"][string(load_i)]["connections"])
        for phase_i = 1:n_phases
            phase_connection = math["load"][string(load_i)]["connections"][phase_i]
            total_load_kW_perphase_sc[phase_connection] += math["load"][string(load_i)]["pd"][phase_i]
            total_load_kVAr_perphase_sc[phase_connection] += math["load"][string(load_i)]["qd"][phase_i]
        end
    end

    global total_load_kW_perphase = vcat(total_load_kW_perphase, total_load_kW_perphase_sc)
    global total_load_kVAr_perphase = vcat(total_load_kVAr_perphase, total_load_kVAr_perphase_sc)

    if solve_OPF_per_scenario == true
        println("Solving OPF for scenario #",sc)

        if flexibility_mode == "no_flexibility"
            for gen_i in keys(math["gen"])
                if math["gen"][gen_i]["name"] != "_virtual_gen.voltage_source.source"
                    delete!(math["gen"], gen_i)
                end
            end

        elseif flexibility_mode == "random_and_fixed"
            println("Note: the generated random P-Q limits for flexible are fixed, using ε-bounded constraints to avoid solver infeasibility")

        elseif flexibility_mode == "minimum_P_and_Q"
            for gen_i in keys(math["gen"])
                if math["gen"][gen_i]["name"] != "_virtual_gen.voltage_source.source"
                    math["gen"][gen_i]["pg"] = [gen_lim_Pmin]
                    math["gen"][gen_i]["pmin"] = [gen_lim_Pmin - gen_ε]
                    math["gen"][gen_i]["pmax"] = [gen_lim_Pmin + gen_ε]
                    math["gen"][gen_i]["qg"] = [gen_lim_Qmin]
                    math["gen"][gen_i]["qmin"] = [gen_lim_Qmin - gen_ε]
                    math["gen"][gen_i]["qmax"] = [gen_lim_Qmin + gen_ε]
                end
            end

        elseif flexibility_mode == "maximum_P_and_Q"
            for gen_i in keys(math["gen"])
                if math["gen"][gen_i]["name"] != "_virtual_gen.voltage_source.source"
                    math["gen"][gen_i]["pg"] = [gen_lim_Pmax]
                    math["gen"][gen_i]["pmin"] = [gen_lim_Pmax - gen_ε]
                    math["gen"][gen_i]["pmax"] = [gen_lim_Pmax + gen_ε]
                    math["gen"][gen_i]["qg"] = [gen_lim_Qmax]
                    math["gen"][gen_i]["qmin"] = [gen_lim_Qmax - gen_ε]
                    math["gen"][gen_i]["qmax"] = [gen_lim_Qmax + gen_ε]
                end
            end

        elseif flexibility_mode == "maximum_PQ_in_one_phase"
            global phase_to_use = 1 # <-- only this phase will provide flexibility
            for gen_i in keys(math["gen"])
                if math["gen"][gen_i]["name"] != "_virtual_gen.voltage_source.source"
                    if math["gen"][gen_i]["connections"][1] == phase_to_use
                        math["gen"][gen_i]["pg"] = [gen_lim_Pmax]
                        math["gen"][gen_i]["pmin"] = [gen_lim_Pmax - gen_ε]
                        math["gen"][gen_i]["pmax"] = [gen_lim_Pmax + gen_ε]
                        math["gen"][gen_i]["qg"] = [gen_lim_Qmax]
                        math["gen"][gen_i]["qmin"] = [gen_lim_Qmax - gen_ε]
                        math["gen"][gen_i]["qmax"] = [gen_lim_Qmax + gen_ε]
                    else
                        delete!(math["gen"], gen_i)
                    end
                end
            end

        elseif flexibility_mode == "minimum_PQ_in_one_phase"
            global phase_to_use = 1 # <-- only this phase will provide flexibility
            for gen_i in keys(math["gen"])
                if math["gen"][gen_i]["name"] != "_virtual_gen.voltage_source.source"
                    if math["gen"][gen_i]["connections"][1] == phase_to_use
                        math["gen"][gen_i]["pg"] = [gen_lim_Pmin]
                        math["gen"][gen_i]["pmin"] = [gen_lim_Pmin - gen_ε]
                        math["gen"][gen_i]["pmax"] = [gen_lim_Pmin + gen_ε]
                        math["gen"][gen_i]["qg"] = [gen_lim_Qmin]
                        math["gen"][gen_i]["qmin"] = [gen_lim_Qmin - gen_ε]
                        math["gen"][gen_i]["qmax"] = [gen_lim_Qmin + gen_ε]
                    else
                        delete!(math["gen"], gen_i)
                    end
                end
            end

        elseif flexibility_mode == "random_min_PQ_in_one_phase"
            for gen_i in keys(math["gen"])
                if math["gen"][gen_i]["name"] != "_virtual_gen.voltage_source.source"
                    if math["gen"][gen_i]["connections"][1] == phase_to_use
                        flex_multiplier = rand(Uniform(0,1))

                        math["gen"][gen_i]["pg"] = [gen_lim_Pmin*flex_multiplier]
                        math["gen"][gen_i]["pmin"] = [gen_lim_Pmin*flex_multiplier - gen_ε]
                        math["gen"][gen_i]["pmax"] = [gen_lim_Pmin*flex_multiplier + gen_ε]
                        math["gen"][gen_i]["qg"] = [gen_lim_Qmin*flex_multiplier]
                        math["gen"][gen_i]["qmin"] = [gen_lim_Qmin*flex_multiplier - gen_ε]
                        math["gen"][gen_i]["qmax"] = [gen_lim_Qmin*flex_multiplier + gen_ε]
                    else
                        delete!(math["gen"], gen_i)
                    end
                end
            end

        end

        pm = instantiate_mc_model(math, ACPUPowerModel, build_mc_opf)

        solution_opf = optimize_model!(pm, optimizer = solver)

        if solution_opf["termination_status"] != MOI.LOCALLY_SOLVED
            println()
            printstyled("WARNING: the OPF did not converge!"; color = :red)
            println()
            printstyled("termination_status: ",solution_opf_0["termination_status"]; color = :red)
            global nonconverging_solutions = true
        else
            # println()
            # printstyled("Solver status OK: Ipopt converged for this OPF!"; color = :green)
            # println()
        end

        # Calculating VUF:
        vuf_allbuses = []
        for bus_i = 1:length(math["bus"])
            vm_var_i = pm.var[:it][:pmd][:nw][0][:vm][bus_i]
            va_var_i = pm.var[:it][:pmd][:nw][0][:va][bus_i]
            vuf_allbuses = vcat(vuf_allbuses,VUF_calculation(vm_var_i, va_var_i))
        end
        vuf_max = maximum(vuf_allbuses)
        global max_VUF_per_scenario = vcat(max_VUF_per_scenario, vuf_max)
    end

end
end # time

if nonconverging_solutions == false
        println()
        printstyled("Solver status OK: Ipopt converged for all simulations!"; color = :green)
        println()
else
        println()
        printstyled("Solver status WARNING: Ipopt did not converge for some simulations!"; color = :red)
        println()
end

sorted_max_VUF_per_scenario = sort(max_VUF_per_scenario, dims=1, rev=true)



# # Saving the resulting VUFs for scenarios:
results_dir = joinpath(@__DIR__, "..", "results")
if flexibility_mode in ["maximum_PQ_in_one_phase", "minimum_PQ_in_one_phase", "random_min_PQ_in_one_phase"]
    save_name = "VUFs_for_scenarios_mode_"*string(flexibility_mode)*string(phase_to_use)*".jld2"
else
    save_name = "VUFs_for_scenarios_mode_"*string(flexibility_mode)*".jld2"
end
save_path = joinpath(results_dir, save_name)
@save save_path VUF_per_scenario=sorted_max_VUF_per_scenario



# # Plot the distribution of loads across phases and scenarios:

using Plots, Plots.PlotMeasures
fz = 16

plt_load_distribution = plot(
    title = "Distribution of loads across phases and scenarios",
    xlabel = "Total active load of the network, kW",
    ylabel = "Total reactive load of the network, kVAr",

    size = (1000,1000), # width and height of the whole plot (in px)
    # aspect_ratio = :equal,

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

for phase_i = 1:3
    scatter!(plt_load_distribution,
            total_load_kW_perphase[:,phase_i],
            total_load_kVAr_perphase[:,phase_i], 
            label = ["Phase A" "Phase B" "Phase C"][phase_i],
            markerstrokewidth = 0,
            markersize = 8,
            alpha = 0.5
    )
end

display(plt_load_distribution)

plt_name = "distribution_of_loads_across_phases_and_scenarios"
savefig("../results/"*plt_name*".png")
# savefig("../results/"*plt_name*".pdf")
# savefig("../results/"*plt_name*".svg")



# Plot the distribution of maximum VUFs across scenarios:

plt_max_vuf = plot(
    title = "VUF duration curve",
    xlabel = "Number of scenarios",
    ylabel = "Maximum VUF, %",

    size = (1000,1000), # width and height of the whole plot (in px)
    # aspect_ratio = :equal,

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
plot!(plt_max_vuf,
    collect(1:length(max_VUF_per_scenario)),
    sorted_max_VUF_per_scenario*100, 
    label = "VUF",
    color = :black,
    w = 4
)

display(plt_max_vuf)

plt_name = "maximum_VUFs_sorted"
savefig("../results/"*plt_name*".png")
# savefig("../results/"*plt_name*".pdf")
# savefig("../results/"*plt_name*".svg")

