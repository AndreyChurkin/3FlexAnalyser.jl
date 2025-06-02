"""
This script is for generating random scenarios for loads and flexible units.
These scenarios can then be used to solve the OPF problems and analyse the resulting voltage unbalances.
For example, use run_and_analyse_scenarios.jl for further analysis.
"""

using PowerModelsDistribution
using JuMP, Ipopt
using CSV
using DataFrames
using JLD2
using Distributions


cd(dirname(@__FILE__))
println(pwd())



# # Select a case to analyse:

# eng = parse_file("../cases/5_bus_case_illustrative/LV_balanced_flex_balanced.dss")
# eng = parse_file("../cases/5_bus_case_illustrative/LV_unbalanced_flex_balanced.dss")
# eng = parse_file("../cases/5_bus_case_illustrative/LV_unbalanced_flex_unbalanced.dss")

eng = parse_file("../cases/221_bus_real_UK_case/Master_221_bus_UK.dss")



eng["settings"]["sbase_default"] = 1 # if = 1, pm.model will be in kW
# eng["settings"]["voltage_scale_factor"] = 1
eng["settings"]["power_scale_factor"] = 1000



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

for gen_i = 1:length(eng["generator"])
    for phase = 1:length(eng["generator"]["g"*string(gen_i)]["pg_ub"])
        eng["generator"]["g"*string(gen_i)]["pg_ub"][phase] = gen_lim_Pmax
        eng["generator"]["g"*string(gen_i)]["pg_lb"][phase] = gen_lim_Pmin
        eng["generator"]["g"*string(gen_i)]["qg_ub"][phase] = gen_lim_Qmax
        eng["generator"]["g"*string(gen_i)]["qg_lb"][phase] = gen_lim_Qmin
    end
end



math = transform_data_model(eng)
bus_lookup_keys = collect(keys(math["bus_lookup"]))


# Select number of scenarios to generaste:
N_scenarios = 1000

# Set the maximum load variation for scenarios (0.2 means ±20%):
load_variation = 0.5


# # Start generating scenarios with random loads and flexible power outputs:

math_all_scenarios = Dict{String, Dict}()


for sc = 1:N_scenarios 
    global loads_copy = deepcopy(math["load"])
    global generators_copy = deepcopy(math["gen"])

    for load_i = 1:length(loads_copy)
        # println("math[load][",load_i,"][pd] = ",math["load"][string(load_i)]["pd"])
        # println("math[load][",load_i,"][qd] = ",math["load"][string(load_i)]["qd"])
        # println()

        n_phases = length(loads_copy[string(load_i)]["connections"])
        load_multipliers = rand(Uniform(1-load_variation, 1+load_variation), n_phases)

        # load_multipliers = ones(n_phases) # <-- for testing purposes

        loads_copy[string(load_i)]["pd"] .=  loads_copy[string(load_i)]["pd"] .* load_multipliers
        loads_copy[string(load_i)]["qd"] .=  loads_copy[string(load_i)]["qd"] .* load_multipliers

    end

    gen_ε = 1e-4
    for gen_i = 1:length(generators_copy)[1]
        if generators_copy[string(gen_i)]["name"] != "_virtual_gen.voltage_source.source"
            gen_n_phases = length(generators_copy[string(gen_i)]["connections"])

            gen_P_min = generators_copy[string(gen_i)]["pmin"]
            gen_P_max = generators_copy[string(gen_i)]["pmax"]
            gen_P_rand = rand(Uniform(gen_P_min[1], gen_P_max[1]), gen_n_phases)
            generators_copy[string(gen_i)]["pg"] = gen_P_rand
            generators_copy[string(gen_i)]["pmin"] = gen_P_rand .- gen_ε
            generators_copy[string(gen_i)]["pmax"] = gen_P_rand .+ gen_ε

            gen_Q_min = generators_copy[string(gen_i)]["qmin"]
            gen_Q_max = generators_copy[string(gen_i)]["qmax"]
            gen_Q_rand = rand(Uniform(gen_Q_min[1], gen_Q_max[1]), gen_n_phases)
            generators_copy[string(gen_i)]["qg"] = gen_Q_rand
            generators_copy[string(gen_i)]["qmin"] = gen_Q_rand .- gen_ε
            generators_copy[string(gen_i)]["qmax"] = gen_Q_rand .+ gen_ε

        end
    end

    math_all_scenarios["scenario_$(sc)"] = Dict(
        "load" => loads_copy,
        "gen" => generators_copy
    )

end

# # Saving the generated load and flexibility scenarios:
# results_dir = joinpath(@__DIR__, "..", "results")
# save_path = joinpath(results_dir, "load_and_flexibility_1000_scenarios_221_bus.jld2")
# @save save_path math_generated_scenarios=math_all_scenarios
