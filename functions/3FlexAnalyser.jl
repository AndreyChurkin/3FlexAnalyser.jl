"""
This is the main code for calculating P-Q flexibility areas, plotting them, and imposing voltage unbalance and phase coordination constraints.
The code reads data about the network and flexibility from /cases in the OpenDSS format (.dss).
Yet, many parameters of the simulations have to be specified in this code to produce correct simulations.
For example, it is necessary to select the case study (.dss file), limits of flexible units, phase (a, b, or c), 
location of flexibility aggregation, number of simulations (points used to estimate P-Q flexibility areas,
voltage unbalance constraints, phase coordination constraints, some plotting parameters, etc.

This code is not perfect and may be further improved. There is still ongoing work and research as we prepare a manuscript for an IEEE Transactions journal.
Nonetheless, as of August 2024, the code works reliably and reasonably fast for the two included case studies.
A more detailed description can be found on GitHub and in the manuscript [LINK!!!!!].

"""

using PowerModelsDistribution
using JuMP, Ipopt

cd(dirname(@__FILE__))

include("../functions/calculate_VUF_a_posteriori.jl")
include("../functions/build_VUF_constraint.jl")
include("../functions/build_VUF_constraint_new_all_bus.jl")
include("../functions/build_phase_coordination_constraints.jl")


# # Select a case to analyse:

# eng = parse_file("../cases/5_bus_case_illustrative/LV_balanced_flex_balanced.dss")
# eng = parse_file("../cases/5_bus_case_illustrative/LV_unbalanced_flex_balanced.dss")
eng = parse_file("../cases/5_bus_case_illustrative/LV_unbalanced_flex_unbalanced.dss")

# eng = parse_file("../cases/221_bus_real_UK_case/Master_221_bus_UK.dss")


# # Activate for adjusting transformer voltage (can be used for the UK case):
# if haskey(eng, "transformer")
#     transformer_keys = collect(keys(eng["transformer"]))
#     for trans = 1:length(eng["transformer"])
#         eng["transformer"][transformer_keys[trans]]["tm_set"][2] = [1.05, 1.05, 1.05]
#     end
# end


# # Some simulation settings used in PowerModelsDistribution.jl
eng["settings"]["sbase_default"] = 1 # if = 1, pm.model will be in kW
# eng["settings"]["voltage_scale_factor"] = 1
eng["settings"]["power_scale_factor"] = 1000


# # Set P-Q limits for flexible generators (in kW and kVAr):

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


# # Select phase to optimise (P-Q flexibility area will be maximised for this phase):
# phase_i = 1
phase_i = 2
# phase_i = 3


# # Define an allowed (feasible) voltage range:
v_ub = 1.10
v_lb = 0.94



# # Select an objective (location/element) for flexibility aggregation:

# aggregation_objective = "source" # <-- select to maximise flexibility for the source generator
aggregation_objective = "line_flow" # <-- select to maximise flexibility for a specific power flow (line)
"""
!! if selecting "line_flow", define the line's number (used in math model) !!
!! If not sure about the exact line number (index), check "math" and "solution_opf_0" !!
"""
# aggregation_line_number = 3 # for testing the 5-bus system (feeder line 1-2)
# aggregation_line_number = 6 # for testing the 5-bus system (line from source bus)
aggregation_line_number = 123 # for 221-bus UK case (line from source bus)


# # Imosing voltage unbalance constraints:
# global impose_vuf_constraints = false # <-- no additional voltage unbalance constraints (P-Q flexibility areas should not be reduced)
global impose_vuf_constraints = true # <-- impose VUF limits (P-Q flexibility areas may be reduced)

# # Set the VUF limit for simulations (if imposing voltage unbalance constraints):
# global vuf_threshold = 0.02 # note that 2% is 0.02
# global vuf_threshold = 0.015
# global vuf_threshold = 0.013
# global vuf_threshold = 0.01
# global vuf_threshold = 0.009
# global vuf_threshold = 0.008
# global vuf_threshold = 0.007
# global vuf_threshold = 0.006
global vuf_threshold = 0.005
# global vuf_threshold = 0.004
# global vuf_threshold = 0.003
# global vuf_threshold = 0.002
# global vuf_threshold = 0.001

# # Specify for which buses VUF constraints should be imposed:
# global all_buses_vuf_constrained = false # <-- if false, VUF constraint can be imposed only for a single bus
global all_buses_vuf_constrained = true # <-- if true, VUF constraints can be imposed for every bus 
"""
Note: Activating "all_buses_vuf_constrained = true" is necessary to impose constraints for multiple specific buses, e.g., using "vuf_constrained_buses" and "exclude_buses_from_vuf_constraints"
"""

# # Set buses to exclude from VUF constraints (e.g., because of transformers' voltage levels)
global exclude_buses_from_vuf_constraints = [] # <-- [] means no exclusions
# global exclude_buses_from_vuf_constraints = [48, 156, 223, 225] # UK case
# global exclude_buses_from_vuf_constraints = [46, 48, 127] # UK case
# global exclude_buses_from_vuf_constraints = collect(1:226) # (test)
# global exclude_buses_from_vuf_constraints = collect(11:226) # (test)
# global exclude_buses_from_vuf_constraints = vcat(collect(1:11),collect(40:226))

# # A fixed set of buses to impose VUF constraints (only these buses will be considered for VUF limits):
global vuf_constrained_buses = [] # <-- no specific buses defined

# # (Used in the 221-bus UK case: buses in the most unbalanced part of the network)
# global vuf_constrained_buses = [
# "bus_36049497_01"
# "bus_36067332_01"
# "bus_36067558_01"
# "bus_36049503_01"
# "bus_36049305"
# "bus_36041228_01"
# "bus_36049000_01"
# ]


# # Introduce phase coordination constraints:
# global impose_phase_coordination_constraints = false # <-- If false, no additional constraints are imposed (P-Q flexibility areas will not be reduced)
global impose_phase_coordination_constraints = true # <-- If true, phase coordination constraints are imposed (P-Q flexibility areas will be reduced significantly)


# # Now, the mathematical model is formulated below:

solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>1)
# solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>1,"tol"=>1e-9)

math = transform_data_model(eng)

# # Locating the source generator (i.e., the flexibility aggregation bus):
multiple_source_check = 0
for gen_i = 1:length(math["gen"])
    if math["gen"][string(gen_i)]["name"] == "_virtual_gen.voltage_source.source"
        global source_gen_i = gen_i
        global source_bus_i = math["gen"][string(gen_i)]["gen_bus"]
        global multiple_source_check += 1
    end
end
if multiple_source_check == 0
    println()
    printstyled("WARNING: No source generator found in the system!"; color = :red)
    println()
elseif multiple_source_check >= 2
    println()
    printstyled("WARNING: Multiple source generators found in the system! Please define a specific bus for flexibility aggregation."; color = :red)
end

# # Activate to make the source bus voltage not fixed to 1.0 pu:
# math["bus"][string(source_bus_i)]["vmin"] = [0.0, 0.0, 0.0]
# math["bus"][string(source_bus_i)]["vmax"] = [Inf, Inf, Inf]
# # delete!(math["bus"][string(source_bus_i)], "vm")
# # delete!(math["bus"][string(source_bus_i)], "va")

# # Define the VUF regulation bus:
global vuf_regulation_bus = 0
for branch_i = 1:length(math["branch"]) # find the bus next to the source bus
    if math["branch"][string(branch_i)]["name"] == "_virtual_branch.voltage_source.source"
        global vuf_regulation_bus = math["branch"][string(branch_i)]["t_bus"]
    end
end
# vuf_regulation_bus = 4 # (test for the 5-bus system)

if length(vuf_constrained_buses) >= 1 # if there is a set of buses for vuf constraints
    global exclude_buses_from_vuf_constraints = collect(1:length(math["bus"])) # exclude all buses first
    global vuf_constrained_bus_numbers = []
    for include_bus = 1:length(vuf_constrained_buses)
        global vuf_constrained_bus_numbers = vcat(vuf_constrained_bus_numbers, math["bus_lookup"][vuf_constrained_buses[include_bus]])
    end
    global indices_to_delete = sort(vuf_constrained_bus_numbers, rev=true)
    for idx in indices_to_delete
        deleteat!(exclude_buses_from_vuf_constraints, idx)
    end
end

pm = instantiate_mc_model(math, ACPUPowerModel, build_mc_opf)
# pm = instantiate_model(math)

pm.data["per_unit"] = false

for i = 1:length(math["bus"]) # define voltage limits
    @constraint(pm.model, v_lb <= pm.var[:it][:pmd][:nw][0][:vm][i][1] <= v_ub)
    @constraint(pm.model, v_lb <= pm.var[:it][:pmd][:nw][0][:vm][i][2] <= v_ub)
    @constraint(pm.model, v_lb <= pm.var[:it][:pmd][:nw][0][:vm][i][3] <= v_ub)
end
for gen_i = 1:length(math["gen"]) # set equal phase outputs for 3-phase generators
    if math["gen"][string(gen_i)]["name"] != "_virtual_gen.voltage_source.source" 
        if length(math["gen"][string(gen_i)]["connections"]) == 3
            printstyled("Note: 3-phase balanced generators are connected to the network"; color = :blue)
            @constraint(pm.model, pm.var[:it][:pmd][:nw][0][:pg][gen_i][1] == pm.var[:it][:pmd][:nw][0][:pg][gen_i][2])
            @constraint(pm.model, pm.var[:it][:pmd][:nw][0][:pg][gen_i][2] == pm.var[:it][:pmd][:nw][0][:pg][gen_i][3])

            @constraint(pm.model, pm.var[:it][:pmd][:nw][0][:qg][gen_i][1] == pm.var[:it][:pmd][:nw][0][:qg][gen_i][2])
            @constraint(pm.model, pm.var[:it][:pmd][:nw][0][:qg][gen_i][2] == pm.var[:it][:pmd][:nw][0][:qg][gen_i][3])
        end
    end
end
if impose_vuf_constraints == true
    if all_buses_vuf_constrained == false
        vm_var_i = pm.var[:it][:pmd][:nw][0][:vm][vuf_regulation_bus]
        va_var_i = pm.var[:it][:pmd][:nw][0][:va][vuf_regulation_bus]
        # vuf_constraint,vuf_constraint1,vuf_constraint2,pm_i = build_vuf_constraint(pm,vm_var_i,va_var_i,vuf_threshold)
        build_vuf_constraint_new(pm, vm_var_i, va_var_i, vuf_threshold)
    else
        N_bus = length(math["bus"])
        vm_var_allbus = pm.var[:it][:pmd][:nw][0][:vm]
        va_var_allbus = pm.var[:it][:pmd][:nw][0][:va]
        build_vuf_constraint_allbus(pm, vm_var_allbus, va_var_allbus, vuf_threshold, N_bus)
    end
end
if impose_phase_coordination_constraints == true
    build_phase_coordination_constraints(pm, phase_i)
end



# # Tests of the opf solutions (not required):
# solution_opf_1 = solve_mc_opf(eng, ACPUPowerModel, Ipopt.Optimizer)

# solution_opf_2 = optimize_model!(pm, optimizer = solver)

# objective_pg_ref_1 = pm.var[:it][:pmd][:nw][0][:pg][source_gen_i][1]
# objective_pg_ref_2 = pm.var[:it][:pmd][:nw][0][:pg][source_gen_i][2]
# objective_pg_ref_3 = pm.var[:it][:pmd][:nw][0][:pg][source_gen_i][3]

# @objective(pm.model, Min, (objective_pg_ref_1 + objective_pg_ref_2 + objective_pg_ref_3))

# solution_opf_3 = optimize_model!(pm, optimizer = solver)




# # Solving the model with no flexibility (initial operating point):

println()
println("Solving the model with no flexibility (initial operating point) ...")
eng_noflex = deepcopy(eng)
delete!(eng_noflex, "generator")
math_noflex = transform_data_model(eng_noflex)

# # Make source bus voltage not fixed to 1.0 pu
# math_noflex["bus"][string(source_bus_i)]["vmin"] = [0.0, 0.0, 0.0]
# math_noflex["bus"][string(source_bus_i)]["vmax"] = [Inf, Inf, Inf]
# # delete!(math_noflex["bus"][string(source_bus_i)], "vm")
# # delete!(math_noflex["bus"][string(source_bus_i)], "va")

pm_noflex = instantiate_mc_model(eng_noflex, ACPUPowerModel, build_mc_opf)

for i = 1:length(math_noflex["bus"]) # define voltage limits
    @constraint(pm_noflex.model, v_lb <= pm_noflex.var[:it][:pmd][:nw][0][:vm][i][1] <= v_ub)
    @constraint(pm_noflex.model, v_lb <= pm_noflex.var[:it][:pmd][:nw][0][:vm][i][2] <= v_ub)
    @constraint(pm_noflex.model, v_lb <= pm_noflex.var[:it][:pmd][:nw][0][:vm][i][3] <= v_ub)
end
for gen_i = 1:length(math_noflex["gen"]) # set equal phase outputs for 3-phase generators
    if math_noflex["gen"][string(gen_i)]["name"] != "_virtual_gen.voltage_source.source" 
        if length(math["gen"][string(gen_i)]["connections"]) == 3
            printstyled("Note: 3-phase balanced generators are connected to the network"; color = :blue)
            @constraint(pm_noflex.model, pm_noflex.var[:it][:pmd][:nw][0][:pg][gen_i][1] == pm_noflex.var[:it][:pmd][:nw][0][:pg][gen_i][2])
            @constraint(pm_noflex.model, pm_noflex.var[:it][:pmd][:nw][0][:pg][gen_i][2] == pm_noflex.var[:it][:pmd][:nw][0][:pg][gen_i][3])

            @constraint(pm_noflex.model, pm_noflex.var[:it][:pmd][:nw][0][:qg][gen_i][1] == pm_noflex.var[:it][:pmd][:nw][0][:qg][gen_i][2])
            @constraint(pm_noflex.model, pm_noflex.var[:it][:pmd][:nw][0][:qg][gen_i][2] == pm_noflex.var[:it][:pmd][:nw][0][:qg][gen_i][3])
        end
    end
end
# # (Note: VUF constraint can be too restrictive for the initital OPF with no flexibility)
# if impose_vuf_constraints == true
#     vm_var_noflex = pm_noflex.var[:it][:pmd][:nw][0][:vm][vuf_regulation_bus]
#     va_var_noflex = pm_noflex.var[:it][:pmd][:nw][0][:va][vuf_regulation_bus]
#     # vuf_constraint,vuf_constraint1,vuf_constraint2,pm_i = build_vuf_constraint(pm_i,vm_var_i,va_var_i,vuf_threshold)
#     build_vuf_constraint_new(pm_noflex, vm_var_noflex, va_var_noflex, vuf_threshold)
# end

solution_opf_0 = optimize_model!(pm_noflex, optimizer = solver)
# # Writing down the results for the initial operating point (P and Q for each phase):
if aggregation_objective == "source"
    global flex_area_results_0 = [
        solution_opf_0["solution"]["gen"]["1"]["pg"][1]  solution_opf_0["solution"]["gen"]["1"]["qg"][1]
        solution_opf_0["solution"]["gen"]["1"]["pg"][2]  solution_opf_0["solution"]["gen"]["1"]["qg"][2]
        solution_opf_0["solution"]["gen"]["1"]["pg"][3]  solution_opf_0["solution"]["gen"]["1"]["qg"][3]
    ]
elseif aggregation_objective == "line_flow"
    global flex_area_results_0 = [
        solution_opf_0["solution"]["branch"][string(aggregation_line_number)]["pf"][1]  solution_opf_0["solution"]["branch"][string(aggregation_line_number)]["qf"][1]
        solution_opf_0["solution"]["branch"][string(aggregation_line_number)]["pf"][2]  solution_opf_0["solution"]["branch"][string(aggregation_line_number)]["qf"][2]
        solution_opf_0["solution"]["branch"][string(aggregation_line_number)]["pf"][3]  solution_opf_0["solution"]["branch"][string(aggregation_line_number)]["qf"][3]
    ]
end

# # Calculating VUF:
vm_var_0 = pm_noflex.var[:it][:pmd][:nw][0][:vm][vuf_regulation_bus]
va_var_0 = pm_noflex.var[:it][:pmd][:nw][0][:va][vuf_regulation_bus]
global vuf_0 = VUF_calculation(vm_var_0, va_var_0) # !!! Note that only VUF for the regulation bus is computed here !!!
global vuf_0_allbuses = []
for bus_i = 1:length(math_noflex["bus"])
    vm_var_0_i = pm_noflex.var[:it][:pmd][:nw][0][:vm][bus_i]
    va_var_0_i = pm_noflex.var[:it][:pmd][:nw][0][:va][bus_i]
    global vuf_0_allbuses = vcat(vuf_0_allbuses,VUF_calculation(vm_var_0_i, va_var_0_i))
end

# # Use this to analyse the OPF solution (dictionary)
if solution_opf_0["termination_status"] != MOI.LOCALLY_SOLVED
        println()
        printstyled("WARNING: Initial OPF with no flexible units did not converge!"; color = :red)
        println()
        printstyled("termination_status: ",solution_opf_0["termination_status"]; color = :red)
end
# solution_opf_0_bus_keys = collect(keys(solution_opf_0["solution"]["bus"]))
# for key = 1:length(solution_opf_0["solution"]["bus"])
#     # println("bus: ",solution_opf_0_bus_keys[key],"  ", solution_opf_0["solution"]["bus"][solution_opf_0_bus_keys[Int(key)]])
#     println("bus: ",solution_opf_0_bus_keys[key]," minimum vm: ", minimum(solution_opf_0["solution"]["bus"][solution_opf_0_bus_keys[key]]["vm"]))
# end


## Building the flexibility areas:

# # Number of intervals:
# K = 2
# K = 3
# K = 5
K = 10
# K = 15
# K = 20 # used for the figures in the paper
# K = 30

@time begin

global flex_area_results = zeros(2, 2) # combinations of P and Q for ref. node
global uncertain_results = hcat([],[]) # questionable results to check (e.g., convergence issues)
global term_status_records = [] # to track solver's convergence

if aggregation_objective == "line_flow" # set the power flow variables' index
    global aggregation_line_index = (aggregation_line_number, 
                                    math["branch"][string(aggregation_line_number)]["f_bus"],
                                    math["branch"][string(aggregation_line_number)]["t_bus"])
end

# # Loop #1: iterating between Qmax and Qmin operational limits
println()
println("Solving the OPF model to find Qmin ...")
# println()
if aggregation_objective == "source"
    @objective(pm.model, Min, pm.var[:it][:pmd][:nw][0][:qg][source_gen_i][phase_i])
elseif aggregation_objective == "line_flow"
    @objective(pm.model, Min, pm.var[:it][:pmd][:nw][0][:q][aggregation_line_index][phase_i])
end    

solution_opf_i = optimize_model!(pm, optimizer = solver)

if aggregation_objective == "source"
    flex_area_results[1,1] = solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i]
    flex_area_results[1,2] = solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]
    global Qmin = solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]
elseif aggregation_objective == "line_flow"
    flex_area_results[1,1] = solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i]
    flex_area_results[1,2] = solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]
    global Qmin = solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]
end
term_status_records = vcat(term_status_records, solution_opf_i["termination_status"])

println()
println("Solving the OPF model to find Qmax ...")
# println()
if aggregation_objective == "source"
    @objective(pm.model, Max, pm.var[:it][:pmd][:nw][0][:qg][source_gen_i][phase_i])
elseif aggregation_objective == "line_flow"
    @objective(pm.model, Max, pm.var[:it][:pmd][:nw][0][:q][aggregation_line_index][phase_i])
end    

solution_opf_i = optimize_model!(pm, optimizer = solver)

if aggregation_objective == "source"
    flex_area_results[2,1] = solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i]
    flex_area_results[2,2] = solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]
    global Qmax = solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]
elseif aggregation_objective == "line_flow"
    flex_area_results[2,1] = solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i]
    flex_area_results[2,2] = solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]
    global Qmax = solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]
end

term_status_records = vcat(term_status_records, solution_opf_i["termination_status"])

global it = 0 # count iterations
for q_interval = range(Qmin, stop = Qmax, length = K)
    global it += 1
    # println()
    println("Solving the OPF model for q_interval # ", it)
    # println()

    local pm_i = instantiate_mc_model(math, ACPUPowerModel, build_mc_opf)
    pm_i.data["per_unit"] = false

    if aggregation_objective == "source" # interval (epsilon) constraints for Q intervals
        @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:qg][source_gen_i][phase_i] <= q_interval + 10^-6)
        @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:qg][source_gen_i][phase_i] >= q_interval - 10^-6)
    elseif aggregation_objective == "line_flow"
        @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:q][aggregation_line_index][phase_i] <= q_interval + 10^-6)
        @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:q][aggregation_line_index][phase_i] >= q_interval - 10^-6)
    end

    for i = 1:length(math["bus"]) # define voltage limits
        @constraint(pm_i.model, v_lb <= pm_i.var[:it][:pmd][:nw][0][:vm][i][1] <= v_ub)
        @constraint(pm_i.model, v_lb <= pm_i.var[:it][:pmd][:nw][0][:vm][i][2] <= v_ub)
        @constraint(pm_i.model, v_lb <= pm_i.var[:it][:pmd][:nw][0][:vm][i][3] <= v_ub)
    end
    for gen_i = 1:length(math["gen"]) # set equal phase outputs for 3-phase generators
        if math["gen"][string(gen_i)]["name"] != "_virtual_gen.voltage_source.source" 
            if length(math["gen"][string(gen_i)]["connections"]) == 3
                printstyled("Note: 3-phase balanced generators are connected to the network"; color = :blue)
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:pg][gen_i][1] == pm_i.var[:it][:pmd][:nw][0][:pg][gen_i][2])
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:pg][gen_i][2] == pm_i.var[:it][:pmd][:nw][0][:pg][gen_i][3])
    
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:qg][gen_i][1] == pm_i.var[:it][:pmd][:nw][0][:qg][gen_i][2])
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:qg][gen_i][2] == pm_i.var[:it][:pmd][:nw][0][:qg][gen_i][3])
            end
        end
    end
    if impose_vuf_constraints == true
        if all_buses_vuf_constrained == false
            local vm_var_i = pm_i.var[:it][:pmd][:nw][0][:vm][vuf_regulation_bus]
            local va_var_i = pm_i.var[:it][:pmd][:nw][0][:va][vuf_regulation_bus]
            # vuf_constraint,vuf_constraint1,vuf_constraint2,pm_i = build_vuf_constraint(pm_i,vm_var_i,va_var_i,vuf_threshold)
            build_vuf_constraint_new(pm_i, vm_var_i, va_var_i, vuf_threshold)
        else
            local N_bus = length(math["bus"])
            local vm_var_allbus = pm_i.var[:it][:pmd][:nw][0][:vm]
            local va_var_allbus = pm_i.var[:it][:pmd][:nw][0][:va]
            build_vuf_constraint_allbus(pm_i, vm_var_allbus, va_var_allbus, vuf_threshold, N_bus)
        end
    end
    if impose_phase_coordination_constraints == true
        build_phase_coordination_constraints(pm_i, phase_i)
    end

    if aggregation_objective == "source"
        @objective(pm_i.model, Min, pm_i.var[:it][:pmd][:nw][0][:pg][source_gen_i][phase_i])
    elseif aggregation_objective == "line_flow"
        @objective(pm_i.model, Min, pm_i.var[:it][:pmd][:nw][0][:p][aggregation_line_index][phase_i])
    end  

    local solution_opf_i = optimize_model!(pm_i, optimizer = solver)

    global term_status_records = vcat(term_status_records, solution_opf_i["termination_status"])

    if aggregation_objective == "source"
        if solution_opf_i["termination_status"] == MOI.LOCALLY_SOLVED
            global flex_area_results = vcat(flex_area_results, [solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i] solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]])
        end
        if solution_opf_i["termination_status"] != MOI.LOCALLY_SOLVED
            global uncertain_results = vcat(uncertain_results, [solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i] solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]])
        end
    elseif aggregation_objective == "line_flow"
        if solution_opf_i["termination_status"] == MOI.LOCALLY_SOLVED
            global flex_area_results = vcat(flex_area_results, [solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i] solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]])
        end
        if solution_opf_i["termination_status"] != MOI.LOCALLY_SOLVED
            global uncertain_results = vcat(uncertain_results, [solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i] solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]])
        end
    end

    if aggregation_objective == "source"
        @objective(pm_i.model, Max, pm_i.var[:it][:pmd][:nw][0][:pg][source_gen_i][phase_i])
    elseif aggregation_objective == "line_flow"
        @objective(pm_i.model, Max, pm_i.var[:it][:pmd][:nw][0][:p][aggregation_line_index][phase_i])
    end  

    local solution_opf_i = optimize_model!(pm_i, optimizer = solver)

    global term_status_records = vcat(term_status_records, solution_opf_i["termination_status"])
    
    if aggregation_objective == "source"
        if solution_opf_i["termination_status"] == MOI.LOCALLY_SOLVED
            global flex_area_results = vcat(flex_area_results, [solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i] solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]])
        end
        if solution_opf_i["termination_status"] != MOI.LOCALLY_SOLVED
            global uncertain_results = vcat(uncertain_results, [solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i] solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]])
        end
    elseif aggregation_objective == "line_flow"
        if solution_opf_i["termination_status"] == MOI.LOCALLY_SOLVED
            global flex_area_results = vcat(flex_area_results, [solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i] solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]])
        end
        if solution_opf_i["termination_status"] != MOI.LOCALLY_SOLVED
            global uncertain_results = vcat(uncertain_results, [solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i] solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]])
        end
    end
end

# # Loop #2: iterating between Pmax and Pmin operational limits
println()
println("Solving the OPF model to find Pmin ...")
# println()
if aggregation_objective == "source"
    @objective(pm.model, Min, pm.var[:it][:pmd][:nw][0][:pg][source_gen_i][phase_i])
elseif aggregation_objective == "line_flow"
    @objective(pm.model, Min, pm.var[:it][:pmd][:nw][0][:p][aggregation_line_index][phase_i])
end  

solution_opf_i = optimize_model!(pm, optimizer = solver)

if aggregation_objective == "source"
    flex_area_results = vcat(flex_area_results, [solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i] solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]])
    global Pmin = solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i]
elseif aggregation_objective == "line_flow"
    flex_area_results = vcat(flex_area_results, [solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i] solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]])
    global Pmin = solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i]
end

term_status_records = vcat(term_status_records, solution_opf_i["termination_status"])

println()
println("Solving the OPF model to find Pmax ...")
# println()
if aggregation_objective == "source"
    @objective(pm.model, Max, pm.var[:it][:pmd][:nw][0][:pg][source_gen_i][phase_i])
elseif aggregation_objective == "line_flow"
    @objective(pm.model, Max, pm.var[:it][:pmd][:nw][0][:p][aggregation_line_index][phase_i])
end  

solution_opf_i = optimize_model!(pm, optimizer = solver)

if aggregation_objective == "source"
    flex_area_results = vcat(flex_area_results, [solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i] solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]])
    global Pmax = solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i]
elseif aggregation_objective == "line_flow"
    flex_area_results = vcat(flex_area_results, [solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i] solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]])
    global Pmax = solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i]
end

term_status_records = vcat(term_status_records, solution_opf_i["termination_status"])

global it = 0 # count iterations
for p_interval = range(Pmin, stop = Pmax, length = K)
    global it += 1
    # println()
    println("Solving the OPF model for p_interval # ", it)
    # println()

    local pm_i = instantiate_mc_model(math, ACPUPowerModel, build_mc_opf)
    pm_i.data["per_unit"] = false

    if aggregation_objective == "source" # interval (epsilon) constraints for P intervals
        @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:pg][source_gen_i][phase_i] <= p_interval + 10^-6)
        @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:pg][source_gen_i][phase_i] >= p_interval - 10^-6)
    elseif aggregation_objective == "line_flow"
        @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:p][aggregation_line_index][phase_i] <= p_interval + 10^-6)
        @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:p][aggregation_line_index][phase_i] >= p_interval - 10^-6)
    end

    for i = 1:length(math["bus"]) # define voltage limits
        @constraint(pm_i.model, v_lb <= pm_i.var[:it][:pmd][:nw][0][:vm][i][1] <= v_ub)
        @constraint(pm_i.model, v_lb <= pm_i.var[:it][:pmd][:nw][0][:vm][i][2] <= v_ub)
        @constraint(pm_i.model, v_lb <= pm_i.var[:it][:pmd][:nw][0][:vm][i][3] <= v_ub)
    end
    for gen_i = 1:length(math["gen"]) # set equal phase outputs for 3-phase generators
        if math["gen"][string(gen_i)]["name"] != "_virtual_gen.voltage_source.source" 
            if length(math["gen"][string(gen_i)]["connections"]) == 3
                printstyled("Note: 3-phase balanced generators are connected to the network"; color = :blue)
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:pg][gen_i][1] == pm_i.var[:it][:pmd][:nw][0][:pg][gen_i][2])
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:pg][gen_i][2] == pm_i.var[:it][:pmd][:nw][0][:pg][gen_i][3])
    
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:qg][gen_i][1] == pm_i.var[:it][:pmd][:nw][0][:qg][gen_i][2])
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:qg][gen_i][2] == pm_i.var[:it][:pmd][:nw][0][:qg][gen_i][3])
            end
        end
    end
    if impose_vuf_constraints == true
        if all_buses_vuf_constrained == false
            local vm_var_i = pm_i.var[:it][:pmd][:nw][0][:vm][vuf_regulation_bus]
            local va_var_i = pm_i.var[:it][:pmd][:nw][0][:va][vuf_regulation_bus]
            # vuf_constraint,vuf_constraint1,vuf_constraint2,pm_i = build_vuf_constraint(pm_i,vm_var_i,va_var_i,vuf_threshold)
            build_vuf_constraint_new(pm_i, vm_var_i, va_var_i, vuf_threshold)
        else
            local N_bus = length(math["bus"])
            local vm_var_allbus = pm_i.var[:it][:pmd][:nw][0][:vm]
            local va_var_allbus = pm_i.var[:it][:pmd][:nw][0][:va]
            build_vuf_constraint_allbus(pm_i, vm_var_allbus, va_var_allbus, vuf_threshold, N_bus)
        end
    end
    if impose_phase_coordination_constraints == true
        build_phase_coordination_constraints(pm_i, phase_i)
    end

    if aggregation_objective == "source"
        @objective(pm_i.model, Min, pm_i.var[:it][:pmd][:nw][0][:qg][source_gen_i][phase_i])
    elseif aggregation_objective == "line_flow"
        @objective(pm_i.model, Min, pm_i.var[:it][:pmd][:nw][0][:q][aggregation_line_index][phase_i])
    end  

    local solution_opf_i = optimize_model!(pm_i, optimizer = solver)

    global term_status_records = vcat(term_status_records, solution_opf_i["termination_status"])

    if aggregation_objective == "source"
        if solution_opf_i["termination_status"] == MOI.LOCALLY_SOLVED
            global flex_area_results = vcat(flex_area_results, [solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i] solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]])
        end
        if solution_opf_i["termination_status"] != MOI.LOCALLY_SOLVED
            global uncertain_results = vcat(uncertain_results, [solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i] solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]])
        end
    elseif aggregation_objective == "line_flow"
        if solution_opf_i["termination_status"] == MOI.LOCALLY_SOLVED
            global flex_area_results = vcat(flex_area_results, [solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i] solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]])
        end
        if solution_opf_i["termination_status"] != MOI.LOCALLY_SOLVED
            global uncertain_results = vcat(uncertain_results, [solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i] solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]])
        end
    end

    if aggregation_objective == "source"
        @objective(pm_i.model, Max, pm_i.var[:it][:pmd][:nw][0][:qg][source_gen_i][phase_i])
    elseif aggregation_objective == "line_flow"
        @objective(pm_i.model, Max, pm_i.var[:it][:pmd][:nw][0][:q][aggregation_line_index][phase_i])
    end 

    local solution_opf_i = optimize_model!(pm_i, optimizer = solver)

    global term_status_records = vcat(term_status_records, solution_opf_i["termination_status"])

    if aggregation_objective == "source"
        if solution_opf_i["termination_status"] == MOI.LOCALLY_SOLVED
            global flex_area_results = vcat(flex_area_results, [solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i] solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]])
        end
        if solution_opf_i["termination_status"] != MOI.LOCALLY_SOLVED
            global uncertain_results = vcat(uncertain_results, [solution_opf_i["solution"]["gen"][string(source_gen_i)]["pg"][phase_i] solution_opf_i["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]])
        end
    elseif aggregation_objective == "line_flow"
        if solution_opf_i["termination_status"] == MOI.LOCALLY_SOLVED
            global flex_area_results = vcat(flex_area_results, [solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i] solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]])
        end
        if solution_opf_i["termination_status"] != MOI.LOCALLY_SOLVED
            global uncertain_results = vcat(uncertain_results, [solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i] solution_opf_i["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]])
        end
    end
end

# # Report the validity of the OPF with no flexible units:
if solution_opf_0["termination_status"] != MOI.LOCALLY_SOLVED
    println()
    printstyled("WARNING: Initial OPF with no flexible units did not converge!"; color = :red)
end

# # Check the validity of all solutions:
infeasible_solutions = 0
for check = 1:size(term_status_records)[1]
        if term_status_records[check] != MOI.LOCALLY_SOLVED
            global infeasible_solutions += 1
            # global uncertain_results = vcat(uncertain_results,flex_area_results[check,:]')
        end
end
if infeasible_solutions == 0
        println()
        printstyled("Solver status OK: Ipopt converged for all simulations!"; color = :green)
        println()
else
        println()
        printstyled("Solver status WARNING: Ipopt did not converge for some simulations!"; color = :red)
        println()
end

end # @time

# # # Save the flexibility results to JLD:
using JLD
# jld_name = "221bus_UK_phaseA_no_constraints.jld"
# jld_name = "221bus_UK_phaseA_nocoordination_noVUF.jld"
# jld_name = "221bus_UK_phaseA_nocoordination_VUF0.005.jld"
# jld_name = "221bus_UK_phaseB_no_constraints.jld"
# jld_name = "221bus_UK_phaseB_nocoordination_noVUF.jld"
# jld_name = "221bus_UK_phaseB_nocoordination_VUF0.005.jld"
# jld_name = "221bus_UK_phaseC_no_constraints.jld"
# jld_name = "221bus_UK_phaseC_nocoordination_noVUF.jld"
# jld_name = "221bus_UK_phaseC_nocoordination_VUF0.005.jld"
jld_name = "test_1.jld"


save(jld_name
     , "flex_area_results",flex_area_results
     , "flex_area_results_0",flex_area_results_0
     , "uncertain_results",uncertain_results
     , "term_status_records",term_status_records
)

## plotting the flexibility area:
using Plots, Plots.PlotMeasures
using LazySets, Polyhedra
using ConcaveHull

# reverse_results = false
reverse_results = true # use "true" if the power flow results are negative (to plot them as positive)

if reverse_results == true
    global plot_flex_area_results = -1*flex_area_results
    global plot_uncertain_results = -1*uncertain_results
    global plot_flex_area_results_0 = -1*flex_area_results_0
else
    global plot_flex_area_results = flex_area_results
    global plot_uncertain_results = uncertain_results
    global plot_flex_area_results_0 = flex_area_results_0
end

points = N -> [plot_flex_area_results[i,:] for i in 1:N]
v = points(size(plot_flex_area_results)[1])
# hull = convex_hull(v)

# zoom_out = 1.0 # kVA
zoom_out = 2.0 # kVA

N_flex_gen = length(eng["generator"]) # number of flexible generators

# font_size = 30
font_size = 26

if phase_i == 1
    global boundary_color = palette(:tab10)[1]
elseif phase_i == 2
    global boundary_color = palette(:tab10)[2]
elseif phase_i == 3
    global boundary_color = palette(:tab10)[3]
end

plt = plot(
            # VPolygon(hull),
            alpha=0.25,
            lw = 3,
            # linecolor =
            # color=palette(:tab10)[1],
            # color=palette(:tab10)[2],
            # color=palette(:tab10)[3],
            # color = boundary_color,
            linealpha = 0.7,
            fontfamily = "Courier",
            # size = (2000,2000),
            size = (1200,1200),
            # xlim=(minimum(plot_flex_area_results[:,1])-zoom_out, maximum(plot_flex_area_results[:,1])+zoom_out),
            # ylim=(minimum(plot_flex_area_results[:,2])-zoom_out, maximum(plot_flex_area_results[:,2])+zoom_out),

            # # xlim = (10,31), # 5-bus case
            # # ylim = (-0.5, 21), # 5-bus case

            # xlim = (8,33), # 5-bus case
            # ylim = (-2.5, 23), # 5-bus case

            # xlim = (plot_flex_area_results_0[phase_i,1] - gen_lim_Pmax*N_flex_gen, plot_flex_area_results_0[phase_i,1] - gen_lim_Pmin*N_flex_gen),
            # ylim = (plot_flex_area_results_0[phase_i,2] - gen_lim_Qmax*N_flex_gen, plot_flex_area_results_0[phase_i,2] - gen_lim_Qmin*N_flex_gen),

            xlim = (plot_flex_area_results_0[phase_i,1] - 30, plot_flex_area_results_0[phase_i,1] + 30), # 221-bus UK case
            ylim = (plot_flex_area_results_0[phase_i,2] - 30, plot_flex_area_results_0[phase_i,2] + 30),

            xlabel = "P, kW", ylabel = "Q, kVAr",
            xtickfontsize = font_size, ytickfontsize = font_size,
            xguidefontsize = font_size, yguidefontsize = font_size, legendfontsize = font_size,
            foreground_color_legend = nothing, legend = false,
            framestyle = :box, margin = 20mm, left_margin = 50mm, minorgrid = :true,
            aspect_ratio=:equal
)

neighbours = 1 # The smoothness of the concave hull can be adjusted by increasing the number of neighbours used
c_hull = concave_hull(v,neighbours) # get a concave hull - to compute its area later
plot!(plt, c_hull, color = boundary_color)

scatter!(plt, plot_flex_area_results[:,1], plot_flex_area_results[:,2],
            # markersize = 8,
            markersize = 5,
            # markerstrokewidth = 1,
            markercolor = :black
)

scatter!(plt, plot_uncertain_results[:,1], plot_uncertain_results[:,2],
            # markersize = 12,
            markersize = 8,
            # markerstrokewidth = 1,
            markercolor = :red
)

scatter!(plt, [plot_flex_area_results_0[phase_i,1]], [plot_flex_area_results_0[phase_i,2]],
            markersize = 20,
            # markersize = 25,
            markershape = :cross,
            markercolor = :black
)

display(plt)

println()
println("c_hull_area (kVA^2):")
c_hull_area = ConcaveHull.area(c_hull) # <-- to avoid conflicting functions
println(c_hull_area)

# figure_name = "5bus_test2"
# figure_name = "5bus_bal_flex_bal"
# figure_name = "5bus_unbal_flex_bal_phaseC"
# figure_name = "5bus_unbal_flex_unbal_phaseB_noVUF"
# figure_name = "5bus_unbal_flex_unbal_phaseB_VUF0.001"
# figure_name = "5bus_unbal_flex_unbal_phaseA_noVUF"
# figure_name = "5bus_unbal_flex_unbal_nocoordination_phaseA_noVUF"
# figure_name = "5bus_unbal_flex_unbal_nocoordination_phaseA_VUF0.001"
# figure_name = "5bus_unbal_flex_unbal_phaseB_noVUF"
# figure_name = "5bus_unbal_flex_unbal_nocoordination_phaseB_noVUF"
# figure_name = "5bus_unbal_flex_unbal_nocoordination_phaseB_VUF0.001"
# figure_name = "5bus_unbal_flex_unbal_phaseС_noVUF"
# figure_name = "5bus_unbal_flex_unbal_nocoordination_phaseС_noVUF"
# figure_name = "5bus_unbal_flex_unbal_nocoordination_phaseС_VUF0.001"

# figure_name = "221bus_UK_phaseA_no_constraints"
# figure_name = "221bus_UK_phaseA_nocoordination_noVUF"
# figure_name = "221bus_UK_phaseA_nocoordination_VUF0.005"
# figure_name = "221bus_UK_phaseB_no_constraints"
# figure_name = "221bus_UK_phaseB_nocoordination_noVUF"
# figure_name = "221bus_UK_phaseB_nocoordination_VUF0.005"
# figure_name = "221bus_UK_phaseC_no_constraints"
# figure_name = "221bus_UK_phaseC_nocoordination_noVUF"
# figure_name = "221bus_UK_phaseC_nocoordination_VUF0.005"
figure_name = "test_1"




savefig(figure_name*".svg")
savefig(figure_name*".png")
savefig(figure_name*".pdf")