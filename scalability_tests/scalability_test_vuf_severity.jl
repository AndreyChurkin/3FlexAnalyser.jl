"""
Scalability test: how does OPF solve time scale with VUF constraint severity?

All 12 flexible units are always active. VUF constraints are imposed at 7 buses 
(matching the 221-bus UK case study in the paper). The loop varies vuf_threshold
across the values defined in vuf_thresholds_to_test.

Only OPFs that converge (LOCALLY_SOLVED or ALMOST_LOCALLY_SOLVED) contribute to
the timing statistics. An interim P-Q area plot is saved after each threshold.

Only the 221-bus case study 'Master_221_bus_UK.dss' is considered here.

Final timing results are saved to:
  results/scalability_tests/scalability_vuf_severity_all_times.csv
  results/scalability_tests/scalability_vuf_severity_statistics.csv
  
"""

using PowerModelsDistribution
using JuMP, Ipopt
using Statistics
using CSV, DataFrames
using Plots, Plots.PlotMeasures
using StatsPlots
using LazySets, Polyhedra
using ConcaveHull

cd(dirname(@__FILE__))

include("../functions/calculate_VUF_a_posteriori.jl")
include("../functions/build_VUF_constraint.jl")
include("../functions/build_VUF_constraint_all_bus.jl")
include("../functions/build_phase_coordination_constraints.jl")


# ============================================================
# SIMULATION SETTINGS
# ============================================================

# K is the number of intervals per sweep direction (Q sweep and P sweep)
# Each interval requires 2 OPF to be solved (Pmin + Pmax, or Qmin + Qmax)
# Total timed OPF solves per threshold value = 4 * K
# The 4 extremes (Qmin, Qmax, Pmin, Pmax) are solved separately and not timed
# Example: K = 25  →  100 timed OPF solves per threshold value
K = 25

# Phase to optimise (1 = A, 2 = B, 3 = C):
phase_i = 1

# P-Q limits for all flexible generators (kW and kVAr):
gen_lim_Pmax =  5.0
gen_lim_Pmin = -5.0
gen_lim_Qmax =  5.0
gen_lim_Qmin = -5.0

# Voltage limits (pu):
v_ub = 1.10
v_lb = 0.94

# Aggregation objective: "source" or "line_flow" (see 3FlexAnalyser.jl for more details)
aggregation_objective = "line_flow"
# Line index in the math model for the 221-bus UK case (feeder head line from source bus)
# Note: branch numbering is topology-based and does not change across VUF threshold values
aggregation_line_number = 123

# VUF constraints are always imposed in this scalability test:
global impose_vuf_constraints = true
global all_buses_vuf_constrained = true
# Buses for which VUF constraints are imposed (matches the 221-bus UK case study in the paper)
global vuf_constrained_buses = [
    "bus_36049497_01",
    "bus_36067332_01",
    "bus_36067558_01",
    "bus_36049503_01",
    "bus_36049305",
    "bus_36041228_01",
    "bus_36049000_01"
]
# exclude_buses_from_vuf_constraints is computed from vuf_constrained_buses below (once, outside the loop)
global exclude_buses_from_vuf_constraints = []   # required global for build_vuf_constraint_allbus

# VUF threshold values to test (in per-unit; e.g. 0.01 = 1.0%):
# List from loosest to tightest so easier cases run first
# vuf_thresholds_to_test = [0.015, 0.014, 0.013, 0.012, 0.011, 0.010, 0.009, 0.008, 0.007, 0.006, 0.005]
vuf_thresholds_to_test = [0.01]


# Phase coordination constraints (keep false in the scalability tests):
# Note: build_phase_coordination_constraints also requires globals `math` and `source_gen_i`
global impose_phase_coordination_constraints = false

# Total number of flexible units (all units always active in this test):
N_units_total = 12

# Percentile bands can be shown around the median in the timing summary plot
# List bands outermost to innermost; band_alphas sets fill opacity for each (same order).
# Example: (0, 100) = full range (min to max), (25, 75) = interquartile range.
percentile_bands = [(0, 100), (10, 90), (25, 75)]
band_alphas      = [0.12,     0.18,     0.28    ] # can be used for percentile plots (ribbons)

# ============================================================



solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

# Boundary colour by phase (matching 3FlexAnalyser.jl):
if phase_i == 1
    boundary_color = palette(:tab10)[1]
elseif phase_i == 2
    boundary_color = palette(:tab10)[2]
else
    boundary_color = palette(:tab10)[3]
end

font_size = 26
zoom_out  = 2.0   # kVA, extra margin around the P-Q area in interim plots

# reverse_results --> Use 'true' to show negative aggregated power as positive
# The negative values appear in the 221-bus case due to the branch's indexing (opposite direction to the flow)
reverse_results = true


# --- Parse and configure the network once (topology does not change across VUF thresholds) ---
println()
println("Parsing network and configuring generators ...")

eng = parse_file("../cases/221_bus_real_UK_case/Master_221_bus_UK.dss")
eng["settings"]["sbase_default"] = 1
eng["settings"]["power_scale_factor"] = 1000

# All N_units_total generators are always active in this test:
for g in 1:N_units_total
    for ph in 1:length(eng["generator"]["g$g"]["pg_ub"])
        eng["generator"]["g$g"]["pg_ub"][ph] = gen_lim_Pmax
        eng["generator"]["g$g"]["pg_lb"][ph] = gen_lim_Pmin
        eng["generator"]["g$g"]["qg_ub"][ph] = gen_lim_Qmax
        eng["generator"]["g$g"]["qg_lb"][ph] = gen_lim_Qmin
    end
end

global math = transform_data_model(eng)   # global so build_phase_coordination_constraints can access it

# Compute exclude_buses_from_vuf_constraints: all buses except the ones in vuf_constrained_buses
# This mirrors the logic in 3FlexAnalyser.jl
if length(vuf_constrained_buses) >= 1
    global exclude_buses_from_vuf_constraints = collect(1:length(math["bus"]))
    global vuf_constrained_bus_numbers = []
    for include_bus = 1:length(vuf_constrained_buses)
        global vuf_constrained_bus_numbers = vcat(vuf_constrained_bus_numbers, math["bus_lookup"][vuf_constrained_buses[include_bus]])
    end
    global indices_to_delete = sort(vuf_constrained_bus_numbers, rev=true)
    for idx in indices_to_delete
        deleteat!(exclude_buses_from_vuf_constraints, idx)
    end
end

# Locate source generator and bus:
global source_gen_i = 0   # global so build_phase_coordination_constraints can access it
local  source_bus_i = 0
for g in 1:length(math["gen"])
    if math["gen"][string(g)]["name"] == "_virtual_gen.voltage_source.source"
        global source_gen_i = g
        source_bus_i = math["gen"][string(g)]["gen_bus"]
    end
end

# VUF regulation bus (bus adjacent to the source):
vuf_regulation_bus = 0
for br in 1:length(math["branch"])
    if math["branch"][string(br)]["name"] == "_virtual_branch.voltage_source.source"
        vuf_regulation_bus = math["branch"][string(br)]["t_bus"]
    end
end

# Aggregation line index tuple:
if aggregation_objective == "line_flow"
    aggregation_line_index = (
        aggregation_line_number,
        math["branch"][string(aggregation_line_number)]["f_bus"],
        math["branch"][string(aggregation_line_number)]["t_bus"]
    )
end

n_buses = length(math["bus"])
println("  Network parsed: $n_buses buses, $N_units_total flexible units")


# --- Compute no-flex initial operating point (once, topology does not change) ---

println()
println("Computing no-flex initial operating point ...")

eng_noflex = parse_file("../cases/221_bus_real_UK_case/Master_221_bus_UK.dss")
eng_noflex["settings"]["sbase_default"] = 1
eng_noflex["settings"]["power_scale_factor"] = 1000
delete!(eng_noflex, "generator")

pm_noflex = instantiate_mc_model(eng_noflex, ACPUPowerModel, build_mc_opf)
math_noflex = transform_data_model(eng_noflex)
for i in 1:length(math_noflex["bus"])
    @constraint(pm_noflex.model, v_lb <= pm_noflex.var[:it][:pmd][:nw][0][:vm][i][1] <= v_ub)
    @constraint(pm_noflex.model, v_lb <= pm_noflex.var[:it][:pmd][:nw][0][:vm][i][2] <= v_ub)
    @constraint(pm_noflex.model, v_lb <= pm_noflex.var[:it][:pmd][:nw][0][:vm][i][3] <= v_ub)
end
sol_noflex = optimize_model!(pm_noflex, optimizer = solver)
print("  no-flex:  status = ")
printstyled(string(sol_noflex["termination_status"]); color = sol_noflex["termination_status"] in (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED) ? :green : :red)
println()

if sol_noflex["termination_status"] ∉ (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
    printstyled("WARNING: no-flex OPF did not converge — cross marker will be at (0,0).\n"; color = :red)
    P0 = 0.0
    Q0 = 0.0
else
    if aggregation_objective == "source"
        src_key = first(k for (k,v) in sol_noflex["solution"]["gen"]
                        if haskey(v, "pg"))
        P0 = sol_noflex["solution"]["gen"][src_key]["pg"][phase_i]
        Q0 = sol_noflex["solution"]["gen"][src_key]["qg"][phase_i]
    else
        P0 = sol_noflex["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i]
        Q0 = sol_noflex["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]
    end
end
println("  Initial point → P0: $(round(P0,digits=3)) kW   Q0: $(round(Q0,digits=3)) kVAr")


# --- Storage for final timing summary ---
results_vuf_thresholds = Float64[]
results_all_times      = Vector{Float64}[]   # one entry per threshold: all converged OPF times
results_n_converged    = Int[]
results_n_total        = Int[]


for thresh in vuf_thresholds_to_test

    global vuf_threshold = thresh   # required global for build_vuf_constraint_allbus

    vuf_pct = round(vuf_threshold * 100, digits=3)

    println()
    println("="^55)
    println("Testing VUF threshold = $(vuf_pct)%  ($(length(vuf_thresholds_to_test)) thresholds to test in total) ...")
    println("="^55)

    # Helper to extract the aggregation P and Q from a solution:
    function get_PQ(s)
        if aggregation_objective == "source"
            return s["solution"]["gen"][string(source_gen_i)]["pg"][phase_i],
                   s["solution"]["gen"][string(source_gen_i)]["qg"][phase_i]
        else
            return s["solution"]["branch"][string(aggregation_line_number)]["pf"][phase_i],
                   s["solution"]["branch"][string(aggregation_line_number)]["qf"][phase_i]
        end
    end

    # --- Find extremes (Qmin, Qmax, Pmin, Pmax) on a shared base model ---
    # These 4 solves are not timed; they set the sweep ranges and seed flex_area_results

    pm = instantiate_mc_model(math, ACPUPowerModel, build_mc_opf)
    pm.data["per_unit"] = false

    for i in 1:length(math["bus"])
        @constraint(pm.model, v_lb <= pm.var[:it][:pmd][:nw][0][:vm][i][1] <= v_ub)
        @constraint(pm.model, v_lb <= pm.var[:it][:pmd][:nw][0][:vm][i][2] <= v_ub)
        @constraint(pm.model, v_lb <= pm.var[:it][:pmd][:nw][0][:vm][i][3] <= v_ub)
    end
    if impose_vuf_constraints
        if all_buses_vuf_constrained
            build_vuf_constraint_allbus(pm, pm.var[:it][:pmd][:nw][0][:vm],
                                        pm.var[:it][:pmd][:nw][0][:va],
                                        vuf_threshold, length(math["bus"]))
        else
            build_vuf_constraint(pm, pm.var[:it][:pmd][:nw][0][:vm][vuf_regulation_bus],
                                 pm.var[:it][:pmd][:nw][0][:va][vuf_regulation_bus], vuf_threshold)
        end
    end
    if impose_phase_coordination_constraints
        build_phase_coordination_constraints(pm, phase_i)
    end

    skip_combination = false

    # Qmin:
    if aggregation_objective == "source"
        @objective(pm.model, Min, pm.var[:it][:pmd][:nw][0][:qg][source_gen_i][phase_i])
    else
        @objective(pm.model, Min, pm.var[:it][:pmd][:nw][0][:q][aggregation_line_index][phase_i])
    end
    sol = optimize_model!(pm, optimizer = solver)
    print("  Qmin:     status = ")
    printstyled(string(sol["termination_status"]); color = sol["termination_status"] in (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED) ? :green : :red)
    println()
    if sol["termination_status"] ∉ (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
        printstyled("  WARNING: Qmin did not converge for vuf_threshold=$(vuf_pct)%. Skipping.\n"; color = :red)
        skip_combination = true
    end
    Qmin_P, Qmin = get_PQ(sol)

    if !skip_combination
        # Qmax:
        if aggregation_objective == "source"
            @objective(pm.model, Max, pm.var[:it][:pmd][:nw][0][:qg][source_gen_i][phase_i])
        else
            @objective(pm.model, Max, pm.var[:it][:pmd][:nw][0][:q][aggregation_line_index][phase_i])
        end
        sol = optimize_model!(pm, optimizer = solver)
        print("  Qmax:     status = ")
        printstyled(string(sol["termination_status"]); color = sol["termination_status"] in (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED) ? :green : :red)
        println()
        if sol["termination_status"] ∉ (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
            printstyled("  WARNING: Qmax did not converge for vuf_threshold=$(vuf_pct)%. Skipping.\n"; color = :red)
            skip_combination = true
        end
        Qmax_P, Qmax = get_PQ(sol)
    end

    if !skip_combination
        # Pmin:
        if aggregation_objective == "source"
            @objective(pm.model, Min, pm.var[:it][:pmd][:nw][0][:pg][source_gen_i][phase_i])
        else
            @objective(pm.model, Min, pm.var[:it][:pmd][:nw][0][:p][aggregation_line_index][phase_i])
        end
        sol = optimize_model!(pm, optimizer = solver)
        print("  Pmin:     status = ")
        printstyled(string(sol["termination_status"]); color = sol["termination_status"] in (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED) ? :green : :red)
        println()
        if sol["termination_status"] ∉ (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
            printstyled("  WARNING: Pmin did not converge for vuf_threshold=$(vuf_pct)%. Skipping.\n"; color = :red)
            skip_combination = true
        end
        Pmin, Pmin_Q = get_PQ(sol)
    end

    if !skip_combination
        # Pmax:
        if aggregation_objective == "source"
            @objective(pm.model, Max, pm.var[:it][:pmd][:nw][0][:pg][source_gen_i][phase_i])
        else
            @objective(pm.model, Max, pm.var[:it][:pmd][:nw][0][:p][aggregation_line_index][phase_i])
        end
        sol = optimize_model!(pm, optimizer = solver)
        print("  Pmax:     status = ")
        printstyled(string(sol["termination_status"]); color = sol["termination_status"] in (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED) ? :green : :red)
        println()
        if sol["termination_status"] ∉ (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
            printstyled("  WARNING: Pmax did not converge for vuf_threshold=$(vuf_pct)%. Skipping.\n"; color = :red)
            skip_combination = true
        end
        Pmax, Pmax_Q = get_PQ(sol)
    end

    if skip_combination
        continue
    end

    println("  Extremes → Pmin: $(round(Pmin,digits=3))  Pmax: $(round(Pmax,digits=3))  Qmin: $(round(Qmin,digits=3))  Qmax: $(round(Qmax,digits=3))")

    # Seed flex_area_results with the 4 extreme points:
    flex_area_results = [Qmin_P  Qmin
                         Qmax_P  Qmax
                         Pmin    Pmin_Q
                         Pmax    Pmax_Q]
    uncertain_results = Matrix{Float64}(undef, 0, 2)

    # --- Timed interval sweeps ---
    opf_times  = Float64[]
    n_converged = 0
    n_total     = 0

    # Loop 1: Q sweep — fix Q interval, solve for Pmin and Pmax
    for q_interval in range(Qmin, stop = Qmax, length = K)
        for direction in [:Min, :Max]

            local pm_i = instantiate_mc_model(math, ACPUPowerModel, build_mc_opf)
            pm_i.data["per_unit"] = false

            for i in 1:length(math["bus"])
                @constraint(pm_i.model, v_lb <= pm_i.var[:it][:pmd][:nw][0][:vm][i][1] <= v_ub)
                @constraint(pm_i.model, v_lb <= pm_i.var[:it][:pmd][:nw][0][:vm][i][2] <= v_ub)
                @constraint(pm_i.model, v_lb <= pm_i.var[:it][:pmd][:nw][0][:vm][i][3] <= v_ub)
            end
            if impose_vuf_constraints
                if all_buses_vuf_constrained
                    build_vuf_constraint_allbus(pm_i, pm_i.var[:it][:pmd][:nw][0][:vm],
                                                pm_i.var[:it][:pmd][:nw][0][:va],
                                                vuf_threshold, length(math["bus"]))
                else
                    build_vuf_constraint(pm_i, pm_i.var[:it][:pmd][:nw][0][:vm][vuf_regulation_bus],
                                         pm_i.var[:it][:pmd][:nw][0][:va][vuf_regulation_bus], vuf_threshold)
                end
            end
            if impose_phase_coordination_constraints
                build_phase_coordination_constraints(pm_i, phase_i)
            end

            if aggregation_objective == "source"
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:qg][source_gen_i][phase_i] <= q_interval + 1e-6)
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:qg][source_gen_i][phase_i] >= q_interval - 1e-6)
                if direction == :Min
                    @objective(pm_i.model, Min, pm_i.var[:it][:pmd][:nw][0][:pg][source_gen_i][phase_i])
                else
                    @objective(pm_i.model, Max, pm_i.var[:it][:pmd][:nw][0][:pg][source_gen_i][phase_i])
                end
            else
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:q][aggregation_line_index][phase_i] <= q_interval + 1e-6)
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:q][aggregation_line_index][phase_i] >= q_interval - 1e-6)
                if direction == :Min
                    @objective(pm_i.model, Min, pm_i.var[:it][:pmd][:nw][0][:p][aggregation_line_index][phase_i])
                else
                    @objective(pm_i.model, Max, pm_i.var[:it][:pmd][:nw][0][:p][aggregation_line_index][phase_i])
                end
            end

            t_start = time()
            local sol_i = optimize_model!(pm_i, optimizer = solver)
            elapsed = time() - t_start

            p_val, q_val = get_PQ(sol_i)
            n_total += 1
            status_ok = sol_i["termination_status"] in (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
            print("    Q-sweep OPF #$n_total:  status = ")
            printstyled(string(sol_i["termination_status"]); color = status_ok ? :green : :red)
            println("   time = $(round(elapsed, digits=3)) s")
            if status_ok
                push!(opf_times, elapsed)
                n_converged += 1
                flex_area_results = vcat(flex_area_results, [p_val  q_val])
            else
                uncertain_results = vcat(uncertain_results, [p_val  q_val])
            end

        end
    end

    # Loop 2: P sweep — fix P interval, solve for Qmin and Qmax
    for p_interval in range(Pmin, stop = Pmax, length = K)
        for direction in [:Min, :Max]

            local pm_i = instantiate_mc_model(math, ACPUPowerModel, build_mc_opf)
            pm_i.data["per_unit"] = false

            for i in 1:length(math["bus"])
                @constraint(pm_i.model, v_lb <= pm_i.var[:it][:pmd][:nw][0][:vm][i][1] <= v_ub)
                @constraint(pm_i.model, v_lb <= pm_i.var[:it][:pmd][:nw][0][:vm][i][2] <= v_ub)
                @constraint(pm_i.model, v_lb <= pm_i.var[:it][:pmd][:nw][0][:vm][i][3] <= v_ub)
            end
            if impose_vuf_constraints
                if all_buses_vuf_constrained
                    build_vuf_constraint_allbus(pm_i, pm_i.var[:it][:pmd][:nw][0][:vm],
                                                pm_i.var[:it][:pmd][:nw][0][:va],
                                                vuf_threshold, length(math["bus"]))
                else
                    build_vuf_constraint(pm_i, pm_i.var[:it][:pmd][:nw][0][:vm][vuf_regulation_bus],
                                         pm_i.var[:it][:pmd][:nw][0][:va][vuf_regulation_bus], vuf_threshold)
                end
            end
            if impose_phase_coordination_constraints
                build_phase_coordination_constraints(pm_i, phase_i)
            end

            if aggregation_objective == "source"
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:pg][source_gen_i][phase_i] <= p_interval + 1e-6)
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:pg][source_gen_i][phase_i] >= p_interval - 1e-6)
                if direction == :Min
                    @objective(pm_i.model, Min, pm_i.var[:it][:pmd][:nw][0][:qg][source_gen_i][phase_i])
                else
                    @objective(pm_i.model, Max, pm_i.var[:it][:pmd][:nw][0][:qg][source_gen_i][phase_i])
                end
            else
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:p][aggregation_line_index][phase_i] <= p_interval + 1e-6)
                @constraint(pm_i.model, pm_i.var[:it][:pmd][:nw][0][:p][aggregation_line_index][phase_i] >= p_interval - 1e-6)
                if direction == :Min
                    @objective(pm_i.model, Min, pm_i.var[:it][:pmd][:nw][0][:q][aggregation_line_index][phase_i])
                else
                    @objective(pm_i.model, Max, pm_i.var[:it][:pmd][:nw][0][:q][aggregation_line_index][phase_i])
                end
            end

            t_start = time()
            local sol_i = optimize_model!(pm_i, optimizer = solver)
            elapsed = time() - t_start

            p_val, q_val = get_PQ(sol_i)
            n_total += 1
            status_ok = sol_i["termination_status"] in (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
            print("    P-sweep OPF #$n_total:  status = ")
            printstyled(string(sol_i["termination_status"]); color = status_ok ? :green : :red)
            println("   time = $(round(elapsed, digits=3)) s")
            if status_ok
                push!(opf_times, elapsed)
                n_converged += 1
                flex_area_results = vcat(flex_area_results, [p_val  q_val])
            else
                uncertain_results = vcat(uncertain_results, [p_val  q_val])
            end

        end
    end

    println("  Converged: $n_converged / $n_total timed OPFs")
    if !isempty(opf_times)
        println("  Mean OPF time: $(round(mean(opf_times), digits=3)) s  |  Std: $(round(std(opf_times), digits=3)) s")
    end

    # --- Interim P-Q area plot for this vuf_threshold ---

    # Apply sign reversal if the branch index direction gives negative flows:
    if reverse_results
        plot_flex    = -1 .* flex_area_results
        plot_uncertain = size(uncertain_results, 1) > 0 ? -1 .* uncertain_results : uncertain_results
        plot_P0, plot_Q0 = -P0, -Q0
    else
        plot_flex      = flex_area_results
        plot_uncertain = uncertain_results
        plot_P0, plot_Q0 = P0, Q0
    end

    pts    = [plot_flex[i, :] for i in 1:size(plot_flex, 1)]
    c_hull = concave_hull(pts, 1)

    plt_area = plot(
        alpha       = 0.25,
        lw          = 3,
        linealpha   = 0.7,
        fontfamily  = "Courier",
        size        = (1200, 1200),

        # Static limits used in the "Exposing Barriers" paper:
        xlim        = (plot_P0 - 30, plot_P0 + 30),
        ylim        = (plot_Q0 - 30, plot_Q0 + 30),

        xlabel      = "P, kW",
        ylabel      = "Q, kVAr",
        title          = "VUF=$(vuf_pct)%  |  K=$K  |  Phase=$phase_i",
        titlefontsize  = font_size - 4,
        titlefontweight = :bold,
        xtickfontsize  = font_size, ytickfontsize = font_size,
        xguidefontsize = font_size, yguidefontsize = font_size,
        legendfontsize = font_size,
        foreground_color_legend = nothing,
        legend         = false,
        framestyle     = :box,
        margin         = 20mm,
        top_margin     = 5mm,
        left_margin    = 50mm,
        minorgrid      = true,
        aspect_ratio   = :equal
    )

    plot!(plt_area, c_hull, color = boundary_color)

    scatter!(plt_area, plot_flex[:, 1], plot_flex[:, 2],
             markersize = 5, markercolor = :black)

    if size(plot_uncertain, 1) > 0
        scatter!(plt_area, plot_uncertain[:, 1], plot_uncertain[:, 2],
                 markersize = 8, markercolor = :red)
    end

    scatter!(plt_area, [plot_P0], [plot_Q0],
             markersize = 20, markershape = :cross, markercolor = :black)

    vuf_fname = replace(string(vuf_pct), "." => "_")
    area_fname = "../results/scalability_tests/scalability_vuf_severity_area_vuf_$(vuf_fname)"
    savefig(plt_area, area_fname * ".png") # <-- use to save each interim P-Q area plot
    savefig(plt_area, area_fname * ".pdf") # <-- use to save each interim P-Q area plot
    display(plt_area)
    println("  Area plot saved: scalability_vuf_severity_area_vuf$(vuf_fname)pct.png/.pdf")

    push!(results_vuf_thresholds, vuf_threshold)
    push!(results_all_times,      copy(opf_times))
    push!(results_n_converged,    n_converged)
    push!(results_n_total,        n_total)

end


# --- Upsert helper: load existing CSV, remove rows for vuf_threshold in current run, append new rows ---
function upsert_csv(filepath, df_new, key_col)
    if isfile(filepath)
        df_existing = CSV.read(filepath, DataFrame)
        filter!(row -> !(row[key_col] in df_new[!, key_col]), df_existing)
        df_combined = vcat(df_existing, df_new)
    else
        df_combined = df_new
    end
    sort!(df_combined, key_col)
    CSV.write(filepath, df_combined)
end

# --- Save all individual OPF times ---
df_all = DataFrame(
    vuf_threshold = vcat([fill(results_vuf_thresholds[i], length(results_all_times[i]))
                          for i in 1:length(results_vuf_thresholds)]...),
    time_s        = vcat(results_all_times...)
)
path_all = "../results/scalability_tests/scalability_vuf_severity_all_times.csv"
upsert_csv(path_all, df_all, :vuf_threshold)
println()
println("All OPF times saved/updated: results/scalability_tests/scalability_vuf_severity_all_times.csv")

# --- Save statistics summary (mean, median, std, percentiles per threshold) ---
med_times  = [isempty(t) ? NaN : median(t) for t in results_all_times]
mean_times = [isempty(t) ? NaN : mean(t)   for t in results_all_times]
std_times  = [isempty(t) ? NaN : std(t)    for t in results_all_times]

df_stats = DataFrame(
    vuf_threshold = results_vuf_thresholds,
    mean_time_s   = mean_times,
    median_time_s = med_times,
    std_time_s    = std_times,
    n_converged   = results_n_converged,
    n_total_opfs  = results_n_total
)
for (lo, hi) in percentile_bands
    df_stats[!, "p$(lo)_time_s"] = [isempty(t) ? NaN : quantile(t, lo/100) for t in results_all_times]
    df_stats[!, "p$(hi)_time_s"] = [isempty(t) ? NaN : quantile(t, hi/100) for t in results_all_times]
end
path_stats = "../results/scalability_tests/scalability_vuf_severity_statistics.csv"
upsert_csv(path_stats, df_stats, :vuf_threshold)
println("Statistics saved/updated: results/scalability_tests/scalability_vuf_severity_statistics.csv")



# --- Final timing summary plot (violin) ---
font_size_summary = 22

# Use VUF threshold as percentage for x-axis labels:
results_vuf_pct = results_vuf_thresholds .* 100

# Flatten data for violin: repeat vuf_pct label for every OPF time in that combination:
violin_labels = vcat([fill(results_vuf_pct[i], length(results_all_times[i]))
                      for i in 1:length(results_vuf_thresholds)]...)
violin_times  = vcat(results_all_times...)

plt_timing = violin(violin_labels, violin_times,
    fillcolor      = :lightgrey,
    linecolor      = :grey,
    label          = "OPF solve times",
    xlabel         = "VUF limit (%)",
    ylabel         = "OPF solve time (s)",
    size           = (1200, 600),
    framestyle     = :box,
    margin         = 10mm,
    left_margin    = 15mm,
    xticks         = results_vuf_pct,
    legend         = :topleft,
    fontfamily     = "Courier",
    xtickfontsize  = font_size_summary,
    ytickfontsize  = font_size_summary,
    xguidefontsize = font_size_summary,
    yguidefontsize = font_size_summary,
    legendfontsize = font_size_summary - 4
)

# Dashed line connecting medians:
plot!(plt_timing, results_vuf_pct, med_times,
      color     = :black,
      lw        = 2,
      linestyle = :dash,
      label     = false)

# Median markers on top:
scatter!(plt_timing, results_vuf_pct, med_times,
         color       = :black,
         markersize  = 8,
         markershape = :circle,
         label       = "Median")

savefig(plt_timing, "../results/scalability_tests/scalability_vuf_severity.png")
savefig(plt_timing, "../results/scalability_tests/scalability_vuf_severity.pdf")
display(plt_timing)
println("Timing summary plot saved.")
