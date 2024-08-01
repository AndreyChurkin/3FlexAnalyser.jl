"""
This is an independent script that can be used for analysing case studies and various OPF solutions.
Specifically, it analyses the topology of a system and outputs graph files that can be further used in Gephi software.
It also solves the OPF problem, calculates voltage unbalanced, plots the voltage distribution, prints the total load of the network, etc.
This can be a useful step prior to computing the P-Q flexibility limits in 3FlexAnalyser.jl
"""

using PowerModelsDistribution
using JuMP, Ipopt
using CSV
using DataFrames

cd(dirname(@__FILE__))
println(pwd())

include("../functions/calculate_VUF_a_posteriori.jl")


# # Select a case to analyse:

# eng = parse_file("../cases/5_bus_case_illustrative/LV_balanced_flex_balanced.dss")
# eng = parse_file("../cases/5_bus_case_illustrative/LV_unbalanced_flex_balanced.dss")
eng = parse_file("../cases/5_bus_case_illustrative/LV_unbalanced_flex_unbalanced.dss")

# eng = parse_file("../cases/221_bus_real_UK_case/Master_221_bus_UK.dss")


eng["settings"]["sbase_default"] = 1 # if = 1, pm.model will be in kW
# eng["settings"]["voltage_scale_factor"] = 1
eng["settings"]["power_scale_factor"] = 1000

solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>1,"tol"=>1e-9)

# # Remove flexible generators (to analyse the initial network operation without flexibility):
eng_noflex = deepcopy(eng)
delete!(eng_noflex, "generator")
math_noflex = transform_data_model(eng_noflex)
bus_lookup_keys = collect(keys(math_noflex["bus_lookup"]))


# # Use this part to find specific lines for further analysis:

# global find_bus_1 = "bus_36085587_36085583"
# # global find_bus_1 = "bus_36085579"
# # global find_bus_1 = "_virtual_gen.voltage_source.source"

# for check_bus = 1:length(bus_lookup_keys)
#     if bus_lookup_keys[check_bus] == find_bus_1
#         global find_bus_1_index = check_bus
#     end
# end
# for check_branch = 1:length(math_noflex["branch"])
#     if math_noflex["branch"][string(check_branch)]["f_bus"] == find_bus_1_index
#         println()
#         printstyled("Bus "*find_bus_1*" found in connections of line: "*string(check_branch); color = :blue)
#     elseif math_noflex["branch"][string(check_branch)]["t_bus"] == find_bus_1_index
#         println()
#         printstyled("Bus "*find_bus_1*" found in connections of line: "*string(check_branch); color = :blue)
#     end
# end


pm_noflex = instantiate_mc_model(eng_noflex, ACPUPowerModel, build_mc_opf)
solution_opf_0 = optimize_model!(pm_noflex, optimizer = solver)

if solution_opf_0["termination_status"] != MOI.LOCALLY_SOLVED
    println()
    printstyled("WARNING: Initial OPF with no flexible units did not converge!"; color = :red)
    println()
    printstyled("termination_status: ",solution_opf_0["termination_status"]; color = :red)
else
    println()
    printstyled("Solver status OK: Ipopt converged this OPF!"; color = :green)
    println()
end


# # Calculating VUF:
global vuf_0_allbuses = []
for bus_i = 1:length(math_noflex["bus"])
    vm_var_0_i = pm_noflex.var[:it][:pmd][:nw][0][:vm][bus_i]
    va_var_0_i = pm_noflex.var[:it][:pmd][:nw][0][:va][bus_i]
    global vuf_0_allbuses = vcat(vuf_0_allbuses,VUF_calculation(vm_var_0_i, va_var_0_i))
end

vuf_0_allbuses_sorted = sort(vuf_0_allbuses, rev=true)
vuf_0_allbuses_sorted_indices = Any[]
for bus = 1:length(vuf_0_allbuses_sorted)
    global vuf_0_allbuses_sorted_indices = vcat(vuf_0_allbuses_sorted_indices, findfirst(isequal(vuf_0_allbuses_sorted[bus]), vuf_0_allbuses))
end

# # Identify buses with the highest voltage unbalances:
check_N_top_unbalances = 5
# check_N_top_unbalances = 20

vuf_0_top_unbalance_buses_names = Any[]
for bus = 1:check_N_top_unbalances
    global vuf_0_top_unbalance_buses_names = vcat(vuf_0_top_unbalance_buses_names, bus_lookup_keys[vuf_0_allbuses_sorted_indices[bus]])
end


global get_nodes = ["Id" "Label" "kVA" "min_Vm_kV"]
global total_load_kW = 0
global total_load_kVAr = 0
global total_load_kVA = 0
global get_arcs = ["Source" "Target" "Weight"]
global total_load_kW_perphase = [0.0, 0.0, 0.0]
global total_load_kVAr_perphase = [0.0, 0.0, 0.0]

eng_bus_keys = collect(keys(eng["bus"]))
eng_load_keys = collect(keys(eng["load"]))
eng_line_keys = collect(keys(eng["line"]))
solution_opf_0_bus_keys = collect(keys(solution_opf_0["solution"]["bus"]))
solution_opf_0_bus_keys_sorted = sort(solution_opf_0_bus_keys, lt = (x, y) -> parse(Int, x) < parse(Int, y))


for bus = 1:length(eng["bus"])
    bus_ID = bus
    bus_label = eng_bus_keys[bus]
    bus_load = 0
    # global get_nodes = vcat(get_nodes,[bus bus_label bus_load minimum(solution_opf_0["solution"]["bus"][solution_opf_0_bus_keys[bus]]["vm"])])
    global get_nodes = vcat(get_nodes,[bus_ID bus_label bus_load 1.0])
end

# # (print nodes for analysis)
# for i = 1:length(get_nodes[:,1])
#     println(get_nodes[i,:])
# end

# # (export voltages for existing buses)
# for i = 1:length(solution_opf_0["solution"]["bus"])
#     findfirst(isequal("407"), get_nodes[:,2])
# end

for load = 1:length(eng["load"])
    load_kVA = sqrt(sum(eng["load"][eng_load_keys[load]]["pd_nom"])^2 + sum(eng["load"][eng_load_keys[load]]["qd_nom"])^2)
    connected_to_bus = eng["load"][eng_load_keys[load]]["bus"]
    bus_in_get_nodes = findfirst(isequal(connected_to_bus), get_nodes[:,2])
    # get_nodes[bus_in_get_nodes,3] = load_kVA
    get_nodes[bus_in_get_nodes,3] = get_nodes[bus_in_get_nodes,3] + load_kVA

    global total_load_kW += sum(eng["load"][eng_load_keys[load]]["pd_nom"])
    global total_load_kVAr += sum(eng["load"][eng_load_keys[load]]["qd_nom"])
    global total_load_kVA += sqrt(sum(eng["load"][eng_load_keys[load]]["pd_nom"])^2 + sum(eng["load"][eng_load_keys[load]]["qd_nom"])^2)

    if length(eng["load"][eng_load_keys[load]]["connections"]) >= 3
        for phase = 1:3
            global add_load_kW = eng["load"][eng_load_keys[load]]["pd_nom"][phase]
            global add_load_kVAr = eng["load"][eng_load_keys[load]]["qd_nom"][phase]
            global total_load_kW_perphase[phase] += add_load_kW
            global total_load_kVAr_perphase[phase] += add_load_kVAr
        end
    else
        global add_load_kW = eng["load"][eng_load_keys[load]]["pd_nom"][1]
        global add_load_kVAr = eng["load"][eng_load_keys[load]]["qd_nom"][1]
        global total_load_kW_perphase[eng["load"][eng_load_keys[load]]["connections"][1]] += add_load_kW
        global total_load_kVAr_perphase[eng["load"][eng_load_keys[load]]["connections"][1]] += add_load_kVAr
    end
end

for line = 1:length(eng["line"])
    from_bus_position = findfirst(isequal(eng["line"][eng_line_keys[line]]["f_bus"]), get_nodes[:,2])
    to_bus_postion = findfirst(isequal(eng["line"][eng_line_keys[line]]["t_bus"]), get_nodes[:,2])

    from_bus = get_nodes[from_bus_position,1]
    to_bus = get_nodes[to_bus_postion,1]
    global get_arcs = vcat(get_arcs,[from_bus to_bus 1])
end


# # Printing some information about the OPF solution:
println()
println("total_load_kW: ",total_load_kW)
println("total_load_kVAr: ",total_load_kVAr)
println()
println("total_load_kW_perphase: ",total_load_kW_perphase')
println("total_load_kVAr_perphase: ",total_load_kVAr_perphase')
println()
# # Use this part to print power flow through a specific line
# # e.g., branch [123] in the 221-bus UK system (supplying the network), or branch [4] in the 5-bus system
println("solution_opf_0[solution][branch][4][pf]: ")
println(solution_opf_0["solution"]["branch"]["4"]["pf"])
println("solution_opf_0[solution][branch][4][qf]: ")
println(solution_opf_0["solution"]["branch"]["4"]["qf"])
println()
# # Use this part to print power of a generator (or source):
println("solution_opf_0[solution][gen][1][pg]: ")
println(solution_opf_0["solution"]["gen"]["1"]["pg"])
println("solution_opf_0[solution][gen][1][qg]: ")
println(solution_opf_0["solution"]["gen"]["1"]["qg"])
println()
println("The highest VUF among all buses:: ")
println(maximum(vuf_0_allbuses))
println("The highest VUF bus index: ",findfirst(x -> x == maximum(vuf_0_allbuses), vuf_0_allbuses))
println()


# # Saving topology files in CSV - can be used in Gephi software
CSV.write("../results/Arcs.csv", DataFrame(get_arcs,:auto))
CSV.write("../results/Nodes.csv", DataFrame(get_nodes,:auto))
CSV.write("../results/vuf_0_allbuses.csv", DataFrame(Column1 = vuf_0_allbuses))


# Let's plot the distribution of voltages and voltage unbalances:

using Plots, Plots.PlotMeasures
fz = 18

plt_vuf = plot(
    title = "Distribution of voltage unbalances",
    xlabel = "bus #",
    ylabel = "VUF, %",

    size = (1000,1000), # width and height of the whole plot (in px)
    # aspect_ratio = :equal,

    xtickfontsize=fz, ytickfontsize=fz,
    fontfamily = "Courier", 
    titlefontsize = fz,
    xguidefontsize = fz,
    yguidefontsize = fz,

    # legend = false,
    
    framestyle = :box,
    margin = 10mm,

    minorgrid = :true,
)
plot!(plt_vuf,
    collect(1:length(vuf_0_allbuses)),
    vuf_0_allbuses*100, 
    label = "VUF"
)

display(plt_vuf)
savefig("../results/plt_vuf.png")
savefig("../results/plt_vuf.pdf")
savefig("../results/plt_vuf.svg")


plt_vm_per_phase = plot(
    title = "Voltage distribution (per bus and phase)",
    xlabel = "bus #",
    ylabel = "Voltage magnitude, p.u.",

    size = (1000,1000), # width and height of the whole plot (in px)
    # aspect_ratio = :equal,

    xtickfontsize=fz, ytickfontsize=fz,
    fontfamily = "Courier", 
    titlefontsize = fz,
    xguidefontsize = fz,
    yguidefontsize = fz,

    # legend = false,
    
    framestyle = :box,
    margin = 10mm,
    
    minorgrid = :true,
)

# # Use order from solution_opf_0:
# values_Vm_phaseA = map(dict -> get(dict, "vm", nothing)[1], values(solution_opf_0["solution"]["bus"]))
# values_Vm_phaseB = map(dict -> get(dict, "vm", nothing)[2], values(solution_opf_0["solution"]["bus"]))
# values_Vm_phaseC = map(dict -> get(dict, "vm", nothing)[3], values(solution_opf_0["solution"]["bus"]))

# # Use sorted numbering of buses:
solution_opf_0_sorted = Dict()
for add_bus = 1:length(solution_opf_0["solution"]["bus"])
    solution_opf_0_sorted[string(add_bus)] = solution_opf_0["solution"]["bus"][solution_opf_0_bus_keys_sorted[add_bus]]
end
values_Vm_phaseA = Float64[]
values_Vm_phaseB = Float64[]
values_Vm_phaseC = Float64[]
for sorted_bus = 1:length(solution_opf_0["solution"]["bus"])
    global values_Vm_phaseA = vcat(values_Vm_phaseA, solution_opf_0["solution"]["bus"][solution_opf_0_bus_keys_sorted[sorted_bus]]["vm"][1])
    global values_Vm_phaseB = vcat(values_Vm_phaseB, solution_opf_0["solution"]["bus"][solution_opf_0_bus_keys_sorted[sorted_bus]]["vm"][2])
    global values_Vm_phaseC = vcat(values_Vm_phaseC, solution_opf_0["solution"]["bus"][solution_opf_0_bus_keys_sorted[sorted_bus]]["vm"][3])

end

plot!(plt_vm_per_phase,
    collect(1:length(solution_opf_0["solution"]["bus"])),
    values_Vm_phaseA,
    label = "phase A"
)
plot!(plt_vm_per_phase,
    collect(1:length(solution_opf_0["solution"]["bus"])),
    values_Vm_phaseB,
    label = "phase B"
)
plot!(plt_vm_per_phase,
    collect(1:length(solution_opf_0["solution"]["bus"])),
    values_Vm_phaseC,
    label = "phase C"
)

display(plt_vm_per_phase)
savefig("../results/plt_vm_per_phase.png")
savefig("../results/plt_vm_per_phase.pdf")
savefig("../results/plt_vm_per_phase.svg")