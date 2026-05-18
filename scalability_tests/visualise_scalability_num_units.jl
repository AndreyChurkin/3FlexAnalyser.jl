"""
Visualise scalability test results: OPF solve time vs number of flexible units.

Reads the CSV with all OPF times produced by scalability_test_num_units.jl and creates
the violin plot without re-running any OPF simulations.

Input:  results/scalability_tests/scalability_num_units_all_times.csv (or a backup copy)
Output: results/scalability_tests/scalability_num_units.png / .pdf

"""

using CSV, DataFrames
using Statistics
using Plots, Plots.PlotMeasures
using StatsPlots

cd(dirname(@__FILE__))

# ============================================================
# SETTINGS
# ============================================================

# Path to the all-times CSV (change to a backup filename if needed):

# csv_path = "../results/scalability_tests/scalability_num_units_all_times.csv"
csv_path = "../results/scalability_tests/scalability_num_units_all_times_2026_05_17.csv"


# ============================================================


# --- Load data ---
df = CSV.read(csv_path, DataFrame)
sort!(df, :n_units)

all_n_units = sort(unique(df.n_units))

# Reconstruct time vectors per n_units (same structure as 'results_all_times' in the test script):
results_n_units   = all_n_units
results_all_times = [df[df.n_units .== n, :time_s] for n in all_n_units]

med_times = [isempty(t) ? NaN : median(t) for t in results_all_times]

println("Loaded $(nrow(df)) OPF times for $(length(all_n_units)) unit combination(s): $all_n_units")


# --- Violin plot ---
font_size_summary = 22

violin_labels = vcat([fill(results_n_units[i], length(results_all_times[i]))
                      for i in 1:length(results_n_units)]...)
violin_times  = vcat(results_all_times...)

plt_timing = violin(violin_labels, violin_times,
    fillcolor      = :lightgrey,
    linecolor      = :grey,
    lw             = 2,
    label          = "Distribution",
    xlabel         = "Number of flexible units",
    ylabel         = "OPF solve time, s",
    size           = (1200, 600),
    framestyle     = :box,
    margin         = 10mm,
    left_margin    = 15mm,
    xticks         = all_n_units,
    legend         = :topleft,
    fontfamily     = "Courier",
    xtickfontsize  = font_size_summary,
    ytickfontsize  = font_size_summary,
    xguidefontsize = font_size_summary,
    yguidefontsize = font_size_summary,
    legendfontsize = font_size_summary - 4,

    ylim           = (0, 300)
)

plot!(plt_timing, results_n_units, med_times,
      color     = :black,
      lw        = 2,
      linestyle = :dash,
      label     = false)

scatter!(plt_timing, results_n_units, med_times,
         color       = :black,
         markersize  = 8,
         markershape = :circle,
         label       = "Median")

for (i, (n, m)) in enumerate(zip(results_n_units, med_times))
    annotate!(plt_timing, n, m + 14, text("$(round(m, digits=1))", font(font_size_summary-8, "Courier", :bold), :center))
end


savefig(plt_timing, "../results/scalability_tests/scalability_num_units.png")
savefig(plt_timing, "../results/scalability_tests/scalability_num_units.pdf")
display(plt_timing)
println("Plot saved to results/scalability_tests/scalability_num_units.png / .pdf")
