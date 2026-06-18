"""
Visualise scalability test results: OPF solve time vs number of VUF-constrained locations

Reads the CSV with all OPF times produced by scalability_test_constrained_locations.jl
and creates the violin plot without re-running any OPF simulations.

Input:  results/scalability_tests/scalability_constrained_locations_all_times.csv (or a backup copy)
Output: results/scalability_tests/scalability_constrained_locations.png / .pdf

"""

using CSV, DataFrames
using Statistics

ENV["GKSFONTPATH"] = "C:\\Windows\\Fonts"
using Plots, Plots.PlotMeasures
using StatsPlots

cd(dirname(@__FILE__))

# ============================================================
# SETTINGS
# ============================================================

# Path to the all-times CSV (change to a backup filename if needed):
csv_path = "../results/scalability_tests/scalability_constrained_locations_all_times.csv"


# ============================================================


# --- Load data ---
df = CSV.read(csv_path, DataFrame)
sort!(df, :n_locations)

all_n_locations = sort(unique(df.n_locations))

# Reconstruct time vectors per n_locations (same structure as 'results_all_times' in the test script):
results_n_locations = all_n_locations
results_all_times   = [df[df.n_locations .== n, :time_s] for n in all_n_locations]

med_times = [isempty(t) ? NaN : median(t) for t in results_all_times]

println("Loaded $(nrow(df)) OPF times for $(length(all_n_locations)) location tests: $all_n_locations")


# --- Violin plot ---
font_size_summary = 22

violin_labels = vcat([fill(results_n_locations[i], length(results_all_times[i]))
                      for i in 1:length(results_n_locations)]...)
violin_times  = vcat(results_all_times...)

plt_timing = violin(violin_labels, violin_times,
    fillcolor      = :lightgrey,
    linecolor      = :grey,
    lw             = 2,
    label          = "Distribution",
    xlabel         = "Number of VUF-constrained locations",
    ylabel         = "OPF solve time (s)",
    size           = (1200, 600),
    framestyle     = :box,
    margin         = 10mm,
    left_margin    = 15mm,
    xticks         = all_n_locations,
    legend         = :topleft,
    fontfamily     = "Courier",
    xtickfontsize  = font_size_summary,
    ytickfontsize  = font_size_summary,
    xguidefontsize = font_size_summary,
    yguidefontsize = font_size_summary,
    legendfontsize = font_size_summary - 4,
    # ylim           = (0, 300)
)

plot!(plt_timing, results_n_locations, med_times,
      color     = :black,
      lw        = 2,
      linestyle = :dash,
      label     = false)

scatter!(plt_timing, results_n_locations, med_times,
         color       = :black,
         markersize  = 8,
         markershape = :circle,
         label       = "Median")

for (i, m) in enumerate(med_times)
    annotate!(plt_timing, results_n_locations[i], m + 14,
              text("$(round(m, digits=1))", font(font_size_summary-8, "Courier"), :center))
end

savefig(plt_timing, "../results/scalability_tests/scalability_constrained_locations.png")
savefig(plt_timing, "../results/scalability_tests/scalability_constrained_locations.pdf")
display(plt_timing)
println("Plot saved to results/scalability_tests/scalability_constrained_locations.png / .pdf")
