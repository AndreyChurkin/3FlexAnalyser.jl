"""
Visualise scalability test results: OPF solve time vs severity of VUF constraints

Reads the CSV with all OPF times produced by scalability_test_vuf_severity.jl and creates
the violin plot without re-running any OPF simulations.

Input:  results/scalability_tests/scalability_vuf_severity_all_times.csv (or a backup copy)
Output: results/scalability_tests/scalability_vuf_severity.png / .pdf

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

# csv_path = "../results/scalability_tests/scalability_vuf_severity_all_times.csv"
csv_path = "../results/scalability_tests/scalability_vuf_severity_all_times_2026_05_20.csv"


# ============================================================


# --- Load data ---
df = CSV.read(csv_path, DataFrame)
sort!(df, :vuf_threshold)

all_n_thresholds = sort(unique(df.vuf_threshold), rev=true)

# Reconstruct time vectors per n_thresholds (same structure as 'results_all_times' in the test script):
results_n_thresholds   = all_n_thresholds .* 100
results_all_times = [df[df.vuf_threshold .== n, :time_s] for n in all_n_thresholds]

med_times = [isempty(t) ? NaN : median(t) for t in results_all_times]

println("Loaded $(nrow(df)) OPF times for $(length(all_n_thresholds)) VUF thresholds: $all_n_thresholds")


# --- Violin plot ---
font_size_summary = 22

n_thresholds = length(results_n_thresholds)
# tick_labels  = [string(round(v, digits=2)) for v in results_n_thresholds]
tick_labels  = [isodd(i) ? string(round(v, digits=2)) : "" for (i,v) in enumerate(results_n_thresholds)]

# Use integer positions 1…N so violins are evenly spaced regardless of threshold spacing:
violin_labels = vcat([fill(i, length(results_all_times[i])) for i in 1:n_thresholds]...)
violin_times  = vcat(results_all_times...)

plt_timing = violin(violin_labels, violin_times,
    fillcolor      = :lightgrey,
    linecolor      = :grey,
    lw             = 2,
    label          = "Distribution",
    xlabel         = "VUF limit, %",
    ylabel         = "OPF solve time, s",
    size           = (1200, 600),
    framestyle     = :box,
    margin         = 10mm,
    left_margin    = 15mm,
    xticks         = (1:n_thresholds, tick_labels),
    legend         = :topleft,
    fontfamily     = "Courier",
    xtickfontsize  = font_size_summary,
    ytickfontsize  = font_size_summary,
    xguidefontsize = font_size_summary,
    yguidefontsize = font_size_summary,
    legendfontsize = font_size_summary - 4,

    ylim           = (0, 300)
)

plot!(plt_timing, 1:n_thresholds, med_times,
      color     = :black,
      lw        = 2,
      linestyle = :dash,
      label     = false)

scatter!(plt_timing, 1:n_thresholds, med_times,
         color       = :black,
         markersize  = 8,
         markershape = :circle,
         label       = "Median")

for (i, m) in enumerate(med_times)
    annotate!(plt_timing, i, m + 14, text("$(round(m, digits=1))", font(font_size_summary-8, "Courier New Bold"), :center))
end


savefig(plt_timing, "../results/scalability_tests/scalability_vuf_severity.png")
savefig(plt_timing, "../results/scalability_tests/scalability_vuf_severity.pdf")
display(plt_timing)
println("Plot saved to results/scalability_tests/scalability_vuf_severity.png / .pdf")
