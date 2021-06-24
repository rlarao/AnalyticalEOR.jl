using AnalyticalEOR, Plots

swr = 0.2
sor = 0.2

kro0 = 1.0
krw0 = 0.5
nw = 3.0
no = 1.0

kr = RelPerms(swr, sor, krw0, kro0, nw, no)

μw = 1.0
μo = 3.0

si = 0.20
sj = 0.80


wf = solve_waterflooding(si, sj, kr, μw, μo)


tracer = solve_tracer(wf)

plot_tracer_fw(wf, tracer)


plot_sw_profile(wf, tracer, 0.1)
animate_sw_profile(wf, tracer)