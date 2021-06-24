using AnalyticalEOR

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


plot_fw(wf)
plot_sw_profile(wf, 0.05, 0.2, 0.3, 0.5)








