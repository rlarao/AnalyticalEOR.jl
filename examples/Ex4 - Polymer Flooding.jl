using AnalyticalEOR, Roots, Revise


kr_ww = RelPerms(swr=0.2,
                sor=0.2,
                krw0=0.5,
                kro0=1.0,
                nw=3.0,
                no=2.0)

kr_ow = RelPerms(swr=0.2,
                sor=0.2,
                krw0=0.245,
                kro0=0.54,
                nw=2.111,
                no=3.706)

μw = 1.0
μo = 3.0

si = 0.2
sj = 0.80

D = 1.5

pf = solve_polymerflood(si, sj, kr_ow, kr_ww, D, μw, μo)
plot_fw(pf)

wf = solve_waterflooding(0.2, sj, kr_ow, μw, μo)
plot_fw(wf)

plot_sw_profile(wf, 0.5)


tracer = solve_tracer(pf)

plot_sw_profile(pf, tracer, 0.2)

animate_sw_profile(pf, tracer)

plot_fw(pf)
