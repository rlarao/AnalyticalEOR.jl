using AnalyticalEOR, Roots


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

si = 0.55
sj = 0.80

D = 0.2

pf = solve_polymerflood(si, sj, kr_ow, kr_ww, D, μw, μo)
tracer = solve_tracer(pf)

plot_sw_profile(pf, tracer, 0.2)

animate_sw_profile(pf, tracer)