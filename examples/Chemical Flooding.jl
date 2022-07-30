using Revise, AnalyticalEOR, Plots, Roots, ForwardDiff


kr = RelPerms(swr=0.2,
                sor=0.2,
                krw0=0.2,
                kro0=0.5,
                nw= 2.0,
                no=3.5)

μw = 1.0
μo = 5.0

sj = 0.
si = 0.80

wf = WaterFlooding(si=si, sj=sj, kr=kr, μw=μw, μo=μo)

wf_sol = solve_wf(wf, si, sj)
plot_fractional_flow(wf, wf_sol)

s, v = get_saturation_speeds(wf_sol)

plot_sat_profile(wf_sol, 0.4)
animate_sat_profile(wf_sol, 1)

## Chemical Flooding
kr1 = RelPerms(swr=0.2,
                sor=0.2,    
                krw0=0.5,
                kro0=1.0,
                nw=3.0,
                no=2.0)

kr2 = RelPerms(swr=0.2,
                sor=0.2,    
                krw0=0.2,
                kro0=1.0,
                nw=3.0,
                no=2.0)
            
kr_list = [kr2, kr1]


ChemicalFlooding(si = 0.2,
                 sj = 0.8,
                 krs = [kr1 kr2],
                 μw = 1.0,
                 μo = 5.0
                 )

f_list = [s -> fractional_flow(s, kr, μw, μo) for kr in kr_list]
df_list = [s -> ForwardDiff.derivative.(s -> f(s), s) for f in f_list]
ddf_list = [s -> ForwardDiff.derivative.(s -> df(s), s) for df in df_list]

D_list = [0.3, 0.1]

sol = []

sj2 = 0.8

si = 0.2
sj = 0.8

for (i,D) in enumerate(D_list)
    if i == length(D_list)
        sj2, pf_sol = solve_chemical_flooding(f_list[i], df_list[i], D_list[i], sj2, wf.f)
        wf_sol = solve_wf(wf, si, sj2)
        sol = vcat(sol, pf_sol, wf_sol)
    else
        sj2, pf_sol = solve_chemical_flooding(f_list[i], df_list[i], D_list[i], sj2, f_list[i+1])
        sol = vcat(sol, pf_sol)
        
    end
end



plot_sat_profile(sol, 0.3)




D = 0.3
f = f_list[1]
f3 = f_list[2]
df = df_list[1]

Δfc(s, D) = f(s) / (s - D)
si = 0.2
sj = 0.8

shock_is_admissible(Δfc, f, -D, sj)

s1 = find_shock_spreading(Δfc, df, -D, sj)
m = Δfc(s1, -D)
y(s) = m .* (s .+ D)

plot(s, f(s), label=false)
scatter!([-D s1 sj], [0 f(s1) f(sj)], labels=false)
plot!([-D, s1], [0, y(s1)])


s2 = find_cross_f(y, f3, -D, s1)


plot(s, f(s), label=false)
plot!(s, f3(s), label=false)
scatter!([-D s1 s2 sj], [0 f(s1) f3(s2) f(sj)], labels=false)
plot!([-D, s1], [0, y(s1)])

D2 = 0.1
f4 = wf.f

Δfc2(s, D) = f3(s) / (s - D)
m2 = Δfc2(s2, -D2)
y2(s) = m2 .* (s .+ D2)


plot(s, f(s), label=false)
plot!(s, f3(s), label=false)
plot!(s, f4(s), label=false)
scatter!([-D2 -D s1 s2 sj], [0 0 f(s1) f3(s2) f(sj)], labels=false)
plot!([-D, s1], [0, y(s1)])

plot!([-D2, s2], [0, y2(s2)], labels=false)


find_cross_f(y2, f4, -D, s2)




