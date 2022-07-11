using Revise, AnalyticalEOR, Plots, Roots, ForwardDiff


kr = RelPerms(swr=0.2,
                sor=0.2,
                krw0=0.2,
                kro0=0.5,
                nw= 3.5,
                no=3.5)

μw = 1.0
μo = 5.0


sj = 0.80
si = 0.20

wf = solve_waterflooding(si, sj, kr, μw, μo)
plot_fw(wf)


sj = 0.20
si = 0.80

wf = solve_waterflooding(si, sj, kr, μw, μo)
plot_fw(wf)







fw(sw) = fractional_flow(sw, kr, μw, μo)
fo(sw) = 1 .- fw(sw)
dfo(sw) = ForwardDiff.derivative.(sw -> fo(sw), sw)
dfw(sw) = fw_derivative(sw, kr, μw, μo)
d2fw(sw) = fw_derivative2(sw, kr, μw, μo)
d3fw(sw) = fw_derivative3(sw, kr, μw, μo)
Δfw(sw) = (fw(sw) - fw(si)) / (sw - si) 

zeros = find_zeros(s -> d2fw(s), si, sj)

d2fw.(zeros)
d3fw.(zeros)

zeros[2]  ≈ sj

s = collect(si:0.01:sj)

plot(s, fw(s), lims=(0,1))
plot(s, fo(s), lims=(0,1))
plot(s, dfo(s), xlims=(0,1))

plot(s, dfw(s), xlims=(0,1))
plot(s, -d2fw(s), xlims=(0,1), color=:red)
plot(s, d3fw(s), xlims=(0,1), color=:orange)
plot!(s, -d2sw(s), xlims=(0,1), color=:green)

fw_derivative2(s, kr, μw, μo)

S̃  = find_zeros(sw -> dfw(sw) - Δfw.(sw), si, sj)

plot(s, dfw(s))