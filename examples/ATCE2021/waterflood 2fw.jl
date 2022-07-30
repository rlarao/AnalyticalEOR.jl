using AnalyticalEOR
using Plots


kr_ww = RelPerms(swr=0.11,
                sor=0.10,
                krw0=0.2,
                kro0=0.8,
                nw=3.,
                no=2.)

kr_ow = RelPerms(swr=0.11,
                sor=0.16,
                krw0=0.4,
                kro0=0.5,
                nw=2.,
                no=3.)

krw_ww(s) = water_rel_perm(s, kr_ww)
kro_ww(s) = oil_rel_perm(s, kr_ww)

krw_ow(s) = water_rel_perm(s, kr_ow)
kro_ow(s) = oil_rel_perm(s, kr_ow)


plot(s, krw_ww.(s), color=:dodgerblue2, lw=3, alpha=0.7)
plot!(s, kro_ww.(s), color=:dodgerblue2, lw=2.5, alpha=0.7, ls=:dash)
plot!(lims=(0,1), xlabel="Water saturation, s", ylabel="Relative permeability",
size=(450, 400), legend=false)

plot!(s, krw_ow.(s), color=:orangered3, lw=3, alpha=0.7 )
plot!(s, kro_ow.(s), color=:orangered3, lw=3, alpha=0.7, ls=:dash )


coo_ww = 0.002098
coo_ow = 0.0818417
coo_mw = 0.01664

ω = (coo_mw - coo_ww) / (coo_ow - coo_ww)

sor(ω) = (1-ω) * kr_ww.sor + ω * kr_ow.sor
krw0(ω) = (1-ω) * kr_ww.krw0 + ω * kr_ow.krw0
kro0(ω) = (1-ω) * kr_ww.kro0 + ω * kr_ow.kro0
no(ω) = (1-ω) * kr_ww.no + ω * kr_ow.no
nw(ω) = (1-ω) * kr_ww.nw + ω * kr_ow.nw


kr_mw = RelPerms(swr=0.11,
                sor=sor(ω),
                krw0=krw0(ω),
                kro0=kro0(ω),
                nw=nw(ω),
                no=no(ω))

krw_mw(s) = water_rel_perm(s, kr_mw)
kro_mw(s) = oil_rel_perm(s, kr_mw)



plot!(s, krw_mw.(s), color=:green, lw=3, alpha=0.7 )
plot!(s, kro_mw.(s), color=:green, lw=3, alpha=0.7, ls=:dash )

μw = 1.0
μo = 5.0


fww(s) = fractional_flow(s, kr_ww, μw, μo)
fow(s) = fractional_flow(s, kr_ow, μw, μo)
fmw(s) = fractional_flow(s, kr_mw, μw, μo)




si = 0.11
sj = 1 - 0.16

fww(sj)


σj = σ[2] - 1
σq = σ[2] - 1
σp = σ[4] - 1


dfww(s) = fw_derivative(s, kr_ww, μw, μo)
dfow(s) = fw_derivative(s, kr_ow, μw, μo)
dfmw(s) = fw_derivative(s, kr_mw, μw, μo)


# * Step 1 - Calculate vsj W3
v1(s) = fww.(s) ./ (s + σj)
s1 = fzero(s -> v1(s) - dfww(s), sj)
vsj = v1(sj)


# * Step2 - Calculate s1 and vs1
s1 = fzero(s -> (s + σj) * vsj - fow(s), 0.5)
vs1 = v2(s1)


# * Step3 - Calculate s1 and vs1 W2
v2(s) = fow.(s) ./ (s + σp)
s2 = fzero(s -> v2(s) - fow(s), 0.45)

vs2 = (fow(s1) - fow(s4)) / (s1 - s4)

# * Step 5 Initial shock
Δfow(s) = (fow(s) - fow(si)) / (s - si)
s4  = fzero(s -> dfow(s) - Δfow.(s), 0.3)
vs4 = Δfow(s4)

# * Step 6 Tracer
Δc3 = fow(s5) / s5

s = collect(0:0.01:1)
plot(s, fww.(s), color=:dodgerblue2, lw=3, alpha=0.7)
plot!(s, fow.(s), color=:orangered3, lw=3, alpha=0.7)
plot!(xlim=(0,1))


plot!([sj s1 s4 si], [fww(sj) fow(s1) fow(s4) fow(si)], seriestype=:scatter, ms=3,
        color=:black)
plot!([-σj, sj], [0, fww(sj)], color=:black, ls=:dash, lw=1)
plot!([-σp, s1], [0, fow(s1)], color=:black, ls=:dash, lw=1)
plot!([s4, s1], [fow(s4), fow(s1)], color=:black, ls=:dash, lw=1)
plot!([si, s4], [0, fow(s4)], color=:black, ls=:dash, lw=1)
plot!([0, s5], [0, fow(s5)], color=:black, ls=:dash, lw=1)

plot!(xlabel="Water saturation, s", ylabel="Water fractional flow, f", legend=false,
xlim=(0, 1),  size=(450,400))
#  xlim=(0., 1), ylim=(0, 1))

plot!(ann=(0.84, 0.96, "J"))
plot!(ann=(0.58, 0.94, "1"))
plot!(ann=(0.42, 0.73, "2"))
# plot!(ann=(0.41, 0.71, "4"))
plot!(ann=(0.11, 0.05, "I"))
plot!(ann=(0.09, 0.26, L"\mathcal{W_1}"), annotationfontsize=12)
plot!(ann=(0.09, 0.77, L"\mathcal{W_2}"), annotationfontsize=12)
plot!(ann=(0.09, 0.97, L"\mathcal{W_3}"), annotationfontsize=12)


plot!([0.6, 0.7], [fww(0.6), fww(0.6)], color=:black, lw=0.5)
plot!([0.29, 0.7], [fow(0.29), fow(0.29)], color=:black, lw=0.5)

plot!(ann=(0.72, fww(0.6), L"\textbf{c}^J, \ \omega \equal 1"), annotationfontsize=10, annotationhalign=:left)
plot!(ann=(0.72, fow(0.29), L"\textbf{c}^P, \ \omega \equal 0"), annotationfontsize=10, annotationhalign=:left)

savefig("one_fractional_flow.svg")

t= 0.35
n = 10
s = vcat([sj, sj, s1, s1 ,s1, s4, s4, si , si])
x = Vector{Float64}(undef, length(s))

x[1] = 0.0
x[2] = vsj * t
x[3] = vsj * t 
x[4] = vs1 * t 
x[5] = vs2 * t 
x[6] = vs2 * t 
x[7] = vs4 * t 
x[8] = vs4 * t 
x[9] = 1.0

x[5:5+n-1] = dfow(s[5:5+n-1]) * t
x[n+4] = vs4 * t
x[n+5] =  vs4 * t
x[n+6 ] = 1.0

x
s[5:5+n-1]

plot(x, s, lw=3, color=:dodgerblue2, fill=(0,0.07))
plot!(xlabel="Dimensionless distance, x", ylabel="Water saturation, w", xlim=(0,1), ylim=(0,1), legend=false,
size=(550,400))
hline!([sj, s1, s4 , si], color=:black, lw=0.5, ls=:dash)

plot!(ann=(0.5, sj+0.03, L"s_J"))
plot!(ann=(0.5, s1+0.01, L"s_1"))
plot!(ann=(0.5, s4-0.03, L"s_2"))
plot!(ann=(0.5, si+0.03, L"s_I"))
plot!(ann=(0.87, 0.95, L"t \equal 0.35 \ PV"))

plot!(ann=(x1+0.03, 0.25, L"\mathcal{W_3}"), annotationfontsize=10)
plot!(ann=(x2+0.03, 0.25, L"\mathcal{W_2}"), annotationfontsize=10)
plot!(ann=(x3+0.03, 0.25, L"\mathcal{W_1}"), annotationfontsize=10)
savefig("sw_profile_2fw.svg")


vline!([x3 x2 x1], color=:black, lw=0.5, ls=:dash)
x3 = Δc3 *t
x2 = vs2 * t
x1 = vsj * t







# * Step 1 - Calculate vsj W3
v1(s) = fww.(s) ./ (s + σj)
s1 = fzero(s -> v1(s) - dfww(s), sj)
vsj = v1(sj)


# * Step2 - Calculate s1 and vs1
s1 = fzero(s -> (s + σj) * vsj - fow(s), 0.5)
vs1 = v2(s1)


# * Step3 - Calculate s1 and vs1 W2
v2(s) = fow.(s) ./ (s + σp)
s2 = fzero(s -> v2(s) - fow(s), 0.45)

vs2 = (fow(s1) - fow(s4)) / (s1 - s4)

# * Step 5 Initial shock
Δfow(s) = (fow(s) - fow(si)) / (s - si)
s4  = fzero(s -> dfow(s) - Δfow.(s), 0.3)
vs4 = Δfow(s4)

# * Step 6 Tracer
Δc3 = fow(s5) / s5

s = collect(0:0.01:1)
plot(s, fww.(s), color=:dodgerblue2, lw=3, alpha=0.7)
plot!(s, fow.(s), color=:orangered3, lw=3, alpha=0.7)
plot!(xlim=(0,1))


plot!([sj s1 s4 si], [fww(sj) fow(s1) fow(s4) fow(si)], seriestype=:scatter, ms=3,
        color=:black)
plot!([-σj, sj], [0, fww(sj)], color=:black, ls=:dash, lw=1)
plot!([-σp, s1], [0, fow(s1)], color=:black, ls=:dash, lw=1)
plot!([s4, s1], [fow(s4), fow(s1)], color=:black, ls=:dash, lw=1)
plot!([si, s4], [0, fow(s4)], color=:black, ls=:dash, lw=1)
plot!([0, s5], [0, fow(s5)], color=:black, ls=:dash, lw=1)

plot!(xlabel="Water saturation, s", ylabel="Water fractional flow, f", legend=false,
xlim=(0, 1),  size=(450,400))
#  xlim=(0., 1), ylim=(0, 1))

plot!(ann=(0.84, 0.96, "J"))
plot!(ann=(0.58, 0.94, "1"))
plot!(ann=(0.42, 0.73, "2"))
# plot!(ann=(0.41, 0.71, "4"))
plot!(ann=(0.11, 0.05, "I"))
plot!(ann=(0.09, 0.26, L"\mathcal{W_1}"), annotationfontsize=12)
plot!(ann=(0.09, 0.77, L"\mathcal{W_2}"), annotationfontsize=12)
plot!(ann=(0.09, 0.97, L"\mathcal{W_3}"), annotationfontsize=12)


plot!([0.6, 0.7], [fww(0.6), fww(0.6)], color=:black, lw=0.5)
plot!([0.29, 0.7], [fow(0.29), fow(0.29)], color=:black, lw=0.5)

plot!(ann=(0.72, fww(0.6), L"\textbf{c}^J, \ \omega \equal 1"), annotationfontsize=10, annotationhalign=:left)
plot!(ann=(0.72, fow(0.29), L"\textbf{c}^P, \ \omega \equal 0"), annotationfontsize=10, annotationhalign=:left)

savefig("one_fractional_flow_rarefaction.svg")

t= 0.3
n = 10
s = vcat([sj, sj, s1, s1], collect(range(s1,s4, length=n)), [ s4, si , si])
x = Vector{Float64}(undef, length(s))

x[1] = 0.0
x[2] = vsj * t
x[3] = vsj * t 
x[4] = vs1 * t 
x[5] = vs2 * t 
x[5:5+n-1] = dfow(s[5:5+n-1]) * t
x[n+5] = vs4 * t 
x[n+6] = vs4 * t 
x[n+7] = 1.0





plot(x, s, lw=3, color=:dodgerblue2, fill=(0,0.07))
plot!(xlabel="Dimensionless distance, x", ylabel="Water saturation, w", xlim=(0,1), ylim=(0,1), legend=false,
size=(550,400))
hline!([sj, s1, s4 , si], color=:black, lw=0.5, ls=:dash)

plot!(ann=(0.5, sj+0.03, L"s_J"))
plot!(ann=(0.5, s1+0.01, L"s_1"))
plot!(ann=(0.5, s4-0.03, L"s_2"))
plot!(ann=(0.5, si+0.03, L"s_I"))
plot!(ann=(0.87, 0.95, L"t \equal 0.35 \ PV"))

plot!(ann=(x1+0.03, 0.25, L"\mathcal{W_3}"), annotationfontsize=10)
plot!(ann=(x2+0.03, 0.25, L"\mathcal{W_2}"), annotationfontsize=10)
plot!(ann=(x3+0.03, 0.25, L"\mathcal{W_1}"), annotationfontsize=10)
savefig("sw_profile_2fw_rarecation.svg")


vline!([x3 x2 x1], color=:black, lw=0.5, ls=:dash)
x3 = Δc3 *t
x2 = vs2 * t
x1 = vsj * t
x[1:]