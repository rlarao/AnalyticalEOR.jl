using AnalyticalEOR
using Plot
using LaTeXStrings


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
s = collect(0:0.01:1)

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
coo_mw1 = 0.01066
coo_mw2 = 0.00954

ω1 = (coo_mw1 - coo_ww) / (coo_ow - coo_ww)
ω2 = (coo_mw2 - coo_ww) / (coo_ow - coo_ww)

sor(ω) = (1-ω) * kr_ww.sor + ω * kr_ow.sor
krw0(ω) = (1-ω) * kr_ww.krw0 + ω * kr_ow.krw0
kro0(ω) = (1-ω) * kr_ww.kro0 + ω * kr_ow.kro0
no(ω) = (1-ω) * kr_ww.no + ω * kr_ow.no
nw(ω) = (1-ω) * kr_ww.nw + ω * kr_ow.nw


kr_mw1 = RelPerms(swr=0.11,
                sor=sor(ω1),
                krw0=krw0(ω1),
                kro0=kro0(ω1),
                nw=nw(ω1),
                no=no(ω1))

kr_mw2 = RelPerms(swr=0.11,
                sor=sor(ω2),
                krw0=krw0(ω2),
                kro0=kro0(ω2),
                nw=nw(ω2),
                no=no(ω2))



krw_mw1(s) = water_rel_perm(s, kr_mw1)
kro_mw1(s) = oil_rel_perm(s, kr_mw1)

krw_mw2(s) = water_rel_perm(s, kr_mw2)
kro_mw2(s) = oil_rel_perm(s, kr_mw2)



plot!(s, krw_mw.(s), color=:green, lw=3, alpha=0.7 )
plot!(s, kro_mw.(s), color=:green, lw=3, alpha=0.7, ls=:dash )

μw = 1.0
μo = 5.0

fww(s) = fractional_flow(s, kr_ww, μw, μo)
fow(s) = fractional_flow(s, kr_ow, μw, μo)
fmw1(s) = fractional_flow(s, kr_mw1, μw, μo)
fmw2(s) = fractional_flow(s, kr_mw2, μw, μo)

si = 0.11
sj = 1 - kr_mw2.sor

σq = 5.264485
σp = 0.83538

dfww(s) = fw_derivative(s, kr_ww, μw, μo)
dfow(s) = fw_derivative(s, kr_ow, μw, μo)
dfmw1(s) = fw_derivative(s, kr_mw1, μw, μo)
dfmw2(s) = fw_derivative(s, kr_mw2, μw, μo)


# * Step 1 - Calculate vsj
v1(s) = fmw1.(s) ./ (s + σp)
s1 = fzero(s -> v1(s) - dfmw1(s), 0.8)
vs1 = v1(s1)

# * Step2 - Calculate s1 and vs1
s2 = fzero(s -> (s + σp) * vs1 - fow(s), 0.5)
v2(s) = (fow(s2)) /(s2 + σq)

# * Step 5 Initial shock
Δfow(s) = (fow(s) - fow(si)) / (s - si)
s3  = fzero(s -> dfow(s) - Δfow.(s), 0.3)
vs3 = (fow(s3) - fow(si)) / (s3 - si)


# * Step 6 Tracer
Δc3 = fow(s4) / s4

s = collect(0:0.01:1)

plot(s, fmw2.(s), color=:green, lw=3, alpha=0.7, xlim=(0,1), legend=false)
plot!(s, fow.(s), color=:orangered3, lw=3, alpha=0.7)
# plot!(s, fmw1.(s), color=:green, lw=3, alpha=0.7)

plot!([sj s1 s2 s3 si], [fmw2(sj) fmw2(s1) fow(s2) fow(s3) fow(si)], seriestype=:scatter, ms=4,color=:black)
plot!([-σp, s1], [0, fmw2(s1)], color=:black, ls=:dash, lw=1)
plot!([si, s3], [0, fow(s3)], color=:black, ls=:dash, lw=1)
plot!([0, s3], [0, fow(s4)], color=:black, ls=:dash, lw=1)


plot!(xlabel="Water saturation, s", ylabel="Water fractional flow, f", legend=false,
xlim=(0, 1),  size=(450,420))


# plot!(xlim=(0.65, 0.91), ylim=(0.85,1.005), aspect_ratio=:equal, xlabel="Water saturation, s" )
# plot!(ann=(0.66, 0.885, L"\mathcal{W_2}"), annotationfontsize=12)
# plot!(ann=(0.66, 0.92, L"\mathcal{W_3}"), annotationfontsize=12)

plot!(ann=(0.89, 0.96, "J"))
plot!(ann=(0.76, 0.925, "1"))
plot!(ann=(0.4, 0.8, "2"))

plot!(ann=(0.11, 0.05, "I"))

plot!(ann=(0.05, 0.18, L"\mathcal{W_1}"), annotationfontsize=12)
plot!(ann=(0.05, 0.6, L"\mathcal{W_2}"), annotationfontsize=12)
# plot!(ann=(0.05, 0.65, L"\mathcal{W_3}"), annotationfontsize=12)


plot!([0.45, 0.5], [fmw1(0.45), fmw1(0.45)], color=:black, lw=0.5, label=false)
# plot!([0.39, 0.5], [fmw(0.39), fmw(0.39)], color=:black, lw=0.5, label=false)
plot!([0.24, 0.5], [fow(0.24), fow(0.24)], color=:black, lw=0.5, label=false)

plot!(ann=(0.52, fmw1(0.46), L"\textbf{c}^P, \ \omega \equal 0.10"), annotationfontsize=10, annotationhalign=:left)
# plot!(ann=(0.52, fmw(0.40), L"\textbf{c}^P, \ \omega \equal 0.04"), annotationfontsize=10, annotationhalign=:left)
plot!(ann=(0.52, fow(0.25), L"\textbf{c}^I, \ \omega \equal 1"), annotationfontsize=10, annotationhalign=:left)
plot!(title="b) Glycine 1 wt%")

rectangle(x1, x2, y1, y2) = Shape([x1,x2,x2,x1], [y1,y1, y2, y2])
plot!(rectangle(0.65, 0.91, 0.85, 1.005), opacity=0.1, color=:gray)

savefig("fractional_flow_sol_gly1.svg")



t= 0.35

n = 10

s = vcat(sj, collect(range(sj,s1, length=n)),
        [s2, s2, si, si ])
x = Vector{Float64}(undef, length(s))

x[1] = 0.0
x[2:n+1] = dfmw1(s[2:n+1]) * t
x[n+2] = vs1 * t 
x[n+3] = vs3 * t 
x[n+4] = vs3 * t 
x[n+5] = 1.0 


plot(y, s_wf, lw=2, color=:orangered3, fill=(0,0.03),ls=:dash, label="Waterflooding", xlim=(0,1))
plot!(x, s, lw=3, color=:forestgreen, fill=(0,0.05), label="Glycine Inj")
plot!(xlabel="Dimensionless distance", ylabel="Water saturation, s", xlim=(0,1), ylim=(0,1), legend=(0.8, 0.6),
size=(550,420), title="a) Glycine 1 wt%", titlelocation=:left)
hline!([sj, s1, s2, si], color=:black, lw=0.5, ls=:dash, label=false)


x[1] = 0


vs1
vs3
sj
si
s1
s2

function get_saturation_g1(t)
        sj = 0.8944005607966523
        si = 0.11
        s1 = 0.7518598702831422
        s2 = 0.40538250237934736

        vs1 = 0.6095418616330046
        vs3 = 2.5605429823967616

        n = 10

        s = vcat(sj, collect(range(sj,s1, length=n)),
                [s2, s2, si, si ])
        x = Vector{Float64}(undef, length(s))

        x[1] = 0.0
        x[2:n+1] = dfmw1(s[2:n+1]) * t
        x[n+2] = vs1 * t 
        x[n+3] = vs3 * t 
        x[n+4] = vs3 * t 
        x[n+5] = 1.0 

        return min.(x, 1), s
end








plot!(ann=(0.5, sj+0.03, L"s_J"))
plot!(ann=(0.5, s1+0.01, L"s_1"))
plot!(ann=(0.5, s2+0.02, L"s_2"))
plot!(ann=(0.5, si+0.01, L"s_I"))
plot!(ann=(0.87, 0.95, L"t \equal 0.35 \ PV"))

x3 = Δc3 *t
x2 = vs1 * t
x1 = vs1 * t
# plot!(ann=(x1+0.03, 0.25, L"\mathcal{W_3}"), annotationfontsize=10)
plot!(ann=(x2+0.03, 0.25, L"\mathcal{W_2}"), annotationfontsize=10)
plot!(ann=(x3+0.03, 0.25, L"\mathcal{W_1}"), annotationfontsize=10)


vline!([x3 x2 x1], color=:black, lw=0.5, ls=:dash, label=false)


savefig("sw_profile_gly1.svg")

v5(s) = (fow.(s) - fow.(si)) ./ (s - si)
s5 = fzero(s -> v5(s) - dfow(s), 0.5) 
vs5 = v5(s5)

s = collect(0:0.01:1)
plot(s, fow.(s), color=:orangered3, lw=3, alpha=0.7, xlim=(0,1), legend=false)
plot!(xlabel="Water saturation, s", ylabel="Water fractional flow, f", legend=false,
xlim=(0, 1),  size=(450,420))

sj2 = 1-kr_ow.sor
plot!([sj2 s5 si], [fow(sj2) fow(s5) fow(si)], seriestype=:scatter, ms=4,color=:black)
plot!([si, s5], [0, fow(s5)], color=:black, ls=:dash, lw=1)

plot!(ann=(0.84, 0.96, "J"))
plot!(ann=(0.4, 0.8, "1"))
plot!(ann=(0.11, 0.05, "I"))
plot!([0.24, 0.4], [fow(0.24), fow(0.24)], color=:black, lw=0.5, label=false)
plot!(ann=(0.4, fow(0.25), L"\textbf{c}^I, \ \omega \equal 1"), annotationfontsize=10, annotationhalign=:left)
plot!(title="a) Waterflooding")
savefig("Waterflooding_fw.svg")

m = 50
s_wf = vcat(sj,
        collect(range(sj2, s5, length=m)),
        [s5, si, si],
)

y = Vector{Float64}(undef, length(s_wf))

y[1] = 0.
y[2:1+m] = dfow(s_wf[2:1+m]) *t
y[m+2] = vs5 *t
y[m+3] = vs5 *t
y[m+4] = 1

plot!(xlabel="Dimensionless distance, x", ylabel="Water saturation, w", xlim=(0,1), ylim=(0,1), legend=false,
size=(550,400))