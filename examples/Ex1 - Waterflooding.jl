using Revise, AnalyticalEOR, Plots, Roots, ForwardDiff


kr = RelPerms(swr=0.2,
                sor=0.2,
                krw0=0.2,
                kro0=0.5,
                nw= 2.0,
                no=3.5)

μw = 1.0
μo = 5.0

sj = 0.40
si = 0.20

wf = WaterFlooding(si=si,
                    sj=sj,
                    kr=kr,
                    μw=μw,
                    μo=μo
                    )

wf.sol



f(s) = fractional_flow(s, kr, μw, μo)
df(s) = ForwardDiff.derivative.(s -> f(s), s) 
d2f(s, si, sj) = ForwardDiff.derivative.(s -> df(s), s) * sign(sj-si)


s = collect(range(si, 0.8; length=1000))[2:end-1]
plot(s, f(s))

s = collect(range(si, sj; length=1000))[2:end-1]
plot!(s, f(s), color=:green)



s = collect(range(si, 0.8; length=1000))[2:end-1]
plot(s, df(s))

s = collect(range(si, sj; length=1000))[2:end-1]
plot!(s, df(s), color=:green)


s = collect(range(si, 0.8; length=1000))[2:end-1]
plot(s, d2f(s, si, sj))

s = collect(range(si, sj; length=1000))[2:end-1]
plot!(s, d2f(s, si, sj), color=:green)

plot(s, df(s))
plot(s, d2f(s, si, sj))

all(d2f(s, si, sj) .>=0)
all(d2f(s, si, sj) .<=0)

get_saturation_speeds(wf)

s, v = get_sat_speeds(wf)

plot_sat_profile(wf, 0.2)
animate_sat_profile(wf, 0.5, dt=0.005)
plot_fractional_flow(wf)


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


f(s) = fractional_flow(s, kr_ww, μw, μo)
fow(s) = fractional_flow(s, kr_ow, μw, μo)

s = collect(range(kr.swr, 1-kr.sor; length=101))

plot(s, f(s))
plot!(s, fow(s))
plot!(xlims=(-0.5,1))

D = 0.3
Δfc(s) = f(s) / (s + D)


sj = 0.8

scatter!([sj], [f(sj)])

s = collect(range(-D, sj; length=101))
m = Δfc(sj)
plot!(s, m .* (s .+D))

s = collect(range(si, sj; length=101))
plot(s, f(s) - m .* (s .+ D))

n_zeros= find_zeros(s -> f(s) .- m .* (s .+ D), si, sj)

df(s) = ForwardDiff.derivative.(s -> f(s), s)

if length(n_zeros) == 1
    m = Δfc(sj)
    plot!(s, m .* (s .+D))
    s1 = find_zeros(s -> fow(s) .- m .* (s .+ D), si, sj)[1]
    sol = Dict((sj, s1) => (:shock, m))
    
else
    s1 = find_zeros(s -> -df(s) + Δfc.(s), kr_ww.swr * 1.001, (1 - kr_ww.sor)*0.999)[1]
    m = Δfc(s1)

    s2 = find_zeros(s -> fow(s) .- m .* (s .+ D), kr_ow.swr * 1.011, (1 - kr_ow.sor) * 0.999)[1]

    sol = Dict((sj, s1) => (:spreading, df),
                (s1, s2) => (:shock, m)
                ) 
    
    wf_sol = solve_waterflooding(si, s2, fow)
    
end

wf = WaterFlooding()

si
s2
wf = WaterFlooding(si=si,
                    sj=0.,
                    kr=kr_ow,
                    μw=μw,
                    μo=μo
                    )

plot_fractional_flow(wf)



sol[1][1]
scatter!([s2], [fow(s2)])


dfow(s) = ForwardDiff.derivative.(s -> fow(s), s)
d2f(s, si, sj) = ForwardDiff.derivative.(s -> dfow(s), s) 
d3f(s, si, sj) = ForwardDiff.derivative.(s -> d2f(s, si, sj), s) 



s = collect(range(si, s2; length=1000 ))
all(d2f(s, si, s2) .< 0)
all(d2f(s, si, s2) .> 0)

plot(s, f)
plot(s, d2f(s, si, s2))
plot(s, d3f(s, si, s2))
solve_waterflooding(0.2, 0.45, fow)

