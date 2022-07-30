using Revise, AnalyticalEOR, Plots, Roots

# Inputs
begin
    kr = RelPerms(swr=0.2,
                    sor=0.2,
                    krw0=0.2,
                    kro0=0.5,
                    nw= 2.0,
                    no=3.5)

    μw = 1.0
    μo = 5.0

    si = 0.2
    sj = 0.8

    wf = WaterFlooding(
                        kr=kr,
                        μw=μw,
                        μo=μo
                        )
end

# Solution
begin
    si = 0.2
    sj = 0.8
    sol = solve_wf(wf, si, sj)
    plot_fractional_flow(wf, sol)
end

# Saturation Profile Plot
begin
    t = 0.1
    plot_sat_profile(sol, t)
end

# Saturation Profile Animation
begin
    max_t = 1
    animate_sat_profile(sol, max_t)
end

# Recovery Factor
f = wf.f
df = wf.df
d²f = wf.d²f

Δf(s, si) = (f(s) - f(si)) / (s - si)

s = collect(range(si, sj; length=1000))[2:end-1]

if all(d²f(s) .* sign(sj-si) .< 0)
    sol = [(sj, si, :spreading, df, f)]

# Can get there through a shock?
elseif shock_is_admissible(Δf, f, si, sj) && any(d²f(s) .* sign(sj - si) .< 0)
    m = Δf(sj, si)
    y(s) = m .* (s .- si)
    sol = [(sj, si, :shock, Δf(sj, si), f)]

# Otherwise, got the get there through a spreading wave and a shock
else
    s1 = find_shock_spreading(Δf, df, si, sj, wf.kr.swr, wf.kr.sor)

    m = Δf(s1, si)
    y(s) = m .* (s .- si)

    sol = [(sj, s1, :spreading, df, f),
            (s1, si, :shock, Δf(s1, si), f)
    ] 
end



plot_fractional_flow(wf, sol)


plot(s, d²f(s).*sign(sj - si))
plot(s, f)