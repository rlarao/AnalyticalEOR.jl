# function solve_waterflooding(si, sj, kr, μw, μo)
#     if si <= sj
#         fw(sw) = fractional_flow(sw, kr, μw, μo)
#         dfw(sw) = fw_derivative(sw, kr, μw, μo)
#         Δfw(sw) = (fw(sw) - fw(si)) / (sw - si) 

#         S̃  = find_zeros(sw -> dfw(sw) - Δfw.(sw), si, sj)
#         if isempty(S̃)
#             S̃ = NaN
#         else
#             S̃ = S̃[end]
#         end
#             WaterFlooding(S̃, si, sj, kr, μw, μo)
#     else
#         fw(sw) = 1-fractional_flow(sw, kr, μw, μo)
#         dfw(sw) = fw_derivative(sw, kr, μw, μo)
#         Δfw(sw) = (1-fw(sw) - fw(sj)) / (1-sw - sj) 

#         S̃  = find_zeros(sw -> dfw(sw) - Δfw.(sw), si, sj)
#         if isempty(S̃)
#             S̃ = NaN
#         else
#             S̃ = S̃[end]
#         end
#             WaterFlooding(S̃, si, sj, kr, μw, μo)
#     end
# end

function get_saturation_speeds(sol)
    s = []
    v = []

    for (s1, s2, wave, speed, f) in sol

        if wave == :spreading
            s_range = collect(range(s1, s2, length=20))
            s = vcat(s, s_range)
            v = vcat(v, speed(s_range))
        elseif wave == :shock
            s = vcat(s, [s1, s2])
            v = vcat(v, [speed, speed])

        end
    end

    insert!(s,1, s[1])
    insert!(v,1, 0.)

    push!(s, s[end])
    push!(v, 10000.)

    s = convert(Array{Float64, 1}, s)
    v = convert(Array{Float64, 1}, v)

    return s, v
end


# function saturation_profile(wf::WaterFlooding, time::Real)
#     si = wf.Si
#     sj = wf.Sj
#     s̃ = wf.S̃

#     fw(sw) = fractional_flow(sw, wf.kr, wf.μw, wf.μo)
#     dfw(sw) = fw_derivative(sw, wf.kr, wf.μw, wf.μo)
#     Δfw(sw) = (fw(sw) - fw(si)) / (sw - si) 

#     if !isnan(s̃)
#         s = vcat([si,si], collect(s̃:0.01:sj), [sj])
#         x = Vector{Float64}(undef, length(s))
#         x[1] = 1.0
#         x[2:3] .= Δfw(s̃) * time
#         x[findall(s .> s̃)] = dfw(s[s .> s̃]) * time
#         x[end] = 0.0
#     else
#         s = vcat([si], collect(si:0.01:sj), [sj])
#         x = Vector{Float64}(undef, length(s))
#         x[1] = 1.0
#         x[2:end-1] = dfw.(s[2:end-1]) * time
#         x[end] = 0.0
#     end

#     return x, s
# end




# struct WaterFlooding{T <: Real}
#     S̃::T
#     Si::T
#     Sj::T
#     RelPerms::RelPerms
#     μw::T
#     μo::T

#     function WaterFlooding(Si::T, Sj::T, krs::RelPerms,  μw::T, μo::T) where {T <: Real}
#         μw, μo = μ
#         df(s) = fw_derivative(sw, krs, μw, μo)    
#         Δfw(sw) = (fw(sw) - fw(si)) / (sw - si) 
#         S̃  = find_zero(sw -> dfw(sw) - Δfw.(sw), (si + 0.0001, sj))
#         new(S̃, Si, Sj, krs, μw, μo)
#     end
# end

function find_sat_shock_spreading(Δf, df, si, sj, swr, sor) 
    min_bound = minimum([swr, 1 - sor]) + 0.001
    max_bound = maximum([swr, 1 - sor]) - 0.001

    roots = find_zeros(s -> df(s) - Δf.(s, si), min_bound, max_bound)
    return deleteat!(roots, findall(x->x≈si || x≈sj, roots))[1]
end

function shock_is_admissible(Δf, f, si, sj)
    m = Δf(sj, si)
    y(s, si) = m .* (s .- si) .+ f(si)

    # min_bound = minimum([si,sj]) + 0.001
    # max_bound = maximum([si,sj]) - 0.001
    # roots = find_zeros(s-> y(s, si) .- f(s), min_bound, max_bound)
    # return isempty(roots)
    s = collect(range(si, sj; length=200))[2:end-1]
    res = (y(s,si) .- f(s)) .* sign(sj - si) .< 0

    return ~any(res)
end

function solve_wf(wf::WaterFlooding, si::Float64, sj::Float64)
    f = wf.f
    df = wf.df
    d²f = wf.d²f

    Δf(s, si) = (f(s) - f(si)) / (s - si)

    s = collect(range(si, sj; length=100))[2:end-1]
    m = Δf(sj, si)
    y(s) = m .* (s .- si)

    if all(d²f(s) .* sign(sj-si) .< 0)
        sol = [(sj, si, :spreading, df, f)]

    # Can get there through a shock?
    elseif shock_is_admissible(Δf, f, si, sj) #all(d²f(s) .* sign(sj - si) .< 0)
        sol = [(sj, si, :shock, Δf(sj, si), f)]

    # Otherwise, got the get there through a spreading wave and a shock
    else
        s1 = find_sat_shock_spreading(Δf, df, si, sj, wf.kr.swr, wf.kr.sor)

        sol = [(sj, s1, :spreading, df, f),
                (s1, si, :shock, Δf(s1, si), f)
        ] 
    end
end


function solve_chemical_flooding(f, df, D, sj, f2, swr, sor)
    Δfc(s, D) = f(s) / (s - D)

    if shock_is_admissible(Δfc, f, -D, sj)
        m = Δfc(sj, - D)
        y(s) = m .* (s .+ D)
        s2 = find_cross_f(y, f2, -D, sj)
        sol = [(sj, s2,:shock, m, y)]
    else
        s1 = find_sat_shock_spreading(Δfc, df, -D, sj, swr, sor)
        m = Δfc(s1, -D)
        shock_line(s) = m .* (s .+ D)

        s2 = find_cross_f(shock_line, f2, -D, s1)
        sol = [
                (sj, s1, :spreading, df, f),
                (s1, s2, :shock, m, shock_line)
                ]
    end
        return s2, sol
end


function find_cross_f(y, f, si, sj)
    min_bound = minimum([si,sj]) + 0.001
    max_bound = maximum([si,sj]) - 0.001

    find_zeros(s-> y(s) .- f(s), min_bound, max_bound)[end]
end



function solve_cf(cf::ChemicalFlooding, si::Float64, sj::Float64)
    sol = []

    f = cf.f
    df = cf.df
    wf = cf.wf
    kr = cf.kr

    sj2 = sj

        for (i,D) in enumerate(cf.D)
            if i == length(cf.D)
                sj2, pf_sol = solve_chemical_flooding(f[i], df[i], D, sj2, wf.f, kr[i].swr, kr[i].sor)
                wf_sol = solve_wf(wf, si, sj2)
                sol = vcat(sol, pf_sol, wf_sol)
            else
                sj2, pf_sol = solve_chemical_flooding(f[i], df[i], D, sj2, f[i+1], kr[i].swr, kr[i].sor)
                sol = vcat(sol, pf_sol)
            end
        end

    return sol
end
