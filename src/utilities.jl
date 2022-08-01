function RecoveryFactor(sol, PVs) 
    s, v = get_saturation_speeds(sol)
    OOIP = 1 - s[end]

    x(t) = min.(v .* t, 1.0)

    RF = [1 - CumTrapz(x(t), s) / OOIP for t in PVs]

end

function CumTrapz(x::T, s::T) where {T <: AbstractVector}
    # Initialize Output
    out = similar(x)
    out[1] = 0

    # Iterate over arrays
    n = length(x)
    for i in 2:n
      out[i] = out[i-1] + 0.5*(x[i] - x[i-1])*(s[i] + s[i-1])
    end

    return 1 - out[end]
end