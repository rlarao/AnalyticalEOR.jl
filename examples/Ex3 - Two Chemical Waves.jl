using Revise, AnalyticalEOR

# Inputs
begin
    kr = RelPerms(swr=0.2,
                    sor=0.2,
                    krw0=0.2,
                    kro0=0.5,
                    nw= 2.0,
                    no=3.5)

    μw = 1.0
    μo = 5.

    wf = WaterFlooding(si=0.2,
                        sj=0.8,
                        kr=kr,
                        μw=μw,
                        μo=μo
                        )
end

begin
    kr2 = RelPerms(swr=0.2,
                    sor=0.2,    
                    krw0=0.5,
                    kro0=1.0,
                    nw=3.0,
                    no=2.0)

    kr3 = RelPerms(swr=0.2,
                    sor=0.2,    
                    krw0=0.2,
                    kro0=1.0,
                    nw=3.0,
                    no=2.0)

    cf = ChemicalFlooding(
                    wf = wf,
                    krs = [kr3, kr2],
                    D = [0.3, 0.1]
                    )
end


# Solve
begin
    sol = solve_cf(cf)
    plot_fractional_flow(cf, sol)
end

# Saturation Profile Plot
begin
    t = 0.3
    plot_sat_profile(sol, 0.3)
end


# Saturation Profile Animation
begin
    max_t = 2
    animate_sat_profile(sol, max_t)
end