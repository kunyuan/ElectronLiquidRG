@inline function phase(extT, ninL, noutL, ninR, β) 
    # println(extT)
    tInL, tOutL, tInR, tOutR = extT
    winL, woutL, winR = π * (2ninL + 1) / β, π * (2noutL + 1) / β, π * (2ninR + 1) / β
    woutR = winL + winR - woutL
    return exp(-1im * (tInL * winL - tOutL * woutL + tInR * winR - tOutR * woutR))
end

#Particle-particle diagram, order = 2
#
# kamp,up   kamp2,down
#   |--- u ---|
#   |         |
# up^         ^down
#   |         |
#   |-- KO ---|
#kamp,up  kamp2,down

function _PP(idx, vars, config)
    K, T, Ang, Phi, ExtKidx = vars
    para, kgrid, n = config.userdata
    k = K[1]
    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0
    
    τ1, τ2, τ3 = T[1], T[2], T[3]
    # θ = 2.533
    # φ = 0.0
    
    θ = Ang[1]
    φ = Phi[1]
    extK1 = [sin(0.5θ), cos(0.5θ), 0] * kgrid[ExtKidx[1]]
    extK2 = [-sin(0.5θ), cos(0.5θ), 0] * kgrid[ExtKidx[1]]
    extK3 = [sin(0.5θ)*cos(φ),cos(0.5θ),sin(0.5θ)*sin(φ)] * kgrid[ExtKidx[1]]

    k1 = extK1 + extK2 - k
    q = extK1 - k
    qabs = sqrt(dot(q, q))
    # println(q)
    ϵ = dot(k, k) / (2me) - μ
    G1 = ElectronLiquid.Propagator.green_derive(τ3 - τ1, ϵ, β, 0)
    ϵ = dot(k1, k1) / (2me) - μ
    G21 = ElectronLiquid.Propagator.green_derive(τ3 - τ2, ϵ, β, 0)
    G22 = ElectronLiquid.Propagator.green_derive(τ3 - τ1, ϵ, β, 0)
    Vd = ElectronLiquid.interactionDynamic(para, qabs, τ1, τ2)
    Vs = ElectronLiquid.interactionStatic(para, qabs, τ1, τ2)
    l1 = n[1]
    l2 = n[2]
    l3 = n[3]
    factor = 1.0 / (2π)^3
    factor_legend = sin(θ) / (4π)
    f1 = phase((τ1, τ3, τ2, τ3), l1, l2, l3, β) * factor * factor_legend
    f2 = phase((τ1, τ3, τ1, τ3), l1, l2, l3, β) * factor * factor_legend
    return G1 * (G21 * Vd * f1 + G22 * Vs * f2)
end

#
#   kamp,up       kamp2,down
#       |             |
#       |-- < -\      |
#       |  up   \     |
#       KO      |- u -|
#       |  up   /     | 
#       |-- > -/      |
#       |             |
#   kamp,up       kamp2,down

function _Lver3_direct(idx, vars, config)
    K, T, Ang, Phi, ExtKidx = vars
    para, kgrid, n = config.userdata
    k = K[1]
    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0
    
    τ1, τ2, τ3 = T[1], T[2], T[3]
    # θ = 2.533
    # φ = 0.0
    
    θ = Ang[1]
    φ = Phi[1]
    extK1 = [sin(0.5θ), cos(0.5θ), 0] * kgrid[ExtKidx[1]]
    extK2 = [-sin(0.5θ), cos(0.5θ), 0] * kgrid[ExtKidx[1]]
    extK3 = [sin(0.5θ)*cos(φ),cos(0.5θ),sin(0.5θ)*sin(φ)] * kgrid[ExtKidx[1]]

    k1 = extK3 - extK1 + k
    q = extK1 - k
    qabs = sqrt(dot(q, q))
    # println(q)
    ϵ = dot(k, k) / (2me) - μ
    G1 = ElectronLiquid.Propagator.green_derive(τ2 - τ1, ϵ, β, 0)
    ϵ = dot(k1, k1) / (2me) - μ
    G2d = ElectronLiquid.Propagator.green_derive(τ3 - τ2, ϵ, β, 0)
    G2s = ElectronLiquid.Propagator.green_derive(τ1 - τ2, ϵ, β, 0)
    Vd = ElectronLiquid.interactionDynamic(para, qabs, τ1, τ3)
    Vs = ElectronLiquid.interactionStatic(para, qabs, τ1, τ3)
    l1 = n[1]
    l2 = n[2]
    l3 = n[3]
    factor = 1.0 / (2π)^3
    factor_legend = sin(θ) / (4π)
    f1 = phase((τ1, τ3, τ2, τ2), l1, l2, l3, β) * factor * factor_legend
    f2 = phase((τ1, τ1, τ2, τ2), l1, l2, l3, β) * factor * factor_legend
    return G1 * (G2d * Vd * f1 + G2s * Vs * f2) 
end

#
#   kamp,up   kamp2,down
#         \    /
#            x 
#         /    \
#       /        \
#       |-- < -\   \
#       |  down \    \
#       KO      |- u -|
#       |    up /     | 
#       |-- > -/      |
#       |             |
#   kamp,up       kamp2,down

function _Lver3_exchange(idx, vars, config) 
    K, T, Ang, Phi, ExtKidx = vars
    para, kgrid, n = config.userdata
    k = K[1]
    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0
    
    τ1, τ2, τ3 = T[1], T[2], T[3]
    # θ = 2.533
    # φ = 0.0
    
    θ = Ang[1]
    φ = Phi[1]
    extK1 = [sin(0.5θ), cos(0.5θ), 0] * kgrid[ExtKidx[1]]
    extK2 = [-sin(0.5θ), cos(0.5θ), 0] * kgrid[ExtKidx[1]]
    extK3 = [sin(0.5θ)*cos(φ),cos(0.5θ),sin(0.5θ)*sin(φ)] * kgrid[ExtKidx[1]]

    k1 = extK2 - extK3 + k
    q = extK1 - k
    qabs = sqrt(dot(q, q))
    # println(q)
    ϵ = dot(k, k) / (2me) - μ
    G1 = ElectronLiquid.Propagator.green_derive(τ2 - τ1, ϵ, β, 0)
    ϵ = dot(k1, k1) / (2me) - μ
    G2d = ElectronLiquid.Propagator.green_derive(τ3 - τ2, ϵ, β, 0)
    G2s = ElectronLiquid.Propagator.green_derive(τ1 - τ2, ϵ, β, 0)
    Vd = ElectronLiquid.interactionDynamic(para, qabs, τ1, τ3)
    Vs = ElectronLiquid.interactionStatic(para, qabs, τ1, τ3)
    l1 = n[1]
    l2 = n[2]
    l3 = n[3]
    factor = 1.0 / (2π)^3
    factor_legend = sin(θ) / (4π)
    f1 = phase((τ1, τ2, τ2, τ3), l1, l2, l3, β) * factor * factor_legend
    f2 = phase((τ1, τ2, τ2, τ1), l1, l2, l3, β) * factor * factor_legend
    return G1 * (G2d * Vd * f1 + G2s * Vs * f2)
end

#
# kamp,up             kamp2,down
#     |                    | 
#     |      /- < -\       |
#     |     /   up  \      |
#     |-KO-|         |- u -|
#     |     \   up  /      | 
#     |      \- > -/       |
# kamp,up             kamp2,down

function _Lver3_bubble(idx, vars, config) 
    K, T, Ang, Phi, ExtKidx = vars
    para, kgrid, n = config.userdata
    k = K[1]
    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0
    
    τ1, τ2, τ3 = T[1], T[2], T[3]
    # θ = 2.533
    # φ = 0.0
    
    θ = Ang[1]
    φ = Phi[1]
    extK1 = [sin(0.5θ), cos(0.5θ), 0] * kgrid[ExtKidx[1]]
    extK2 = [-sin(0.5θ), cos(0.5θ), 0] * kgrid[ExtKidx[1]]
    extK3 = [sin(0.5θ)*cos(φ),cos(0.5θ),sin(0.5θ)*sin(φ)] * kgrid[ExtKidx[1]]

    k1 = extK3 - extK1 + k
    q = extK3 - extK1
    qabs = sqrt(dot(q, q))
    # println(q)
    ϵ = dot(k, k) / (2me) - μ
    G1d = ElectronLiquid.Propagator.green_derive(τ2 - τ3, ϵ, β, 0)
    G1s = ElectronLiquid.Propagator.green_derive(τ2 - τ1, ϵ, β, 0)
    ϵ = dot(k1, k1) / (2me) - μ
    G2d = ElectronLiquid.Propagator.green_derive(τ3 - τ2, ϵ, β, 0)
    G2s = ElectronLiquid.Propagator.green_derive(τ1 - τ2, ϵ, β, 0)
    Vd = ElectronLiquid.interactionDynamic(para, qabs, τ3, τ1)
    Vs = ElectronLiquid.interactionStatic(para, qabs, τ3, τ1)
    l1 = n[1]
    l2 = n[2]
    l3 = n[3]
    factor = 1.0 / (2π)^3
    factor_legend = sin(θ) / (4π)
    f1 = phase((τ1, τ1, τ2, τ2), l1, l2, l3, β) * factor * factor_legend
    return  (G1d * G2d * Vd + G1s * G2s * Vs) * f1
end

function _measure_ver3(idx, vars, obs, weight, config) # for the mcmc algorithm
    # K, T, N, Ang, ExtKidx = vars
    # ki = ExtKidx[1]
    ki = vars[5][1]
    obs[idx][ki] += weight
end

function bLambda_MCMC(para::ParaMC;
    kamp=[para.kF],
    n=[0, 0, 0], 
    neval=1e6,
    config=nothing,
    integrand=_PP, # or _Lver3_direct, _Lver3_exchange, _Lver3_bubble
    filename=nothing,
    kwargs...
)
    UEG.MCinitialize!(para, false)
    
    kF = para.kF
    kgrid = kamp
    K = MCIntegration.FermiK(3, kF, 0.2 * kF, 10.0 * kF)
    T = MCIntegration.Continuous(0.0, para.β, offset=1)
    T.data[1] = 0.0
    Ang = MCIntegration.Continuous(0.0, π)
    Phi = MCIntegration.Continuous(0.0, 2π)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid)) 
    
    dof = [[1,2,1,1,1],]
    obs = [zeros(ComplexF64, length(kgrid)),]

    if isnothing(config)
        config = Configuration(;
            var=(K, T, Ang, Phi, ExtKidx),
            dof=dof,
            type=ComplexF64, # type of the integrand
            obs=obs,
            userdata=(para, kgrid, n)
        )
    end

    result = integrate(integrand; measure=_measure_ver3, config=config, solver=:mcmc, neval=neval, kwargs...)
    
    if isnothing(result) == false
        avg, std = result.mean[1], result.stdev[1]
        r = measurement.(real(avg), real(std))
        i = measurement.(imag(avg), imag(std))
        # ver3 = Complex.(r, i)
        ver3 = r
        # println("ver3=$(ver3)")
        # println("filename=$(filename)")
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (kamp, ver3)
            end
        end

        return ver3, result
    else
        return nothing, nothing
    end
end

# function bLambda(para::ParaMC;
#     kamp=[para.kF],
#     n=[0, 0, 0], 
#     neval=1e6,
#     config=nothing,
#     integrand=_PP, # or _Lver3_direct, _Lver3_exchange, _Lver3_bubble
#     filename=nothing,
#     kwargs...
# )
#     UEG.MCinitialize!(para, false)
    
#     kF = para.kF
#     kgrid = kamp
#     K = MCIntegration.FermiK(3, kF, 0.2 * kF, 10.0 * kF)
#     T = MCIntegration.Continuous(0.0, para.β, offset=1)
#     T.data[1] = 0.0
#     Ang = MCIntegration.Continuous(0.0, π)
#     Phi = MCIntegration.Continuous(0.0, 2π)
#     ExtKidx = MCIntegration.Discrete(1, length(kgrid)) 
    
#     dof = [[1,2,1,1,1],]
#     obs = [zeros(ComplexF64, length(kgrid)),]

#     if isnothing(config)
#         config = Configuration(;
#             var=(K, T, Ang, Phi, ExtKidx),
#             dof=dof,
#             type=ComplexF64, # type of the integrand
#             obs=obs,
#             userdata=(para, kgrid, n)
#         )
#     end

#     result = integrate(integrand; measure=_measure_ver3, config=config, solver=:mcmc, neval=neval, kwargs...)
#     if isnothing(result) == false
        
#         avg, std = result.mean[1], result.stdev[1]
#         r = measurement.(real(avg), real(std))
#         i = measurement.(imag(avg), imag(std))
#         # ver3 = Complex.(r, i)
#         ver3 = r
#         return ver3, result
#     else
#         return nothing, nothing
#     end
# end

# function bLambda_MCMC(para::ParaMC;
#     kamp=[para.kF],
#     n=[0, 0, 0], 
#     neval=1e6,
#     config=nothing,
#     integrand=_PP, # or _Lver3_direct, _Lver3_exchange, _Lver3_bubble
#     filename=nothing,
#     kwargs...
# )
#     ver3, result = bLambda(para; kamp, n, neval, config, integrand, filename)
#     if isnothing(ver3) == false
#         for (idx, kΛ) in enumerate(kamp)
#             println("k=$(kΛ/para.kF)kF, ver3=$(ver3)")
#         end

#         if isnothing(filename) == false
#             jldopen(filename, "a+") do f
#                 key = "$(UEG.short(para))"
#                 if haskey(f, key)
#                     @warn("replacing existing data for $key")
#                     delete!(f, key)
#                 end
#                 f[key] = (kamp, ver3)
#             end
#         end

#         return ver3, result
#     end
# end