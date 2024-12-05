function aLambda_avg(para; kamp=[para.kF,], kamp2=kamp, n=[0, 0, 0, 0], theta=[0,], phi=[0,],
    neval=1e6, filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    filter=[NoHartree],    # filter=[NoHartree, NoBubble, Proper],
    channels=[PHr, PHEr, PPr, Alli],
    partition=UEG.partition(para.order),
    transferLoop=nothing, extK=nothing, optimize_level=1,
    verbose=0
)
    # data = Ver4.MC_Spec_Jl(para; kamp=kamp, kamp2=kamp2, theta=theta, phi=phi, n=n, neval=neval, filter=filter, channels=channels, transferLoop=transferLoop)
    ver4, result = Ver4.MC_Spec_Jl(para; kamp=kamp, kamp2=kamp2, theta=theta, phi=phi, n=n, neval=neval, filename=filename,filter=filter, channels=channels, transferLoop=transferLoop)
    return ver4, result
end

# function dcLambdas_avg(para; kΛ, a_s, theta=[0,], neval=1e6)
#     function integrand_c(idx, vars, config)
#         K, N, T, Ang, ExtKidx = vars
#         para ,kgrid, ngrid, thetagrid= config.userdata
        
#         dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0

#         # extK1 = [1.0, 0.0, 0.0] *kgrid[ExtKidx[1]]
#         # extK2 = [0.0, 0.0, 0.0]
#         θ = thetagrid[Ang[1]] 
#         extK1 = [sin(0.5θ),cos(0.5θ),0] * kgrid[ExtKidx[1]]
#         extK2 = [-sin(0.5θ),cos(0.5θ), 0] * kgrid[ExtKidx[1]]
#         k = K[1]
#         k1 = extK1 + extK2 -k
#         # k1 = extK1 + k
#         # k1 = extK1 - k
#         τ1 = 0.0
#         τ2 = T[1]
#         ϵ1 = dot(k,k)/(2.0*me) - μ
#         ϵ2 = dot(k1, k1)/(2.0*me) - μ        
#         G1 = ElectronLiquid.Propagator.green_derive(τ2-τ1, ϵ1, β, 0)
#         G2 = ElectronLiquid.Propagator.green_derive(τ2-τ1, ϵ2, β, 1)

#         factor_d = dot([0.0,2.0*cos(0.5θ),0.0],k1)/me
        
#         # println(factor_d)
#         # println(extK1+extK2)

#         factor = 1.0/(2π)^3

#         l1 = ngrid[N[1]][1]
#         l2 = ngrid[N[1]][2]
#         l3 = ngrid[N[1]][3]
#         f1 = phase(l1,l2,l3,τ1,τ1,τ2,τ2,β)
#         return G1 * G2 *factor *factor_d *f1  
#         # return G1 *(G21*Vd*f1 + G22*Vs*f2)* uΛ *para.NF
#     end

#     function measure(idx, vars, obs, weight, config) # for the mcmc algorithm
#         n_theta = vars[4][1]
#         obs[idx][n_theta] += weight
#     end

#     kgrid = [kΛ]
#     ngrid = [[0, 0, 0]]
#     thetagrid = theta
#     # K = MCIntegration.FermiK(3, kΛ, 0.02 , 30.0)
#     K = MCIntegration.FermiK(3, kΛ, 0.2 * kΛ, 10.0 * kΛ)
#     ExtKidx = MCIntegration.Discrete(1, length(kgrid)) 
#     N = MCIntegration.Discrete(1, length(ngrid))
#     Ang = MCIntegration.Discrete(1, length(thetagrid))
#     T = MCIntegration.Continuous(0.0, para.β)
#     dof = [[1,1,1,1,1],]
#     obs = [zeros(ComplexF64, length(thetagrid))]
#     config = Configuration(;
#             var=(K, N, T, Ang, ExtKidx),
#             dof=dof,
#             type=ComplexF64, # type of the integrand
#             obs=obs,
#             userdata=(para ,kgrid, ngrid, thetagrid)
#     )
#     neval = neval
#     result = integrate(integrand_c; measure=measure, config=config, solver=:mcmc, neval=neval)
#     # result = integrate(integrand_test; measure=measure, config=config, solver=:mcmc, neval=neval)
#     avg, std = result.mean[1], result.stdev[1]
#     r = measurement.(real(avg), real(std))
#     # println(r)
#     data = [0.0, 0.0]
#     Nθ = length(thetagrid)
#     avga_theta = zeros(Float64, Nθ)
#     avge_theta = zeros(Float64, Nθ)
#     for (ti, _theta) in enumerate(thetagrid)
#         avga_theta[ti] = r[ti].val * sin(_theta)
#         avge_theta[ti] = r[ti].err * sin(_theta)
#     end
#     data[1] = Interp.integrate1D(avga_theta, thetagrid)
#     data[2] = Interp.integrate1D(avge_theta, thetagrid)
#     data = data/(2.0*para.NF) *para.kF
    
#     return measurement.(data[1], data[2])
# end

function c_coeff_pp(para, kamp=para.kF, kamp2=para.kF)
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 32)
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * cos(θ) * kamp * kamp2) for θ in θgrid.grid]

    Wp = zeros(ComplexF64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = ElectronGas.Polarization.Ladder0_FiniteTemp(q, 0, para)
    end

    vud = Interp.integrate1D(Wp .* sin.(θgrid.grid), θgrid) / 2
    return 0, -real(vud)
end


function cLambda_Pi(para; kΛ)
    
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 32, 0.001, 32)
    qs = [kΛ * sqrt(2.0*(1 - cos(θ))) for θ in θgrid.grid]
    Wp = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = UEG.polarKW(q, 0, para)
    end
    # println(Wp)
    avgPiu = Ver4.Legrendre(0, Wp, θgrid) 
    # println("Piu=$(avgPiu)")

    # θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 32, 0.001, 32)
    # φgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, 2π], [0.0, 2π], 32, 0.001, 32)
    

    # avg_theta = zeros(Float64, length(θgrid))
    # avg_phi = zeros(Float64, length(φgrid))
    # for (idx, _theta) in enumerate(θgrid)
    #     for (jdx, _phi) in enumerate(φgrid)
    #         extK1 = [sin(0.5_theta), cos(0.5_theta), 0] * kΛ
    #         extK2 = [-sin(0.5_theta), cos(0.5_theta), 0] * kΛ
    #         extK3 = [sin(0.5_theta)*cos(_phi),cos(0.5_theta),sin(0.5_theta)*sin(_phi)] * kΛ
    #         # qabs = sqrt(dot(extK1-extK2, extK1-extK2))
    #         qabs = sqrt(dot(extK2-extK3, extK2-extK3))
    #         avg_phi[jdx] = UEG.polarKW(qabs, 0, para)
            
    #     end
    #     avg_theta[idx] = Interp.integrate1D(avg_phi, φgrid) / (2π)
    # end
    # avgPiu = Ver4.Legrendre(0, avg_theta, θgrid)
    # println("Piu=$(avgPiu)")
   
    avgPiu = avgPiu/para.NF
    return avgPiu
end

function aLambda_PH_avg(para::ParaMC;
    kamp=[para.kF],
    n=[0, 0, 0], 
    neval=1e6,
    config=nothing,
    filename=nothing,
    kwargs...
    )

    function integrand(idx, vars, config) 
        K, T, Ang, Phi, ExtKidx = vars
        para, kgrid, n = config.userdata
        k = K[1]
        dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0
        
        τ1, τ2, τ3, τ4 = T[1], T[2], T[3], T[4]
        # θ = 2.533
        # φ = 0.0
        
        θ = Ang[1]
        φ = Phi[1]
        extK1 = [sin(0.5θ), cos(0.5θ), 0] * kgrid[ExtKidx[1]]
        extK2 = [-sin(0.5θ), cos(0.5θ), 0] * kgrid[ExtKidx[1]]
        extK3 = [sin(0.5θ)*cos(φ),cos(0.5θ),sin(0.5θ)*sin(φ)] * kgrid[ExtKidx[1]]
    
        k1 = extK3 - extK1 + k
        q1 = extK1 - k
        q2 = extK3 - extK1
        
        # println(q)
        ϵ = dot(k, k) / (2me) - μ
        G1d = ElectronLiquid.Propagator.green_derive(τ4 - τ1, ϵ, β, 0)
        G1s = ElectronLiquid.Propagator.green_derive(τ2 - τ1, ϵ, β, 0)
        ϵ = dot(k1, k1) / (2me) - μ
        G2dd = ElectronLiquid.Propagator.green_derive(τ3 - τ4, ϵ, β, 0)
        G2sd = ElectronLiquid.Propagator.green_derive(τ1 - τ4, ϵ, β, 0)
        G2ds = ElectronLiquid.Propagator.green_derive(τ3 - τ2, ϵ, β, 0)
        G2ss = ElectronLiquid.Propagator.green_derive(τ1 - τ2, ϵ, β, 0)
        qabs = sqrt(dot(q1, q1))
        V1d = ElectronLiquid.interactionDynamic(para, qabs, τ3, τ1)
        V1s = ElectronLiquid.interactionStatic(para, qabs, τ3, τ1)
        qabs = sqrt(dot(q2, q2))
        V2d = ElectronLiquid.interactionDynamic(para, qabs, τ2, τ4)
        V2s = ElectronLiquid.interactionStatic(para, qabs, τ2, τ4)
        l1 = n[1]
        l2 = n[2]
        l3 = n[3]
        factor = 1.0 / (2π)^3
        factor_legend = sin(θ) / (4π)
        f1 = phase((τ1, τ3, τ2, τ2), l1, l2, l3, β) * factor * factor_legend
        f2 = phase((τ1, τ1, τ2, τ2), l1, l2, l3, β) * factor * factor_legend
        return  (V1d * (G1d * G2dd * V2d + G1s * G2ds * V2s) * f1 + V1s * (G1d * G2sd * V2d + G1s * G2ss * V2s) * f2) * para.NF
    end

    function measure(idx, vars, obs, weight, config) # for the mcmc algorithm
        # K, T, N, Ang, ExtKidx = vars
        # ki = ExtKidx[1]
        ki = vars[5][1]
        obs[idx][ki] += weight
    end

    UEG.MCinitialize!(para, false)
    
    kF = para.kF
    kgrid = kamp
    K = MCIntegration.FermiK(3, kF, 0.2 * kF, 10.0 * kF)
    T = MCIntegration.Continuous(0.0, para.β, offset=1)
    T.data[1] = 0.0
    Ang = MCIntegration.Continuous(0.0, π)
    Phi = MCIntegration.Continuous(0.0, 2π)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid)) 
    
    dof = [[1,3,1,1,1],]
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

    result = integrate(integrand; measure=measure, config=config, solver=:mcmc, neval=neval, kwargs...)
    if isnothing(result) == false
        avg, std = result.mean[1], result.stdev[1]
        r = measurement.(real(avg), real(std))
        i = measurement.(imag(avg), imag(std))
        # ver3 = Complex.(r, i)
        ver4 = r

        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (kamp, ver4)
            end
        end

        return ver4, result
    end
end

function aLambda_thetaphi_avg(para; kamp=[para.kF,], kamp2=kamp, n=[0, 0, 0, 0], theta=[0,], phi=[0,],
    neval=1e6, filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    filter=[NoHartree],    # filter=[NoHartree, NoBubble, Proper],
    channels=[PHr, PHEr, PPr, Alli],
    partition=UEG.partition(para.order),
    transferLoop=nothing, extK=nothing, optimize_level=1,
    verbose=0
)
    # data = Ver4.MC_Spec_Jl(para; kamp=kamp, kamp2=kamp2, theta=theta, phi=phi, n=n, neval=neval, filter=filter, channels=channels, transferLoop=transferLoop)
    ver4, result = Ver4.MC_Spec_Jl(para; kamp=kamp, kamp2=kamp2, theta=theta, phi=phi, n=n, neval=neval,filter=filter, channels=channels, transferLoop=transferLoop, verbose=-1)
    partition = [(1, 0, 0)]
    if isnothing(ver4) == false && verbose != -1
        # for (p, data) in ver4
        datadict = Dict{eltype(partition),Any}()
        thetagrid = CompositeGrids.SimpleG.Arbitrary(theta)
        phigrid = CompositeGrids.SimpleG.Arbitrary(phi)
        for p in partition
            data = ver4[p]
            data_aΛ = []
            Nθ = length(theta)
            Nφ = length(phi)
            for (ki, k) in enumerate(kamp)
                avg = [0.0, 0.0]
                avga_theta = zeros(Float64, Nθ)
                avge_theta = zeros(Float64, Nθ)
                for (ti, _theta) in enumerate(theta)
                    if Nφ == 1
                        d2 = real(data[2, ti, 1, ki])
                        avga_theta[ti] = d2.val * sin(_theta)
                        avge_theta[ti] = d2.err * sin(_theta)
                    else
                        da_phi = zeros(Float64, Nφ)
                        de_phi = zeros(Float64, Nφ)
                        for (pidx, _phi) in enumerate(phi)
                            d2 = real(data[2, ti, pidx, ki])
                            da_phi[pidx] = d2.val
                            de_phi[pidx] = d2.err
                        end
                        avga_theta[ti] = Interp.integrate1D(da_phi, phigrid) * sin(_theta) / π
                        avge_theta[ti] = Interp.integrate1D(de_phi, phigrid) * sin(_theta) / π
                    end
                end
                avg[1] = Interp.integrate1D(avga_theta, thetagrid) / 2.0
                avg[2] = Interp.integrate1D(avge_theta, thetagrid) / 2.0
                push!(data_aΛ, measurement.(avg[1], avg[2]))
            end
            datadict[p] = data_aΛ
        end

        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (kamp, n, datadict)
            end
        end
    end
    return ver4, result
end