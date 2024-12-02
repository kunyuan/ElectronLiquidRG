function sigma(para; neval=1e6, Λgrid=Λgrid(para.kF), filename=nothing)
    sigma, result = Sigma.MC(para; kgrid=Λgrid, ngrid=[-1, 0, 1], neval=neval, filename=filename)
    return sigma, result
end

function zfactor(data, β, ngrid=[-1, 0]) #assume the data are calculated with [-1, 0, 1]
    if ngrid == [0, 1]
        return (imag(data[3]) - imag(data[2])) / (2π / β)
    elseif ngrid == [-1, 0]
        return (imag(data[2]) - imag(data[1])) / (2π / β)
    else
        error("ngrid = $ngrid not implemented")
    end

    # if ngrid == [0, 1]
    #     return @. (imag(data[3, :]) - imag(data[2, :])) / (2π / β)
    # elseif ngrid == [-1, 0]
    #     return @. (imag(data[2, :]) - imag(data[1, :])) / (2π / β)
    # else
    #     error("ngrid = $ngrid not implemented")
    # end
end

function mu(data)
    return real(data[1, 1])
end

function zCT(para, filename; Fs=fdict[para.rs])
    # println("read Fs = $Fs from $filename")
    f = jldopen(filename, "r")
    # z1 = zeros(Measurement{Float64}, length(Fs), length(Λgrid))
    sw = Dict()
    mu = Dict()
    partition = UEG.partition(para.order)
    
    for p in partition
        sw[p] = MeshArray(Fs; dtype=Measurement{Float64})
        mu[p] = 0.0 # we don't need mu for now
    end
    for (fi, F) in enumerate(Fs)
        println(F)
        _para = get_para(para, F)
        key = UEG.short(_para)
        ngrid, kgrid, sigma = f[key]
        
        for p in partition
            # println(sigma[p][1],",",sigma[p][2],",",sigma[p][3])
            sw[p][fi] = zfactor(sigma[p], _para.β)
        end
    end

    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)
    return dzi, dmu, dz
end

# function get_para(para, Fs; order=para.order)
#     UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=Fs, Fa=-0.0, order=order,
#     mass2=para.mass2, isDynamic=true, isFock=false)
# end

function sigma_oneloop(para;kamp=[para.kF],
    neval=1e6,
    config=nothing,
    filename=nothing,
    kwargs...)

    function phase_sigma(extT, nin, nout, β)
        tIn, tOut = extT
        win, wout = π * (2nin + 1) / β, π * (2nout + 1) / β
        
        return exp(-1im * (tIn * win - tOut * wout))
    end
    function integrand(idx, vars, config) # need to benchmark
        K, T, N, ExtKidx = vars
        para, kgrid, ngrid = config.userdata
        k = K[1]
        dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0
        
        τ1, τ2 = T[1], T[2]
        
        extK = [1.0, 0.0, 0.0] * kgrid[ExtKidx[1]]
        
        q = extK - k
        qabs = sqrt(dot(q, q))

        ϵ = dot(k, k) / (2me) - μ
        G1d = ElectronLiquid.Propagator.green_derive(τ2 - τ1, ϵ, β, 0)
        G1s = ElectronLiquid.Propagator.green_derive(τ1 - τ1, ϵ, β, 0)
        
        Vd = ElectronLiquid.interactionDynamic(para, qabs, τ2, τ1)
        Vs = ElectronLiquid.interactionStatic(para, qabs, τ2, τ1)
        l1 = ngrid[N[1]]
        l2 = ngrid[N[1]]

        factor = 1.0 / (2π)^3
        f1 = phase_sigma((τ1, τ2), l1, l2, β) * factor 
        f2 = phase_sigma((τ1, τ1), l1, l2, β) * factor 
        return  (G1d * Vd * f1 + G1s * Vs * f2)
    end
    
    function measure(idx, vars, obs, weight, config) # for the mcmc algorithm
        K, T, N, ExtKidx = vars
        ki = ExtKidx[1]
        ni = N[1]
        obs[idx][ni, ki] += weight
    end

    UEG.MCinitialize!(para, false)
    partition = UEG.partition(para.order)
    
    kF = para.kF
    kgrid = kamp
    K = MCIntegration.FermiK(3, kF, 0.2 * kF, 10.0 * kF)
    T = MCIntegration.Continuous(0.0, para.β, offset=1)
    T.data[1] = 0.0
    ngrid=[-1, 0, 1]    

    ExtKidx = MCIntegration.Discrete(1, length(kgrid)) 
    N = MCIntegration.Discrete(1, length(ngrid))

    dof = [[1,1,1,1],]
    obs = [zeros(ComplexF64, length(ngrid), length(kgrid)),]

    if isnothing(config)
        config = Configuration(;
            var=(K, T, N, ExtKidx),
            dof=dof,
            type=ComplexF64, # type of the integrand
            obs=obs,
            userdata=(para, kgrid, ngrid)
        )
    end

    result = integrate(integrand; measure=measure, config=config, solver=:mcmc, neval=neval, kwargs...)

    if isnothing(result) == false
        datadict = Dict{eltype(partition),Any}()
        for (o, key) in enumerate(partition)
            avg, std = result.mean[o], result.stdev[o]
            r = measurement.(real(avg), real(std))
            i = measurement.(imag(avg), imag(std))
            data = Complex.(r, i)
            datadict[key] = data
        end

        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (ngrid, kgrid, datadict)
            end
        end

        return datadict, result
    end
end