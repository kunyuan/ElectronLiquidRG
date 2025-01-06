using ElectronLiquidRG
using ElectronLiquid
using CompositeGrids
using LinearAlgebra, Measurements
using JLD2
using GreenFunc
using DelimitedFiles
using DataFrames
using CSV

neval = 1e7
dim = 3
rs = [2.0]
mass2 = [1e-2,]
beta = [40.0]
order = [1,]
isDynamic = true

function get_data_z(para, filename; Fs=[0.0], Λgrid=[para.kF])
    f = jldopen(filename, "r")
    sw = Dict()
    mu = Dict()
    partition = UEG.partition(para.order)
    for p in partition
        sw[p] = MeshArray(Fs, Λgrid; dtype=Measurement{Float64})
        mu[p] = 0.0 # we don't need mu for now
    end
    for (fi, F) in enumerate(Fs)
        for key in keys(f)
            para_data = ParaMC(key)
            para_data_Fs0 = ElectronLiquidRG.get_para(para_data, -0.0)
            if para_data_Fs0 == para
                if isapprox(para_data.Fs, F) == true
                    # println("Fs=$Fs_idx, key=$(key)")
                    _para = ElectronLiquidRG.get_para(para, F)
                    ngrid, kgrid, sigma = f[key]
                    # println(length(kgrid), ", ", length(Λgrid))
                    @assert kgrid ≈ Λgrid "length(kgrid) = $(length(kgrid)), length(Λgrid) = $(length(Λgrid))"
                    for p in partition
                        for (ki, kΛ) in enumerate(Λgrid)
                            sw[p][fi, ki] = ElectronLiquidRG.zfactor(sigma[p][:, ki], _para.β)
                            # println("Fs=$F, k=$(kΛ/para.kF)kF, sigma=$(sigma[p][:, ki]), dzi=$(sw[p][fi, ki])")
                        end
                    end
                end
            end

        end
    end

    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)
    return dzi, dmu, dz
end

for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)
    _para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
    kF = _para.kF
    # Λgrid = [1.0] *_para.kF
    # Λgrid = range(1.0*_para.kF, 3.0*_para.kF,21)
    Λgrid = ElectronLiquidRG.CompositeGrid.LogDensedGrid(:cheb, [1.0 * _para.kF, 8.0 * _para.kF], [_para.kF,], 6, 0.02 * _para.kF, 8)
    cutoff_Λ = 29
    label_Λ = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,28,29,30,31,32,33,34,35,36,37,38,39,40] 
    # label_Λ = cat(label_Λ, [cutoff_Λ:length(Λgrid)]; dims=1)
    println(label_Λ)
    Λgrid_measurement = zeros(Float64, length(label_Λ))
    for (idx, idx_Λ) in enumerate(label_Λ)
        Λgrid_measurement[idx] = Λgrid[idx_Λ]
    end
    println(Λgrid_measurement/kF)
    # # rs = 1.0
    # Fs = SimpleG.Uniform([-0.4,0.0],41).grid 
    # Fs1 = SimpleG.Uniform([-0.99,-0.41], 59).grid 
    # Fs2 = SimpleG.Uniform([-1.5,-1.01], 50).grid 
    # Fs1 = cat(Fs2,Fs1; dims=1)
    # # rs = 5.0
    # Fs = SimpleG.Uniform([-0.99,0.0],100).grid
    # Fs1 = SimpleG.Uniform([-2.4,-1.01], 140).grid
    # # rs = 3.0
    # Fs1 = SimpleG.Uniform([-0.99, -0.0], 100).grid 
    # Fs2 = SimpleG.Uniform([-1.6,-1.01], 60).grid 
    # Fs1 = cat(Fs2,Fs1; dims=1)
    # rs = 2.0
    Fs1 = SimpleG.Uniform([-0.99, -0.0], 100).grid 
    Fs2 = SimpleG.Uniform([-1.5,-1.01], 50).grid 
    Fs1 = cat(Fs2,Fs1; dims=1)

    filename_sigma  = "data/sigma_RG.jld2"
    # filename_sigma  = "data_test/sigma_RG_1.jld2"
    filename_data_zfactor_RG = "data_csv/zfactor_RG_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"
    dsigma_RG = DataFrame()
    dzi, dmu, dz = get_data_z(_para, filename_sigma; Fs=Fs1, Λgrid=Λgrid_measurement)
    # println(dz)
    for (jdx, Fs_jdx) in enumerate(Fs1)
        println("Fs=$(Fs_jdx)")
        colname = "Fs_$(Fs_jdx)"
        dsigma_RG[!, colname] = dz[1][jdx,:]
    end

    # dzi, dmu, dz = get_data_z(_para, filename_sigma; Fs=Fs, Λgrid=Λgrid)
    # for (jdx, Fs_jdx) in enumerate(Fs)
    #     println("Fs=$(Fs_jdx)")
    #     colname = "Fs_$(Fs_jdx)"
    #     data = []
    #     for (idx, idx_Λ) in enumerate(label_Λ)
    #         push!(data, dz[1][jdx, idx_Λ])
    #     end
    #     # if isapprox(Fs_jdx, -0.99)
    #     #     println(data)
    #     # end
    #     dsigma_RG[!, colname] = data
    # end
    
    CSV.write(filename_data_zfactor_RG, dsigma_RG)
end
