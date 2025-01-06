using ElectronLiquidRG
using ElectronLiquid
using CompositeGrids
using LinearAlgebra, Measurements
using JLD2
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


function get_data_aLambda(data, kamp, theta, phi)
    Nθ = length(theta)
    Nφ = length(phi)
    aΛ = [] 
    for (ki, kΛ) in enumerate(kamp)
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
                avga_theta[ti] = Interp.integrate1D(da_phi, phi) * sin(_theta) / π
                avge_theta[ti] = Interp.integrate1D(de_phi, phi) * sin(_theta) / π
            end
        end
        avg[1] = Interp.integrate1D(avga_theta, theta) / 2.0
        avg[2] = Interp.integrate1D(avge_theta, theta) / 2.0
        push!(aΛ, measurement.(avg[1], avg[2]))
    end
    return aΛ
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
    println(length(label_Λ))
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
    
    
    
    
    filename_data_a_PP_PHE = "data_csv/aLambda_PP_PHE_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"
    filename_data_a_PH = "data_csv/aLambda_PH_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"
    da_PP_PHE = DataFrame()
    da_PH = DataFrame()
    # f = jldopen(filename_aPP, "r")
    # println(f)

    # new version
    filename_aPP  = "data/aLambda_PP_PHE.jld2"
    f = jldopen(filename_aPP, "r")
    # println(keys(f))
    for (idx, F) in enumerate(Fs1)
        for key in keys(f)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    println(key)
                    # println("Fs=$(F), key_cluster=$(key)")
                    kamp, n, datadict = f[key]
                    # println(datadict)
                    data = datadict[(1, 0, 0)]
                    # data = data[cutoff_Λ:end]
                    colname = "Fs_$(F)"
                    da_PP_PHE[!, colname] = data
                    println(length(data))
                end
            end
        end
    end

    # # println(keys(f))
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 println(key)
    #                 # println("Fs=$(F), key_cluster=$(key)")
    #                 kamp, n, datadict = f[key]
    #                 # println(datadict)
    #                 data = datadict[(1, 0, 0)]
    #                 data_import = []
    #                 for (idx, idx_Λ) in enumerate(label_Λ)
    #                     push!(data_import, data[idx_Λ])
    #                 end
    #                 # data = data[cutoff_Λ:end]
    #                 colname = "Fs_$(F)"
    #                 da_PP_PHE[!, colname] = data_import
    #             end
    #         end
    #     end
    # end
    CSV.write(filename_data_a_PP_PHE, da_PP_PHE)

    filename_aPH  = "data/aLambda_PH.jld2"
    f = jldopen(filename_aPH, "r")
    for (idx, F) in enumerate(Fs1)
        for key in keys(f)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    # println("Fs=$(F), key_cluster=$(key)")
                    kamp, n, datadict = f[key]
                    # println(datadict)
                    data = datadict[(1, 0, 0)]
                    # data = data[cutoff_Λ:end]
                    colname = "Fs_$(F)"
                    da_PH[!, colname] = data
                end
            end
        end
    end

    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 # println("Fs=$(F), key_cluster=$(key)")
    #                 kamp, n, datadict = f[key]
    #                 # println(datadict)
    #                 data = datadict[(1, 0, 0)]
    #                 data_import = []
    #                 for (idx, idx_Λ) in enumerate(label_Λ)
    #                     push!(data_import, data[idx_Λ])
    #                 end
    #                 # data = data[cutoff_Λ:end]
    #                 colname = "Fs_$(F)"
    #                 da_PH[!, colname] = data_import
    #             end
    #         end
    #     end
    # end
    CSV.write(filename_data_a_PH, da_PH)


    # old version
    # f = jldopen(filename_aPP, "r")
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 # println("Fs=$(F), key_cluster=$(key)")
    #                 kamp, n, theta, phi, ver4 = f[key]
    #                 thetagrid = CompositeGrids.SimpleG.Arbitrary(theta)
    #                 phigrid = CompositeGrids.SimpleG.Arbitrary(phi)
    #                 # println(thetagrid)
    #                 data = ver4[(1, 0, 0)]
    #                 a_PP = get_data_aLambda(data, kamp, thetagrid, phigrid)
    #                 # println(a_PP)
    #                 colname = "Fs_$(F)"
    #                 da_PP_PHE[!, colname] = a_PP
    #             end
    #         end
    #     end
    # end

    # f = jldopen(filename_aPH, "r")
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 # println("Fs=$(F), key_cluster=$(key)")
    #                 kamp, n, theta, phi, ver4 = f[key]
    #                 thetagrid = CompositeGrids.SimpleG.Arbitrary(theta)
    #                 phigrid = CompositeGrids.SimpleG.Arbitrary(phi)
    #                 data = ver4[(1, 0, 0)]
    #                 a_PH = get_data_aLambda(data, kamp, thetagrid, phigrid)
    #                 colname = "Fs_$(F)"
    #                 da_PH[!, colname] = a_PH
    #             end
    #         end
    #     end
    # end

    
    
    # test = DataFrame(CSV.File(filename_data_a_PP_PHE))
    # println(test)
end

