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
rs = [5.0]
mass2 = [1e-2,]
beta = [40.0]
order = [1,]
isDynamic = true

for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)
    _para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
    kF = _para.kF
    # Λgrid = [1.0] *_para.kF
    # Λgrid = range(1.0*_para.kF, 3.0*_para.kF,21)
    Λgrid = ElectronLiquidRG.CompositeGrid.LogDensedGrid(:cheb, [1.0 * _para.kF, 8.0 * _para.kF], [_para.kF,], 6, 0.02 * _para.kF, 8)
    
    # Λgrid = ElectronLiquidRG.CompositeGrid.LogDensedGrid(:gauss, [1.0 * _para.kF, 5 * _para.kF], [_para.kF,], 6, 0.01 * _para.kF, 6)
    # println(size(Λgrid))
    # asgrid = zeros(Float64, length(Λgrid))
    # Fs = ElectronLiquidRG.treelevel_RG(_para, Λgrid, asgrid)
    # println(Fs)
    # Fs = SimpleG.Uniform([-0.99,0.0],100).grid 
    Fs = SimpleG.Uniform([-2.4, -1.81], 60).grid
    # Fs = Fs[1:3] 
    partition = UEG.partition(_para.order)
    # partition = [(1, 0, 0)]
    # println(partition)


    # # import aLambda_Large
    # filename_aPP  = "data_test/aLambda_PP_PHE_Large.jld2"

    # filename_aPP_import  = "data_csv/aLambda_PP_PHE_Large.jld2"
    
    # datadict_import = Dict{eltype(partition),Any}()
    # f = jldopen(filename_aPP, "r")
    # f_ver2 = jldopen(filename_aPP_import, "a+")
    # for (idx, F) in enumerate([-0.0])
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 # println("Fs=$(F), key_cluster=$(key)")
    #                 println(key)
    #                 kamp, n, datadict = f[key]
    #                 # println(datadict)
    #                 # data = datadict[(1, 0, 0)]
    #                 # println(length(kamp), length(datadict[(1,0,0)]))
    #                 # datadict[(1, 0, 0)] = data
    #                 if haskey(f_ver2, key)
    #                     @warn("replacing existing data for $key")
    #                     delete!(f_ver2, key)
    #                 end
    #                 f_ver2[key] = (kamp, n, datadict)
    #             end
    #         end
    #     end
    # end

    # import aLambda 
    filename_aPP  = "data_test/aLambda_PP_PHE.jld2"
    filename_aPH  = "data_test/aLambda_PH.jld2"

    filename_aPP_import  = "data/aLambda_PP_PHE.jld2"
    filename_aPH_import  = "data/aLambda_PH.jld2"
    
    datadict_import = Dict{eltype(partition),Any}()
    f = jldopen(filename_aPP, "r")
    f_ver2 = jldopen(filename_aPP_import, "a+")
    for (idx, F) in enumerate(Fs)
        for key in keys(f)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    # println("Fs=$(F), key_cluster=$(key)")
                    println(key)
                    kamp, n, datadict = f[key]
                    # println(datadict)
                    # data = datadict[(1, 0, 0)]
                    # println(length(kamp), length(datadict[(1,0,0)]))
                    # datadict[(1, 0, 0)] = data
                    if haskey(f_ver2, key)
                        @warn("replacing existing data for $key")
                        delete!(f_ver2, key)
                    end
                    f_ver2[key] = (kamp, n, datadict)
                end
            end
        end
    end

    f = jldopen(filename_aPH, "r")
    f_ver2 = jldopen(filename_aPH_import, "a+")
    for (idx, F) in enumerate(Fs)
        for key in keys(f)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    # println("Fs=$(F), key_cluster=$(key)")
                    kamp, n, datadict = f[key]
                    # println(datadict)
                    # data = datadict[(1, 0, 0)]
                    # println(length(kamp), length(datadict[(1,0,0)]))

                    # datadict[(1, 0, 0)] = data
                    if haskey(f_ver2, key)
                        @warn("replacing existing data for $key")
                        delete!(f_ver2, key)
                    end
                    f_ver2[key] = (kamp, n, datadict)
                end
            end
        end
    end

    # import bLambda
    filename = "data_test/bLambda_PP.jld2"
    filename_import = "data/bLambda_PP.jld2"
    f = jldopen(filename, "r")
    f_ver2 = jldopen(filename_import, "a+")
    for (idx, F) in enumerate(Fs)
        for key in keys(f)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    kamp, ver3 = f[key]
                    # println(length(kamp))
                    if haskey(f_ver2, key)
                        @warn("replacing existing data for $key")
                        delete!(f_ver2, key)
                    end
                    f_ver2[key] = (kamp, ver3)
                end
            end
        end
    end

    # _Lver3_direct channel
    filename = "data_test/bLambda_Lver3_direct.jld2"
    filename_import = "data/bLambda_Lver3_direct.jld2"
    f = jldopen(filename, "r")
    f_ver2 = jldopen(filename_import, "a+")
    for (idx, F) in enumerate(Fs)
        for key in keys(f)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    kamp, ver3 = f[key]
                    if haskey(f_ver2, key)
                        @warn("replacing existing data for $key")
                        delete!(f_ver2, key)
                    end
                    f_ver2[key] = (kamp, ver3)
                end
            end
        end
    end

    # _Lver3_exchange channel
    filename = "data_test/bLambda_Lver3_exchange.jld2"
    filename_import = "data/bLambda_Lver3_exchange.jld2"
    f = jldopen(filename, "r")
    f_ver2 = jldopen(filename_import, "a+")
    for (idx, F) in enumerate(Fs)
        for key in keys(f)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    kamp, ver3 = f[key]
                    if haskey(f_ver2, key)
                        @warn("replacing existing data for $key")
                        delete!(f_ver2, key)
                    end
                    f_ver2[key] = (kamp, ver3)
                end
            end
        end
    end

    # _Lver3_bubble channel
    filename = "data_test/bLambda_Lver3_bubble.jld2"
    filename_import = "data/bLambda_Lver3_bubble.jld2"
    f = jldopen(filename, "r")
    f_ver2 = jldopen(filename_import, "a+")
    for (idx, F) in enumerate(Fs)
        for key in keys(f)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    kamp, ver3 = f[key]
                    if haskey(f_ver2, key)
                        @warn("replacing existing data for $key")
                        delete!(f_ver2, key)
                    end
                    f_ver2[key] = (kamp, ver3)
                end
            end
        end
    end

    # import sigma
    filename_sigma  = "data_test/sigma_RG.jld2"
    filename_sigma_import  = "data/sigma_RG.jld2"
    f = jldopen(filename_sigma, "r")
    f_ver2 = jldopen(filename_sigma_import, "a+")
    for (fi, F) in enumerate(Fs)
        for key in keys(f)
            para_data = ParaMC(key)
            para_data_Fs0 = ElectronLiquidRG.get_para(para_data, -0.0)
            if para_data_Fs0 == _para
                if isapprox(para_data.Fs, F) == true
                    # println("Fs=$Fs_idx, key=$(key)")
                    # _para = ElectronLiquidRG.get_para(para, Fs_idx)
                    ngrid, kgrid, sigma = f[key]
                    println(key)
                    if haskey(f_ver2, key)
                        @warn("replacing existing data for $key")
                        delete!(f_ver2, key)
                    end
                    f_ver2[key] = (ngrid, kgrid, sigma)
                end
            end
        end
    end

    # # import aLambda 
    # filename_aPP  = "data/aLambda_PP_PHE.jld2"
    # filename_aPH  = "data/aLambda_PH.jld2"

    # filename_aPP_1  = "data_rs=5_F-0.3-0.0_Large/aLambda_PP_PHE.jld2"
    # filename_aPH_1  = "data_rs=5_F-0.3-0.0_Large/aLambda_PH.jld2"

    # filename_aPP_import  = "data_test/aLambda_PP_PHE_1.jld2"
    # filename_aPH_import  = "data_test/aLambda_PH_1.jld2"

    # # partition = [(1, 0, 0)]
    # partition = UEG.partition(_para.order)
    # # println(partition)
    # datadict_import = Dict{eltype(partition),Any}()
    # f = jldopen(filename_aPP, "r")
    # f_1 = jldopen(filename_aPP_1, "r")
    # f_ver2 = jldopen(filename_aPP_import, "a+")
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 # println("Fs=$(F), key_cluster=$(key)")
    #                 println(key)
    #                 kamp, n, datadict = f[key]
    #                 kamp_1, n_1, datadict_1 = f_1[key]
    #                 # println(datadict)
    #                 # println(cat(kamp, kamp_1; dims=1)/kF)
    #                 kamp_import = cat(kamp, kamp_1; dims=1)
    #                 datadict_import[(1,0,0)] = cat(datadict[(1,0,0)], datadict_1[(1,0,0)]; dims=1)
    #                 # println(datadict_import)
                    
    #                 if haskey(f_ver2, key)
    #                     @warn("replacing existing data for $key")
    #                     delete!(f_ver2, key)
    #                 end
    #                 f_ver2[key] = (kamp_import, n, datadict_import)
    #             end
    #         end
    #     end
    # end

    # datadict_import = Dict{eltype(partition),Any}()
    # f = jldopen(filename_aPH, "r")
    # f_1 = jldopen(filename_aPH_1, "r")
    # f_ver2 = jldopen(filename_aPH_import, "a+")
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 # println("Fs=$(F), key_cluster=$(key)")
    #                 println(key)
    #                 kamp, n, datadict = f[key]
    #                 kamp_1, n_1, datadict_1 = f_1[key]
    #                 # println(datadict)
    #                 # println(cat(kamp, kamp_1; dims=1)/kF)
    #                 kamp_import = cat(kamp, kamp_1; dims=1)
    #                 datadict_import[(1,0,0)] = cat(datadict[(1,0,0)], datadict_1[(1,0,0)]; dims=1)
    #                 # println(datadict_import)
                    
    #                 if haskey(f_ver2, key)
    #                     @warn("replacing existing data for $key")
    #                     delete!(f_ver2, key)
    #                 end
    #                 f_ver2[key] = (kamp_import, n, datadict_import)
    #             end
    #         end
    #     end
    # end

    # import bLambda
    # filename = "data_test/bLambda_PP.jld2"
    # filename_1 = "data/bLambda_PP.jld2"
    # filename_import = "data_test/bLambda_PP_1.jld2"
    # f = jldopen(filename, "r")
    # f_1 = jldopen(filename_1, "r")
    # f_ver2 = jldopen(filename_import, "a+")
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 kamp, ver3 = f[key]
    #                 kamp_1, ver3_1 = f_1[key]
    #                 kamp_import = cat(kamp, kamp_1; dims=1)
    #                 ver3_import = cat(ver3, ver3_1; dims=1)
    #                 # println(kamp_import/kF)
    #                 # println(ver3_import)
    #                 if haskey(f_ver2, key)
    #                     @warn("replacing existing data for $key")
    #                     delete!(f_ver2, key)
    #                 end
    #                 f_ver2[key] = (kamp_import, ver3_import)
    #             end
    #         end
    #     end
    # end
    # # _Lver3_direct channel
    # filename = "data_test/bLambda_Lver3_direct.jld2"
    # filename_1 = "data/bLambda_Lver3_direct.jld2"
    # filename_import = "data_test/bLambda_Lver3_direct_1.jld2"
    # f = jldopen(filename, "r")
    # f_1 = jldopen(filename_1, "r")
    # f_ver2 = jldopen(filename_import, "a+")
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 kamp, ver3 = f[key]
    #                 kamp_1, ver3_1 = f_1[key]
    #                 kamp_import = cat(kamp, kamp_1; dims=1)
    #                 ver3_import = cat(ver3, ver3_1; dims=1)
    #                 if haskey(f_ver2, key)
    #                     @warn("replacing existing data for $key")
    #                     delete!(f_ver2, key)
    #                 end
    #                 f_ver2[key] = (kamp_import, ver3_import)
    #             end
    #         end
    #     end
    # end

    # # _Lver3_exchange channel
    # filename = "data_test/bLambda_Lver3_exchange.jld2"
    # filename_1 = "data/bLambda_Lver3_exchange.jld2"
    # filename_import = "data_test/bLambda_Lver3_exchange_1.jld2"
    # f = jldopen(filename, "r")
    # f_1 = jldopen(filename_1, "r")
    # f_ver2 = jldopen(filename_import, "a+")
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 kamp, ver3 = f[key]
    #                 kamp_1, ver3_1 = f_1[key]
    #                 kamp_import = cat(kamp, kamp_1; dims=1)
    #                 ver3_import = cat(ver3, ver3_1; dims=1)
    #                 if haskey(f_ver2, key)
    #                     @warn("replacing existing data for $key")
    #                     delete!(f_ver2, key)
    #                 end
    #                 f_ver2[key] = (kamp_import, ver3_import)
    #             end
    #         end
    #     end
    # end

    # # _Lver3_bubble channel
    # filename = "data_test/bLambda_Lver3_bubble.jld2"
    # filename_1 = "data/bLambda_Lver3_bubble.jld2"
    # filename_import = "data_test/bLambda_Lver3_bubble_1.jld2"
    # f = jldopen(filename, "r")
    # f_1 = jldopen(filename_1, "r")
    # f_ver2 = jldopen(filename_import, "a+")
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 kamp, ver3 = f[key]
    #                 kamp_1, ver3_1 = f_1[key]
    #                 kamp_import = cat(kamp, kamp_1; dims=1)
    #                 ver3_import = cat(ver3, ver3_1; dims=1)
    #                 if haskey(f_ver2, key)
    #                     @warn("replacing existing data for $key")
    #                     delete!(f_ver2, key)
    #                 end
    #                 f_ver2[key] = (kamp_import, ver3_import)
    #             end
    #         end
    #     end
    # end

    # # import sigma
    # filename_sigma  = "data_test/sigma_RG.jld2"
    # filename_sigma_1  = "data/sigma_RG.jld2"
    # filename_sigma_import  = "data_test/sigma_RG_1.jld2"
    # datadict_import = Dict{eltype(partition),Any}()
    # f = jldopen(filename_sigma, "r")
    # f_1 = jldopen(filename_sigma_1, "r")
    # f_ver2 = jldopen(filename_sigma_import, "a+")
    # for (fi, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para_data = ParaMC(key)
    #         para_data_Fs0 = ElectronLiquidRG.get_para(para_data, -0.0)
    #         if para_data_Fs0 == _para
    #             if isapprox(para_data.Fs, F) == true
    #                 # println("Fs=$Fs_idx, key=$(key)")
    #                 # _para = ElectronLiquidRG.get_para(para, Fs_idx)
    #                 ngrid, kgrid, sigma = f[key]
    #                 ngrid_1, kgrid_1, sigma_1 = f_1[key]
    #                 kgrid_import = cat(kgrid, kgrid_1; dims=1)
    #                 # println(kgrid_import/kF)
    #                 # println(sigma[(1,0,0)])
    #                 println(size(sigma[(1,0,0)]), size(sigma_1[(1,0,0)]))
    #                 datadict_import[(1,0,0)] = cat(sigma[(1,0,0)], sigma_1[(1,0,0)]; dims=2)
    #                 sigma_import = datadict_import
    #                 # println(sigma_import[(1,0,0)])
    #                 # println(key)
    #                 if haskey(f_ver2, key)
    #                     @warn("replacing existing data for $key")
    #                     delete!(f_ver2, key)
    #                 end
    #                 f_ver2[key] = (ngrid, kgrid_import, sigma_import)
    #             end
    #         end
    #     end
    # end
end