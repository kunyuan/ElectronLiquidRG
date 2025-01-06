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

for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)
    _para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
    kF = _para.kF
    # Λgrid = [1.0] *_para.kF
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

    # _PP channel
    filename = "data/bLambda_PP.jld2"
    filename_data = "data_csv/bLambda_PP_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"
    data_bLambda = DataFrame()
    f = jldopen(filename, "r")
    for (idx, F) in enumerate(Fs1)
        for key in keys(f)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    kamp, ver3 = f[key]
                    colname = "Fs_$(F)"
                    data_bLambda[!, colname] = ver3
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
    #                 kamp, ver3 = f[key]
    #                 data_import = []
    #                 for (idx, idx_Λ) in enumerate(label_Λ)
    #                     push!(data_import, ver3[idx_Λ])
    #                 end
    #                 colname = "Fs_$(F)"
    #                 data_bLambda[!, colname] = data_import
    #             end
    #         end
    #     end
    # end
    CSV.write(filename_data, data_bLambda)

    # _Lver3_direct channel
    filename = "data/bLambda_Lver3_direct.jld2"
    filename_data = "data_csv/bLambda_Lver3_direct_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"
    data_bLambda = DataFrame()
    f = jldopen(filename, "r")
    for (idx, F) in enumerate(Fs1)
        for key in keys(f)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    kamp, ver3 = f[key]
                    colname = "Fs_$(F)"
                    data_bLambda[!, colname] = ver3
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
    #                 kamp, ver3 = f[key]
    #                 data_import = []
    #                 for (idx, idx_Λ) in enumerate(label_Λ)
    #                     push!(data_import, ver3[idx_Λ])
    #                 end
    #                 colname = "Fs_$(F)"
    #                 data_bLambda[!, colname] = data_import
    #             end
    #         end
    #     end
    # end
    CSV.write(filename_data, data_bLambda)

    # _Lver3_exchange channel
    filename = "data/bLambda_Lver3_exchange.jld2"
    filename_data = "data_csv/bLambda_Lver3_exchange_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"
    data_bLambda = DataFrame()
    f = jldopen(filename, "r")
    for (idx, F) in enumerate(Fs1)
        for key in keys(f)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    kamp, ver3 = f[key]
                    colname = "Fs_$(F)"
                    data_bLambda[!, colname] = ver3
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
    #                 kamp, ver3 = f[key]
    #                 data_import = []
    #                 for (idx, idx_Λ) in enumerate(label_Λ)
    #                     push!(data_import, ver3[idx_Λ])
    #                 end
    #                 colname = "Fs_$(F)"
    #                 data_bLambda[!, colname] = data_import
    #             end
    #         end
    #     end
    # end
    CSV.write(filename_data, data_bLambda)

    # _Lver3_bubble channel
    filename = "data/bLambda_Lver3_bubble.jld2"
    filename_data = "data_csv/bLambda_Lver3_bubble_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"
    data_bLambda = DataFrame()
    f = jldopen(filename, "r")
    for (idx, F) in enumerate(Fs1)
        for key in keys(f)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    kamp, ver3 = f[key]
                    colname = "Fs_$(F)"
                    data_bLambda[!, colname] = ver3
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
    #                 kamp, ver3 = f[key]
    #                 data_import = []
    #                 for (idx, idx_Λ) in enumerate(label_Λ)
    #                     push!(data_import, ver3[idx_Λ])
    #                 end
    #                 colname = "Fs_$(F)"
    #                 data_bLambda[!, colname] = data_import
    #             end
    #         end
    #     end
    # end
    CSV.write(filename_data, data_bLambda)
end