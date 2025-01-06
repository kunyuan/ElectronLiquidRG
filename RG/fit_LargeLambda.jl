using ElectronLiquidRG
using ElectronLiquid
using CompositeGrids
using LinearAlgebra, Measurements
using JLD2
using DelimitedFiles
using DataFrames
using CSV
using CurveFit
using PyPlot, PyCall

# pyimport("scienceplots")
style = PyPlot.matplotlib."style"
# style.use(["science", "std-colors"])
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")



neval = 1e7
dim = 3
rs = [1.0,]
mass2 = [1e-2,]
beta = [40.0]
order = [1,]
isDynamic = true

for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)
    _para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
    kF = _para.kF
    # Λgrid = [1.0] *_para.kF
    # Λgrid = ElectronLiquidRG.CompositeGrid.LogDensedGrid(:cheb, [1.0 * _para.kF, 8.0 * _para.kF], [_para.kF,], 6, 0.02 * _para.kF, 8)
    # Λgrid = range(1.0*_para.kF, 3.0*_para.kF,21)
    # Λgrid = ElectronLiquidRG.CompositeGrid.LogDensedGrid(:gauss, [1.0 * _para.kF, 5 * _para.kF], [_para.kF,], 6, 0.01 * _para.kF, 6)
    # println(size(Λgrid))
    # asgrid = zeros(Float64, length(Λgrid))
    # Fs = ElectronLiquidRG.treelevel_RG(_para, Λgrid, asgrid)
    # println(Fs)
    # Fs = SimpleG.Uniform([-0.4,0.0],41).grid 
    # Fs = Fs[1:3] 
    Fs = [-0.0]
    # label_idx = 36
    # println(Λgrid[label_idx]/kF)
    # println(Λgrid[40]/kF)
    
    # aLambda
    # filename_aPP  = "data/aLambda_PP_PHE.jld2"
    # data_DF = DataFrame()
    # f = jldopen(filename_aPP, "r")
    # # println(keys(f))
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
    #                 colname = "Fs_$(F)"
    #                 data_DF[!, colname] = data
    #             end
    #         end
    #     end
    # end

    filename_aPP  = "aLambda_PP_PHE_Large.jld2"
    data_DF = DataFrame()
    f = jldopen(filename_aPP, "r")
    # println(keys(f))
    Λgrid_1 = []
    for (idx, F) in enumerate([-0.0])
        for key in keys(f)
            println(key)
            para = ParaMC(key)
            para0 = ElectronLiquidRG.get_para(para, -0.0)
            if para0 == _para
                # println("true")
                if isapprox(para.Fs, F) == true
                    # println("Fs=$(F), key_cluster=$(key)")
                    kamp, n, datadict = f[key]
                    # println(datadict)
                    Λgrid_1 = kamp
                    data = datadict[(1, 0, 0)]
                    colname = "Fs_$(F)"
                    data_DF[!, colname] = data
                end
            end
        end
    end
    data_2 = zeros(Float64, length(Λgrid_1))
    for idx = 1 : length(Λgrid_1)
        d = data_DF[idx , 1]
        # println(d.val)
        data_2[idx] = d.val
    end
    println(Λgrid_1 )
    # println(data_2)
    # println(Λgrid_1/kF)
    Λgrid_measurement = Λgrid_1
    data_measurement = data_2
    println(Λgrid_measurement/kF)
    println(data_measurement)
    asgrid = zeros(Float64, length(Λgrid_measurement))
    Fs = ElectronLiquidRG.treelevel_RG(_para, Λgrid_measurement, asgrid)
    println("f= $(data_measurement+Fs)")

    # # _PP channel
    # filename = "data/bLambda_PP.jld2"
    # f = jldopen(filename, "r")
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
                
    #             if isapprox(para.Fs, F) == true
    #                 println(F)
    #                 kamp, ver3 = f[key]
    #                 colname = "Fs_$(F)"
    #                 data_DF[!, colname] = ver3
    #             end
    #         end
    #     end
    # end
    # data_1 = zeros(Float64, length(Λgrid))
    # for idx = 1 : length(Λgrid)
    #     d = data_DF[idx , 1]
    #     # println(d.val)
    #     data_1[idx] = d.val
    # end
    # data_1 = data_1[label_idx : length(Λgrid)]
    # println(data_1)

    
    # println(data_measurement)
    # println(length(data_measurement),",", length(Λgrid_measurement))
    # Λgrid_measurement = Λgrid[label_idx : length(Λgrid)] / kF
    # println(Λgrid_measurement)
    # a, b = power_fit(Λgrid_measurement, data_measurement)
    # println(a, ",", b)
    # figure(1)
    # plot(Λgrid_measurement, data_measurement, label="data")
    # xgrid = SimpleG.Uniform([8.0,10.0], 41)
    # ygrid = a* xgrid.^b
    # plot(xgrid, ygrid, label="fit")
    # legend()
    # xlabel("\$\\Lambda/k_F\$")
    # ylabel("\$a_\\Lambda\$")
    # show()
    # close()

    # # _Lver3_direct channel
    # filename = "data/bLambda_Lver3_direct.jld2"
    # filename_data = "bLambda_Lver3_direct.csv"
    # data_bLambda = DataFrame()
    # f = jldopen(filename, "r")
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 kamp, ver3 = f[key]
    #                 colname = "Fs_$(F)"
    #                 data_bLambda[!, colname] = ver3
    #             end
    #         end
    #     end
    # end
    # CSV.write(filename_data, data_bLambda)

    # # _Lver3_exchange channel
    # filename = "data/bLambda_Lver3_exchange.jld2"
    # filename_data = "bLambda_Lver3_exchange.csv"
    # data_bLambda = DataFrame()
    # f = jldopen(filename, "r")
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 kamp, ver3 = f[key]
    #                 colname = "Fs_$(F)"
    #                 data_bLambda[!, colname] = ver3
    #             end
    #         end
    #     end
    # end
    # CSV.write(filename_data, data_bLambda)

    # # _Lver3_bubble channel
    # filename = "data/bLambda_Lver3_bubble.jld2"
    # filename_data = "bLambda_Lver3_bubble.csv"
    # data_bLambda = DataFrame()
    # f = jldopen(filename, "r")
    # for (idx, F) in enumerate(Fs)
    #     for key in keys(f)
    #         para = ParaMC(key)
    #         para0 = ElectronLiquidRG.get_para(para, -0.0)
    #         if para0 == _para
    #             # println("true")
    #             if isapprox(para.Fs, F) == true
    #                 kamp, ver3 = f[key]
    #                 colname = "Fs_$(F)"
    #                 data_bLambda[!, colname] = ver3
    #             end
    #         end
    #     end
    # end
    # CSV.write(filename_data, data_bLambda)
end