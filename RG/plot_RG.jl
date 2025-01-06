using ElectronLiquidRG
using ElectronLiquid
using CompositeGrids
using Measurements
using JLD2
using DelimitedFiles
using DataFrames
using CSV
using PyPlot, PyCall

# pyimport("scienceplots")
style = PyPlot.matplotlib."style"
# style.use(["science", "std-colors"])
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")


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
    # Λgrid = range(1.0*_para.kF, 3.0*_para.kF,21)
    Λgrid = ElectronLiquidRG.CompositeGrid.LogDensedGrid(:cheb, [1.0 * _para.kF, 8.0 * _para.kF], [_para.kF,], 6, 0.02 * _para.kF, 8)
    # Λgrid = ElectronLiquidRG.CompositeGrid.LogDensedGrid(:gauss, [1.0 * _para.kF, 5 * _para.kF], [_para.kF,], 6, 0.01 * _para.kF, 6)
    # println(size(Λgrid))
    # asgrid = zeros(Float64, length(Λgrid))
    # Fs = ElectronLiquidRG.treelevel_RG(_para, Λgrid, asgrid)
    # println(Fs)
    Fs = SimpleG.Uniform([-0.4,0.0],41).grid 
    # Λgrid = ElectronLiquidRG.CompositeGrid.LogDensedGrid(:cheb, [1.0 * _para.kF, 8.0 * _para.kF], [_para.kF,], 12, 0.02 * _para.kF, 6)
    
    # Fs = [-0.2]
    # # Fs = Fs[1:3] 
    # data_aΛ_iter = DataFrame(CSV.File("aLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # data_bΛ_iter = DataFrame(CSV.File("bLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # data_uΛ_iter = DataFrame(CSV.File("uLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # data_fsΛ_iter = DataFrame(CSV.File("fsLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))

    # iter_plot = [1,2,5,8,10]
    # # plot aLambda_iter
    # figure(1)
    # for (jdx,iter) in enumerate(iter_plot)
    #     # println("iter=$(iter)")
    #     aΛ = zeros(Float64, length(Λgrid))
    #     for idx in eachindex(aΛ)
    #         # println("idx=$(idx)")
    #         data = measurement(data_aΛ_iter[idx, iter])
    #         aΛ[idx] = data.val 
    #     end
    #     plot(Λgrid/_para.kF, aΛ, label="iter=$(iter)")
    # end
    # legend()
    # xlabel("\$\\Lambda/k_F\$")
    # ylabel("\$a_\\Lambda\$")
    # savefig("figure/aLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).pdf")
    # close()
    # # plot bLambda_iter
    # figure(2)
    # for (jdx,iter) in enumerate(iter_plot)
    #     # println("iter=$(iter)")
    #     aΛ = zeros(Float64, length(Λgrid))
    #     for idx in eachindex(aΛ)
    #         # println("idx=$(idx)")
    #         data = measurement(data_bΛ_iter[idx, iter])
    #         aΛ[idx] = data.val 
    #     end
        
    #     plot(Λgrid/_para.kF, aΛ, label="iter=$(iter)")
    # end
    # legend()
    # xlabel("\$\\Lambda/k_F\$")
    # ylabel("\$b_\\Lambda\$")
    # savefig("figure/bLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).pdf")
    # close()

    # # plot uLambda_iter
    # figure(3)
    # for (jdx,iter) in enumerate(iter_plot)
    #     # println("iter=$(iter)")
    #     aΛ = zeros(Float64, length(Λgrid))
    #     for idx in eachindex(aΛ)
    #         # println("idx=$(idx)")
    #         data = measurement(data_uΛ_iter[idx, iter])
    #         aΛ[idx] = data.val 
    #     end
        
    #     plot(Λgrid/_para.kF, aΛ, label="iter=$(iter)")
        
    # end
    # legend()
    # xlabel("\$\\Lambda/k_F\$")
    # ylabel("\$u_\\Lambda\$")
    # savefig("figure/uLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).pdf")
    # close()
    # # plot fsLambda_iter
    # figure(4)
    # for (jdx,iter) in enumerate(iter_plot)
    #     # println("iter=$(iter)")
    #     aΛ = zeros(Float64, length(Λgrid))
    #     for idx in eachindex(aΛ)
    #         # println("idx=$(idx)")
    #         data = measurement(data_fsΛ_iter[idx, iter])
    #         aΛ[idx] = data.val 
    #     end
        
    #     plot(Λgrid/_para.kF, aΛ, label="iter=$(iter)")
    # end
    # legend()
    # xlabel("\$\\Lambda/k_F\$")
    # ylabel("\$F_s^\\Lambda\$")
    # savefig("figure/fsLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).pdf")
    # close()
end