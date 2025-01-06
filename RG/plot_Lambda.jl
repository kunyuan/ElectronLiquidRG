using ElectronLiquidRG
using ElectronLiquid
using CompositeGrids
using Measurements
using JLD2
using DelimitedFiles
using DataFrames
using CSV
using PyPlot, PyCall
using LinearAlgebra
using Lehmann
using CurveFit


dim = 3
rs = [5.0,]
mass2 = [1e-2]
beta = [40.0]
order = [1,]
isDynamic = true

function getdata(filename)
    df = DataFrame(CSV.File(filename))
    data = []
    for idx in axes(df, 1)
        data_Λ = zeros(Float64, size(df, 2))
        for jdx in axes(df, 2)
            val_data = measurement(df[idx, jdx]).val
            data_Λ[jdx] = val_data
        end
        push!(data, data_Λ)
    end
    return data
end

function compute_derivative(x, y, finegrid)
    y_prime = zeros(Float64, length(finegrid))
    for (idx, x_grid) in enumerate(finegrid)
        y_prime[idx] = Interp.differentiate1D(y, x, x_grid)
    end
    
    return y_prime
end

# derivative of bLambda/z1Lambda throught fitting piecewise
function derivative_fit(x, y, xgrid, n_fit)
    a_poly = poly_fit(x, y.*(x.^2), n_fit)
    dy_fit = zeros(Float64, length(xgrid))
    y_fit = zeros(Float64, length(xgrid))
    for  (jdx, kΛ) in enumerate(xgrid)
        data = 0.0
        for kdx = 1:n_fit+1
            y_fit[jdx] += a_poly[kdx] * kΛ^(kdx - 3)
            # println(a_poly[kdx] * F^(kdx-1))
        end
        for kdx = 1:n_fit
            data += a_poly[kdx+1] * kΛ^(kdx - 1) *kdx
            # println(a_poly[kdx] * F^(kdx-1))
        end
        dy_fit[jdx] = data / kΛ^2 - 2.0*y_fit[jdx] / kΛ
    end
    return y_fit, dy_fit
end


for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)
    _para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
    kF = _para.kF
    # Λgrid = [1.0] *_para.kF
    # Λgrid = SimpleG.Uniform([1.0*_para.kF, 3.0*_para.kF],21)
    Λgrid = ElectronLiquidRG.CompositeGrid.LogDensedGrid(:cheb, [1.0 * _para.kF, 8.0 * _para.kF], [_para.kF,], 6, 0.02 * _para.kF, 8)
    label_Λ = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,28,29,30,31,32,33,34,35,36,37,38,39,40] 
    # label_Λ = cat(label_Λ, [cutoff_Λ:length(Λgrid)]; dims=1)
    println(length(label_Λ))
    Λgrid_measurement = zeros(Float64, length(label_Λ))
    for (idx, idx_Λ) in enumerate(label_Λ)
        Λgrid_measurement[idx] = Λgrid[idx_Λ]
    end
    println(Λgrid_measurement/kF)
    Λgrid = Λgrid_measurement

    # Fs = SimpleG.Uniform([-0.4,0.0],41)
    Fs = SimpleG.Uniform([-0.8,0.0],81)

    # data_aΛ_iter = DataFrame(CSV.File("csv_ver4/aLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # data_bΛ_iter = DataFrame(CSV.File("csv_ver4/bLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # data_dbΛ_iter = DataFrame(CSV.File("csv_ver4/dbLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # data_z1Λ_iter = DataFrame(CSV.File("csv_ver4/z1Lambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # data_dz1Λ_iter = DataFrame(CSV.File("csv_ver4/dz1Lambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # data_uΛ_iter = DataFrame(CSV.File("csv_ver4/uLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # data_uΛ_Large = DataFrame(CSV.File("csv_ver4/uLambda_Large_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # data_fsΛ_iter = DataFrame(CSV.File("csv_ver4/fsLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # data_uΛ_1_iter = DataFrame(CSV.File("csv_ver4/uLambda_1_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # data_uΛ_2_iter = DataFrame(CSV.File("csv_ver4/uLambda_2_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))

    data_aΛ_iter = DataFrame(CSV.File("csv_test/aLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    data_bΛ_iter = DataFrame(CSV.File("csv_test/bLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    data_dbΛ_iter = DataFrame(CSV.File("csv_test/dbLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    data_z1Λ_iter = DataFrame(CSV.File("csv_test/z1Lambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    data_dz1Λ_iter = DataFrame(CSV.File("csv_test/dz1Lambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    data_uΛ_iter = DataFrame(CSV.File("csv_test/uLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    data_uΛ_Large = DataFrame(CSV.File("csv_test/uLambda_Large_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    data_fsΛ_iter = DataFrame(CSV.File("csv_test/fsLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    data_uΛ_1_iter = DataFrame(CSV.File("csv_test/uLambda_1_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    data_uΛ_2_iter = DataFrame(CSV.File("csv_test/uLambda_2_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))

    iter_plot = [50]
    

    # # test bLambad fitting
    # for (jdx,iter) in enumerate(iter_plot)
    #     NΛ = length(Λgrid)
    #     dbΛ = zeros(Float64, length(Λgrid))
    #     bΛ = zeros(Float64, length(Λgrid))
    #     for idx in eachindex(Λgrid)
    #         # println("idx=$(idx)")
    #         bΛ[idx] =  data_bΛ_iter[idx, iter]
    #         dbΛ[idx] =  data_dbΛ_iter[idx, iter]
    #     end

    #     n_fit_Λ = 4
    #     cutoff_Λ = 16
    #     Λ_fit_1 = SimpleG.Uniform([Λgrid[1], Λgrid[cutoff_Λ]], 501)
    #     Λ_fit_2 = SimpleG.Uniform([Λgrid[cutoff_Λ], Λgrid[end]], 501)
    #     Λ_fit = cat(Λ_fit_1, Λ_fit_2[2:end]; dims=1)
    #     bΛ_fit_1, dbΛ_1 = derivative_fit(Λgrid[1:cutoff_Λ], bΛ[1:cutoff_Λ], Λ_fit_1, n_fit_Λ)
    #     bΛ_fit_2, dbΛ_2 = derivative_fit(Λgrid[cutoff_Λ:end], bΛ[cutoff_Λ:end], Λ_fit_2, n_fit_Λ)
    #     bΛ_fit = cat(bΛ_fit_1, bΛ_fit_2[2:end]; dims=1)
    #     # bΛ_fit , dbΛ_fit = derivative_fit(Λgrid/kF, bΛ, Λ_fit/kF, n_fit_Λ)

    #     dbΛ_interp = zeros(Float64, length(Λgrid))
    #     bΛ_int = ones(Float64, length(Λgrid)) * bΛ[NΛ]
    #     for idx = 1:NΛ-1
    #         Λgrid_integrand = SimpleG.Arbitrary(Λgrid[idx:NΛ]/kF)
    #         bΛ_int[idx] = -Interp.integrate1D(dbΛ[idx:NΛ], Λgrid_integrand) + bΛ_int[NΛ]
    #     end

    #     Λgrid_diff = SimpleG.Arbitrary(Λgrid/kF)
    #     for (idx, kΛ_diff) in enumerate(Λgrid_diff)
    #         dbΛ_interp[idx] = Interp.differentiate1D(bΛ, Λgrid_diff, kΛ_diff)
    #     end

    #     bΛ_int_interp = ones(Float64, length(Λgrid)) * bΛ[NΛ]
    #     for idx = 1:NΛ-1
    #         Λgrid_integrand = SimpleG.Arbitrary(Λgrid[idx:NΛ]/kF)
    #         bΛ_int_interp[idx] = -Interp.integrate1D(dbΛ_interp[idx:NΛ], Λgrid_integrand) + bΛ_int[NΛ]
    #     end
    #     figure(1)
    #     scatter(Λgrid/_para.kF, bΛ, label="data")
    #     plot(Λ_fit/_para.kF, bΛ_fit, label="fit")
    #     legend()
    #     xlabel("\$\\Lambda/k_F\$")
    #     ylabel("\$ b_\\Lambda\$")
    #     xscale("log")
    #     # ylabel("\$\\delta b_\\Lambda\$")
    #     savefig("figure_test/bLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).pdf")
    #     close()
    #     figure(2)
    #     scatter(Λgrid/_para.kF, bΛ, label="data")
    #     plot(Λgrid/_para.kF, bΛ_int, label="int_fit")
    #     # plot(Λgrid/_para.kF, bΛ_int_interp, label="int_interp")
    #     legend()
    #     xlabel("\$\\Lambda/k_F\$")
    #     ylabel("\$ b_\\Lambda\$")
    #     xscale("log")
    #     # ylabel("\$\\delta b_\\Lambda\$")
    #     savefig("figure_test/deltabLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).pdf")
    #     close()
    # end

    # plot bLambda
    for (jdx,iter) in enumerate(iter_plot)
        NΛ = length(Λgrid)
        bΛ = zeros(Float64, length(Λgrid))
        FsΛ = data_fsΛ_iter[:,iter]
        println(FsΛ)
    end
end