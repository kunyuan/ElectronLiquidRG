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

# scipy = pyimport("scipy")


dim = 3
rs = [2.0,]
mass2 = [1e-2,]
beta = [40.0]
order = [1,]
isDynamic = true

function data_spline(filename, Λgrid, Fs)
    df = DataFrame(CSV.File(filename))
    # println(size(df))
    df_spline = []
    for (idx, kΛ) in enumerate(Λgrid)
        data_k = zeros(Float64, length(Fs))
        for (jdx, F) in enumerate(Fs)
            data = measurement(df[idx, jdx])
            data_k[jdx] = data.val
        end
        cs = scipy.interpolate.PchipInterpolator(Fs, data_k)
        # cs = scipy.interpolate.CubicSpline(Fs, data_k)
        push!(df_spline, cs)
    end
    return df_spline
end

function getdata(filename)
    df = DataFrame(CSV.File(filename))
    data = []
    data_e = []
    for idx in axes(df, 1)
        data_Λ = zeros(Float64, size(df, 2))
        data_e_Λ = zeros(Float64, size(df, 2))
        for jdx in axes(df, 2)
            val_data = measurement(df[idx, jdx]).val
            data_Λ[jdx] = val_data
            data_e_Λ[jdx] = measurement(df[idx, jdx]).err
        end
        push!(data, data_Λ)
        push!(data_e, data_e_Λ)
    end
    return data, data_e
end


function data_fit(data, Fs, n_fit)
    coeff_poly = []
    for idx in eachindex(data)
        a_poly = poly_fit(Fs, data[idx], n_fit)
        push!(coeff_poly, a_poly)
    end
    return coeff_poly
end

function getdata_fit(coeff, idx, F, n_fit)
    a_poly = coeff[idx]
    data = 0.0
    for jdx = 1:n_fit+1
        data += a_poly[jdx] * F^(jdx - 1)
    end
    return data
end

function extractdata(data, idx_i, idx_f)
    data_1 = []
    for (idx, data_F) in enumerate(data)
        # println(data_F)
        data_1_F = data_F[idx_i:idx_f]
        # println(data_1_F)
        push!(data_1, data_1_F)
    end
    return data_1
end

function data_fit_piecewise(data, Fs, cutoff, n_fit)
    coeff_poly = []
    Fs_1 = Fs[1:cutoff]
    Fs_2 = Fs[cutoff:end]
    for idx in eachindex(data)
        data_F = data[idx]
        data_F_1 = data_F[1:cutoff]
        data_F_2 = data_F[cutoff:end]
        a_poly = []
        push!(a_poly, poly_fit(Fs_1, data_F_1, n_fit))
        push!(a_poly, poly_fit(Fs_2, data_F_2, n_fit))
        push!(coeff_poly, a_poly)
    end
    return coeff_poly
end

function getdata_fit_piecewise(coeff, idx, F, Fs, cutoff, n_fit)
    # println(coeff[idx])
    if F < Fs[cutoff]
        a_poly = coeff[idx][1]
    else
        a_poly = coeff[idx][2]
    end
    data = 0.0
    for jdx = 1:n_fit+1
        data += a_poly[jdx] * F^(jdx - 1)
    end
    return data
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

    if _rs == 1.0
        Fs = SimpleG.Uniform([-0.99, 0.0], 100)
        Fs1 = SimpleG.Uniform([-1.5, -1.01], 50).grid
        Fs = cat(Fs1, Fs; dims=1)
    elseif _rs == 2.0
        Fs = SimpleG.Uniform([-0.99, 0.0], 100)
        Fs1 = SimpleG.Uniform([-1.5, -1.01], 50).grid
        Fs = cat(Fs1, Fs; dims=1)
    elseif _rs == 3.0  
        Fs = SimpleG.Uniform([-0.99, 0.0], 100)
        Fs1 = SimpleG.Uniform([-1.6, -1.01], 60).grid
        Fs = cat(Fs1, Fs; dims=1)
    elseif _rs == 5.0
        Fs = SimpleG.Uniform([-0.99, 0.0], 100).grid
        Fs1 = SimpleG.Uniform([-2.4, -1.01], 140).grid
        Fs = cat(Fs1, Fs; dims=1)
    end
    # Fs = Fs[1:3] 

    # zfactor 
    data_z1, data_z1_e = getdata("data_csv/zfactor_RG_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    # println(data_z1_e[1])
    # println(data_z1)
    # data_z1 = data_z1[cutoff_Λ : length(Λgrid)]
    # println(data_z1)
    # aLambda 
    data_aΛ_PP_PHE, data_aΛ_PP_PHE_e = getdata("data_csv/aLambda_PP_PHE_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    # data_aΛ_PP_PHE = data_aΛ_PP_PHE[cutoff_Λ : length(Λgrid)]
    data_aΛ_PH, data_aΛ_PH_e = getdata("data_csv/aLambda_PH_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    # data_aΛ_PH = data_aΛ_PH[cutoff_Λ : length(Λgrid)]
    # bLambda 
    data_bΛ_PP, data_bΛ_PP_e = getdata("data_csv/bLambda_PP_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    # data_bΛ_PP = data_bΛ_PP[cutoff_Λ : length(Λgrid)]
    data_bΛ_Lver3_direct, data_bΛ_Lver3_direct_e = getdata("data_csv/bLambda_Lver3_direct_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    # data_bΛ_Lver3_direct = data_bΛ_Lver3_direct[cutoff_Λ : length(Λgrid)]
    data_bΛ_Lver3_exchange, data_bΛ_Lver3_exchange_e = getdata("data_csv/bLambda_Lver3_exchange_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    # data_bΛ_Lver3_exchange = data_bΛ_Lver3_exchange[cutoff_Λ : length(Λgrid)]
    data_bΛ_Lver3_bubble, data_bΛ_Lver3_bubble_e = getdata("data_csv/bLambda_Lver3_bubble_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    # data_bΛ_Lver3_bubble = data_bΛ_Lver3_bubble[cutoff_Λ : length(Λgrid)]
    
    idx_i = 1
    idx_f = length(Fs)
    # idx_i = 1
    # idx_f = 121

    Fs = Fs[idx_i:idx_f]
    data_z1 = extractdata(data_z1, idx_i, idx_f)
    data_z1_e = extractdata(data_z1_e, idx_i, idx_f)
    data_aΛ_PP_PHE = extractdata(data_aΛ_PP_PHE, idx_i, idx_f)
    data_aΛ_PP_PHE_e = extractdata(data_aΛ_PP_PHE_e, idx_i, idx_f)
    data_aΛ_PH = extractdata(data_aΛ_PH, idx_i, idx_f)
    data_aΛ_PH_e = extractdata(data_aΛ_PH_e, idx_i, idx_f)
    data_bΛ_PP = extractdata(data_bΛ_PP, idx_i, idx_f)
    data_bΛ_PP_e = extractdata(data_bΛ_PP_e, idx_i, idx_f)
    data_bΛ_Lver3_direct = extractdata(data_bΛ_Lver3_direct, idx_i, idx_f)
    data_bΛ_Lver3_direct_e = extractdata(data_bΛ_Lver3_direct_e, idx_i, idx_f)
    data_bΛ_Lver3_exchange = extractdata(data_bΛ_Lver3_exchange, idx_i, idx_f)
    data_bΛ_Lver3_exchange_e = extractdata(data_bΛ_Lver3_exchange_e, idx_i, idx_f)
    data_bΛ_Lver3_bubble = extractdata(data_bΛ_Lver3_bubble, idx_i, idx_f)
    data_bΛ_Lver3_bubble_e = extractdata(data_bΛ_Lver3_bubble_e, idx_i, idx_f)

    
    n_fit = 4
    # coeff_z1 = data_fit(data_z1, Fs, n_fit)
    # coeff_aΛ_PP_PHE = data_fit(data_aΛ_PP_PHE, Fs, n_fit)
    # coeff_aΛ_PH = data_fit(data_aΛ_PH, Fs, n_fit)
    # coeff_bΛ_PP = data_fit(data_bΛ_PP, Fs, n_fit)
    # coeff_bΛ_Lver3_direct = data_fit(data_bΛ_Lver3_direct, Fs, n_fit)
    # coeff_bΛ_Lver3_exchange = data_fit(data_bΛ_Lver3_exchange, Fs, n_fit)
    # coeff_bΛ_Lver3_bubble = data_fit(data_bΛ_Lver3_bubble, Fs, n_fit)

    cutoff = 60
    println("F=$(Fs[60])")
    coeff_z1 = data_fit_piecewise(data_z1, Fs, cutoff, n_fit)
    coeff_aΛ_PP_PHE = data_fit_piecewise(data_aΛ_PP_PHE, Fs, cutoff, n_fit)
    coeff_aΛ_PH = data_fit_piecewise(data_aΛ_PH, Fs, cutoff, n_fit)
    coeff_bΛ_PP = data_fit_piecewise(data_bΛ_PP, Fs, cutoff, n_fit)
    coeff_bΛ_Lver3_direct = data_fit_piecewise(data_bΛ_Lver3_direct, Fs, cutoff, n_fit)
    coeff_bΛ_Lver3_exchange = data_fit_piecewise(data_bΛ_Lver3_exchange, Fs, cutoff, n_fit)
    coeff_bΛ_Lver3_bubble = data_fit_piecewise(data_bΛ_Lver3_bubble, Fs, cutoff, n_fit)

    Fs_grid = SimpleG.Uniform([Fs[1], Fs[end]], 501)
    
    idx_Λ = [1,5,10,15,20]
    # idx_Λ = [1]
    coeff = coeff_bΛ_Lver3_bubble
    data = data_bΛ_Lver3_bubble
    data_e = data_bΛ_Lver3_bubble_e
    for idx in idx_Λ
        data_fit_Λ = zeros(Float64, length(Fs_grid))
        data_Λ = data[idx]
        data_e_Λ = data_e[idx]
        println(length(data_e_Λ))
        for (jdx, _F) in enumerate(Fs_grid)
            # data_fit_Λ[jdx] = getdata_fit(coeff, idx, _F, n_fit)
            data_fit_Λ[jdx] = getdata_fit_piecewise(coeff, idx, _F, Fs, cutoff, n_fit)
        end
        
        plot(Fs_grid, data_fit_Λ, label="fit")
        # scatter(Fs, data_Λ, label="data")
        errorbar(Fs, data_Λ, yerr=data_e_Λ, fmt="o", label="data")
        xlabel("Fs")
        # ylabel("\$a_\\Lambda\$")
        # ylabel("z_1")
        # ylabel("\$R\\Pi\$")
        legend()
        show()
        # savefig("aLambda_PP_PHE_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).pdf")
        close()
    end

end