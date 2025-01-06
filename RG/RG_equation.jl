# only export the u/f of the last iter given the different inital condition on ~8kF

using ElectronLiquidRG
using ElectronLiquid
using CompositeGrids
using Measurements
using JLD2
using DelimitedFiles
using DataFrames
using CSV
using LinearAlgebra
using PyPlot, PyCall
using CurveFit
using SpecialFunctions

dim = 3
rs = [2.0,]
mass2 = [1e-2,]
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

function data_fit_F_piecewise(data, Fs, cutoff, n_fit)
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
    return dy_fit
end

function derive_counterterm_R(para, kamp, dz, bLambda_bubble)
    ct = true
    θgrid = ElectronLiquidRG.CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 32)
    φgrid = SimpleG.Uniform([0.1π, π],100)
    avgR_theta = zeros(Float64, length(θgrid))
    avgR_phi = zeros(Float64, length(φgrid))
    avgPi_theta = zeros(Float64, length(θgrid))
    avgPi_phi = zeros(Float64, length(φgrid))
    # avgRPi_theta = zeros(Float64, length(θgrid))
    # avgRPi_phi = zeros(Float64, length(φgrid))
    for (idx, _theta) in enumerate(θgrid)
        for (jdx, _phi) in enumerate(φgrid)
            extK1 = [sin(0.5_theta), cos(0.5_theta), 0] * kamp
            extK2 = [-sin(0.5_theta), cos(0.5_theta), 0] * kamp
            extK3 = [sin(0.5_theta)*cos(_phi),cos(0.5_theta),sin(0.5_theta)*sin(_phi)] * kamp
            # qabs = sqrt(dot(extK1-extK2, extK1-extK2))
            qabs = sqrt(dot(extK1-extK3, extK1-extK3))
            avgR_phi[jdx] = UEG.KOstatic(qabs, para; ct=ct)
            avgPi_phi[jdx] = UEG.polarKW(qabs, 0, para)
            # avgRPi_phi[jdx] = avgR_phi[jdx]*avgPi_phi[jdx]
            
        end
        avgR_theta[idx] = Interp.integrate1D(avgR_phi, φgrid) / (π)
        avgPi_theta[idx] = Interp.integrate1D(avgPi_phi, φgrid) / (π)
        # avgRPi_theta[idx] = Interp.integrate1D(avgRPi_phi, φgrid) / (π)
    end
    avgPi = Ver4.Legrendre(0, avgPi_theta, θgrid)
    avgR = Ver4.Legrendre(0, avgR_theta, θgrid)
    # avgRPi = Ver4.Legrendre(0, avgRPi_theta, θgrid)
    # println("avgRPi=$(avgRPi), bLambda_bubble=$(bLambda_bubble)")


    # qs = [kamp * sqrt(2*(1 - cos(θ))) for θ in θgrid.grid]
    # Wp = zeros(Float64, length(qs))
    # for (qi, q) in enumerate(qs)
    #     Wp[qi] =UEG.KOstatic(q, para; ct=ct)
    # end
    # avgR_ex = Ver4.Legrendre(0, Wp, θgrid)
    
    f1 = para.fs
    # δR = 2.0*dz*avgR - 2*f1*avgRPi - f1^2*avgPi
    δR = 2.0*dz*avgR - 2*f1*2.0*bLambda_bubble - f1^2*avgPi

    # if kamp/para.kF>7.94
    #     println("k=$(kamp/para.kF)kF, avgRPi=$(avgRPi), $(2.0*dz*avgR*para.NF), $(- 2*f1*avgRPi*para.NF), $(- f1^2*avgPi*para.NF)")
    # end
    δR = δR*para.NF
    return δR
end

function derive_z1Λ(para,
    Λgrid,
    Fs,
    coeff_z1,
    Fs_grid,
    cutoff_F,
    n_fit
)
    z1Λ = zeros(Float64, length(Λgrid))
    for (idx, F) in enumerate(Fs)
        data = getdata_fit_piecewise(coeff_z1, idx, F, Fs_grid, cutoff_F, n_fit)
        z1Λ[idx] = data
    end
    return z1Λ
end

function derive_bLambda(para,
    Λgrid,
    Fs,
    coeff_bΛ_PP,
    coeff_bΛ_Lver3_direct,
    coeff_bΛ_Lver3_exchange,
    coeff_bΛ_Lver3_bubble,
    Fs_grid,
    cutoff_F,
    n_fit
)
    function b_tail(para, k)
        B = para.me * (para.e0)^2 * (π - 2.0)
        return B/k
    end
    bΛ = zeros(Float64, length(Λgrid))
    for (idx, F) in enumerate(Fs)
        data = getdata_fit_piecewise(coeff_bΛ_PP, idx, F, Fs_grid, cutoff_F, n_fit)
        bΛ[idx] += 2.0*data
        
        data = getdata_fit_piecewise(coeff_bΛ_Lver3_direct, idx, F, Fs_grid, cutoff_F, n_fit)
        bΛ[idx] += 2.0*data

        data = getdata_fit_piecewise(coeff_bΛ_Lver3_exchange, idx, F, Fs_grid, cutoff_F, n_fit)
        bΛ[idx] += 2.0*data

        data = getdata_fit_piecewise(coeff_bΛ_Lver3_bubble, idx, F, Fs_grid, cutoff_F, n_fit)
        bΛ[idx] = bΛ[idx] - 2.0*data
        # println("Lb=$(data)")
        # if idx == length(Λgrid)
        #     println("Lver3_bubble=$data")
        # end
    end
    # idx = length(Λgrid)
    # println("k=$(Λgrid[idx]/para.kF)kF, F=$(Fs[idx]), bΛ=$(bΛ[idx])")
    return bΛ
end

function derive_aLambda(para,
    Λgrid,
    Fs,
    coeff_aΛ_PP_PHE, 
    coeff_aΛ_PH,
    coeff_z1,
    coeff_bΛ_Lver3_bubble,
    Fs_grid,
    cutoff_F,
    n_fit
) 
    # println(Fsgrid)
    aΛ = zeros(Float64, length(Λgrid))
    for (idx, F) in enumerate(Fs)
        # PP+PHE channel
        data = getdata_fit_piecewise(coeff_aΛ_PP_PHE, idx, F, Fs_grid, cutoff_F, n_fit)
        aΛ[idx] += data
        # if idx == 1
        # println("k=$(Λgrid[idx]/para.kF)")
        # println("PP+PHE=$(data)")
        # end
        # PH channel
        data = getdata_fit_piecewise(coeff_aΛ_PH, idx, F, Fs_grid, cutoff_F, n_fit)
        aΛ[idx] += data
        # println(data_aΛ_PH[idx])
        # println("Fs=$(F)")
        # if idx == 1
        # println("PH=$(data)")
        # end

        # δR
        para0 = ElectronLiquidRG.get_para(para, F)
        kΛ = Λgrid[idx]
        
        data = getdata_fit_piecewise(coeff_z1, idx, F, Fs_grid, cutoff_F, n_fit)
        dz = data

        data = getdata_fit_piecewise(coeff_bΛ_Lver3_bubble, idx, F, Fs_grid, cutoff_F, n_fit)
        bLambda_bubble = data

        δR = derive_counterterm_R(para0, kΛ, dz, bLambda_bubble)
        aΛ[idx] = aΛ[idx] - δR
    end
    return aΛ
end


function fit_LargeLambda(filename, _para)
    data_DF = DataFrame()
    f = jldopen(filename, "r")
    # println(keys(f))
    Λgrid = []
    
    for key in keys(f)
        para = ParaMC(key)
        if para == _para
            F = para.Fs
            println("key=$(key)")
            kamp, n, datadict = f[key]
            # println(datadict)
            Λgrid = kamp
            data = datadict[(1, 0, 0)]
            colname = "Fs_$(F)"
            data_DF[!, colname] = data
        end
    end
    
    data = zeros(Float64, length(Λgrid))
    for idx in eachindex(Λgrid)
        d = data_DF[idx , 1]
        # println(d.val)
        data[idx] = d.val
    end
    
    Λgrid_measurement = Λgrid
    data_measurement = data
    # println(Λgrid_measurement/_para.kF)
    a, b = power_fit(Λgrid_measurement, data_measurement)
    # println(a, ",", b)
    return a, b
end

function derive_aLambda_Large(para,
    Λgrid,
    a_fit,
    b_fit
) 
    aΛ = zeros(Float64, length(Λgrid))
    for (idx, kΛ) in enumerate(Λgrid)
        aΛ[idx] = a_fit * kΛ^b_fit
    end
    return aΛ
end

function derive_dbLambda_Large(para,
    Λgrid,
    B
) 
    bΛ = zeros(Float64, length(Λgrid))
    dbΛ = zeros(Float64, length(Λgrid))
    for (idx, kΛ) in enumerate(Λgrid)
        bΛ[idx] = B / kΛ
        dbΛ[idx] = - B / kΛ^2 * para.kF
    end
    return bΛ, dbΛ
end

function derive_equ_u_Large(para, Λgrid, aΛ, dbΛ, dcΛ, uΛ; mix=0.8) 
    NΛ = length(Λgrid)
    
    duΛ_2 = dbΛ .* uΛ
    duΛ_3 = dcΛ .* uΛ.^2
    duΛ = duΛ_2 + duΛ_3
    # println(duΛ)
    uΛ_derive = zeros(Float64, NΛ)
    uΛ_derive[NΛ] = uΛ[NΛ]
    # println("mix=$(mix)")
    # println(Λgrid[3:NΛ])
    for jdx = 1:NΛ-1
        Λgrid_integrand = SimpleG.Arbitrary(Λgrid[jdx:NΛ])
        uΛ_derive[jdx] = aΛ[jdx] - aΛ[NΛ] - Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand) + uΛ_derive[NΛ]
        # uΛ_derive[jdx] = - Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand) + uΛ_derive[NΛ]
    end
    uΛ_derive = mix * uΛ + (1 - mix) *uΛ_derive
    return uΛ_derive
end

function derive_equ_u_Large_v2(para, Λgrid, aΛ, bΛ, dcΛ, uΛ; mix=0.8) 
    NΛ = length(Λgrid)
    
    uΛ_2 = bΛ .* uΛ
    duΛ = dcΛ .* uΛ.^2

    uΛ_derive = zeros(Float64, NΛ)
    uΛ_derive[NΛ] = uΛ[NΛ]
    # println("mix=$(mix)")
    # println(Λgrid[3:NΛ])
    for jdx = 1:NΛ-1
        Λgrid_integrand = SimpleG.Arbitrary(Λgrid[jdx:NΛ])
        uΛ_derive[jdx] = aΛ[jdx] - aΛ[NΛ] + uΛ_2[jdx] - uΛ_2[NΛ] - Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand) + uΛ_derive[NΛ]
    end
    uΛ_derive = mix * uΛ + (1 - mix) *uΛ_derive
    return uΛ_derive
end

function derive_equ_u_Large_v3(para, Λgrid, aΛ, bΛ, cΛ, uΛ; mix=0.8) 
    NΛ = length(Λgrid)
    
    uΛ_2 = bΛ .* uΛ + cΛ .* uΛ.^2

    uΛ_derive = zeros(Float64, NΛ)
    uΛ_derive[NΛ] = uΛ[NΛ]
    # println("mix=$(mix)")
    # println(Λgrid[3:NΛ])
    for jdx = 1:NΛ-1
        uΛ_derive[jdx] = aΛ[jdx] - aΛ[NΛ] + uΛ_2[jdx] - uΛ_2[NΛ] + uΛ_derive[NΛ]
        # println(jdx,",", uΛ_derive[jdx])
    end
    uΛ_derive = mix * uΛ + (1 - mix) *uΛ_derive
    return uΛ_derive
end

function derive_equ_u_Large_v4(para, Λgrid, aΛ, bΛ, dbΛ, cΛ, dcΛ, uΛ; mix=0.8) 
    NΛ = length(Λgrid)
    duΛ = zeros(Float64, NΛ)
    for idx in eachindex(Λgrid)
        duΛ[idx] = (dbΛ[idx] * uΛ[idx] + dcΛ[idx] * (uΛ[idx])^2) / (1.0 - bΛ[idx] - 2.0 * cΛ[idx] * uΛ[idx])
    end

    uΛ_derive = zeros(Float64, NΛ)
    uΛ_derive[NΛ] = uΛ[NΛ]
    # println("mix=$(mix)")
    # println(Λgrid[3:NΛ])
    for jdx = 1:NΛ-1
        Λgrid_integrand = SimpleG.Arbitrary(Λgrid[jdx:NΛ])
        uΛ_derive[jdx] = aΛ[jdx] - aΛ[NΛ] - Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand) + uΛ_derive[NΛ]
        # println(jdx,",", uΛ_derive[jdx])
    end
    uΛ_derive = mix * uΛ + (1 - mix) *uΛ_derive
    return uΛ_derive
end

function derive_equ_u_Large_analytic(para, Λgrid, B, a_c) 
    kF = para.kF
    coeff_c = 6π * para.NF / para.me / kF
    coeff_B = 2.0 * B / kF
    println("a_c=$(a_c), c=$(coeff_c), B=$(coeff_B)")
    uΛ = zeros(Float64, length(Λgrid))
    for (idx, kΛ) in enumerate(Λgrid)
        uΛ[idx] = coeff_c / (1.0/a_c*exp(coeff_B/kΛ) - kΛ - coeff_B*exp(coeff_B/kΛ)*expinti(-coeff_B/kΛ))
        # uΛ[idx] = coeff_c / (- kΛ - coeff_B*exp(coeff_B/kΛ)*expinti(-coeff_B/kΛ))
        # uΛ[idx] = coeff_c / (1.0/a_c * exp(coeff_B/kΛ) - kΛ - coeff_B*log(coeff_B/kΛ))
        # uΛ[idx] = coeff_c / (- kΛ - coeff_B*log(coeff_B/kΛ))
    end
    return uΛ
end

function test_equ_u_Large(para, Λgrid, aΛ, dbΛ, dcΛ, uΛ)
    NΛ = length(Λgrid)
    duΛ_2 = dbΛ .* uΛ
    duΛ_3 = dcΛ .* uΛ.^2
    duΛ = duΛ_2 + duΛ_3

    Λgrid_difference = SimpleG.Arbitrary(Λgrid)
    duΛ_test = zeros(Float64, NΛ)
    for (idx, kΛ) in enumerate(Λgrid)
        duΛ_test[idx] = Interp.differentiate1D(uΛ, Λgrid_difference, kΛ)
    end
    plot(Λgrid, duΛ-duΛ_test)
    show()
    return duΛ_test
end

function derive_equ_u(para, Λgrid, aΛ, bΛ, dcΛs, dcΛu, z1Λ, uΛ, n_fit_Λ, cutoff_Λ; mix=0.95) 
    NΛ = length(Λgrid)
    dbΛ = zeros(Float64, NΛ)
    dz1Λ = zeros(Float64, NΛ)
    
    dbΛ_1 = derivative_fit(Λgrid[1:cutoff_Λ], bΛ[1:cutoff_Λ], Λgrid[1:cutoff_Λ], n_fit_Λ)
    dbΛ_2 = derivative_fit(Λgrid[cutoff_Λ:end], bΛ[cutoff_Λ:end], Λgrid[cutoff_Λ:end], n_fit_Λ)
    dbΛ = cat(dbΛ_1, dbΛ_2[2:end]; dims=1)
    dz1Λ_1 = derivative_fit(Λgrid[1:cutoff_Λ], z1Λ[1:cutoff_Λ], Λgrid[1:cutoff_Λ], n_fit_Λ)
    dz1Λ_2 = derivative_fit(Λgrid[cutoff_Λ:end], z1Λ[cutoff_Λ:end], Λgrid[cutoff_Λ:end], n_fit_Λ)
    dz1Λ = cat(dz1Λ_1, dz1Λ_2[2:end]; dims=1)
    # x = SimpleG.Arbitrary(Λgrid)
    # dbΛ = compute_derivative(x, bΛ, Λgrid)
    # dz1Λ = compute_derivative(x, z1Λ, Λgrid)
    dcΛ = dcΛs + dcΛu
    # duΛ_2 = dbΛ .* uΛ
    duΛ_2 = (dbΛ + 4.0*dz1Λ) .* uΛ
    duΛ_3 = dcΛ .* uΛ.^2
    duΛ = duΛ_2 + duΛ_3

    uΛ_derive = zeros(Float64, NΛ)
    uΛ_derive[NΛ] = uΛ[NΛ]
    uΛ_derive_1 = zeros(Float64, NΛ)
    uΛ_derive_2 = zeros(Float64, NΛ)
    # println("mix=$(mix)")
    # println(Λgrid[3:NΛ])
    for jdx = 1:NΛ-1
        Λgrid_integrand = SimpleG.Arbitrary(Λgrid[jdx:NΛ])
        uΛ_derive[jdx] = aΛ[jdx] - aΛ[NΛ] - Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand) + uΛ_derive[NΛ]
        uΛ_derive_1[jdx] = aΛ[jdx] - aΛ[NΛ]
        uΛ_derive_2[jdx] = Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand)
        # if jdx == NΛ-1
        #     println("u_1=$(aΛ[jdx] - aΛ[NΛ]), u_2=$(Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand)), du=$(duΛ[jdx]),$(duΛ[NΛ])")
        # end
    end
    # println(uΛ_derive[NΛ])
    println("u_last = $(uΛ[1])\t u_derive = $(uΛ_derive[1])")
    uΛ_derive = mix * uΛ + (1 - mix) *uΛ_derive
    return dbΛ, dz1Λ, uΛ_derive_1, uΛ_derive_2, uΛ_derive
end

function derive_equ_u_v2(para, Λgrid, aΛ, bΛ, cΛs, cΛu, z1Λ, uΛ, n_fit_Λ, cutoff_Λ; mix=0.95)
    NΛ = length(Λgrid)
    deltagamma = zeros(Float64, NΛ)

    dbΛ = zeros(Float64, NΛ)
    dz1Λ = zeros(Float64, NΛ)
    dbΛ_1 = derivative_fit(Λgrid[1:cutoff_Λ], bΛ[1:cutoff_Λ], Λgrid[1:cutoff_Λ], n_fit_Λ)
    dbΛ_2 = derivative_fit(Λgrid[cutoff_Λ:end], bΛ[cutoff_Λ:end], Λgrid[cutoff_Λ:end], n_fit_Λ)
    dbΛ = cat(dbΛ_1, dbΛ_2[2:end]; dims=1)
    dz1Λ_1 = derivative_fit(Λgrid[1:cutoff_Λ], z1Λ[1:cutoff_Λ], Λgrid[1:cutoff_Λ], n_fit_Λ)
    dz1Λ_2 = derivative_fit(Λgrid[cutoff_Λ:end], z1Λ[cutoff_Λ:end], Λgrid[cutoff_Λ:end], n_fit_Λ)
    dz1Λ = cat(dz1Λ_1, dz1Λ_2[2:end]; dims=1)
    duΛ = 2.0*dz1Λ.*uΛ

    for (idx, kΛ) in enumerate(Λgrid)
        deltagamma[idx] = aΛ[idx] + bΛ[idx]*uΛ[idx] + 2.0*z1Λ[idx]*uΛ[idx] + (cΛs[idx] + cΛu[idx]) * (uΛ[idx])^2
        # deltagamma[idx] = aΛ[idx] + bΛ[idx]*uΛ[idx] + (cΛs[idx] + cΛu[idx]) * (uΛ[idx])^2
    end

    uΛ_derive = zeros(Float64, NΛ)
    uΛ_derive_1 = zeros(Float64, NΛ)
    uΛ_derive_2 = zeros(Float64, NΛ)
    uΛ_derive[NΛ] = uΛ[NΛ]
    for jdx = 1:NΛ-1
        Λgrid_integrand = SimpleG.Arbitrary(Λgrid[jdx:NΛ])
        uΛ_derive[jdx] = deltagamma[jdx] - deltagamma[NΛ] - Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand) + uΛ_derive[NΛ]
        uΛ_derive_1[jdx] = deltagamma[jdx] - deltagamma[NΛ]
        uΛ_derive_2[jdx] = Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand)
        # if jdx == NΛ-1
        #     println("u_1=$(aΛ[jdx] - aΛ[NΛ]), u_2=$(Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand)), du=$(duΛ[jdx]),$(duΛ[NΛ])")
        # end
    end
    println("u_last = $(uΛ[1])\t u_derive = $(uΛ_derive[1])")
    uΛ_derive = mix * uΛ + (1 - mix) *uΛ_derive
    return dbΛ, dz1Λ, uΛ_derive_1, uΛ_derive_2, uΛ_derive
end 

function derive_equ_u_v3(para, Λgrid, aΛ, bΛ, dcΛs, dcΛu, z1Λ, uΛ, n_fit_Λ, cutoff_Λ; mix=0.95)
    NΛ = length(Λgrid)
    deltagamma = zeros(Float64, NΛ)

    dbΛ = zeros(Float64, NΛ)
    dz1Λ = zeros(Float64, NΛ)
    dbΛ_1 = derivative_fit(Λgrid[1:cutoff_Λ], bΛ[1:cutoff_Λ], Λgrid[1:cutoff_Λ], n_fit_Λ)
    dbΛ_2 = derivative_fit(Λgrid[cutoff_Λ:end], bΛ[cutoff_Λ:end], Λgrid[cutoff_Λ:end], n_fit_Λ)
    dbΛ = cat(dbΛ_1, dbΛ_2[2:end]; dims=1)
    dz1Λ_1 = derivative_fit(Λgrid[1:cutoff_Λ], z1Λ[1:cutoff_Λ], Λgrid[1:cutoff_Λ], n_fit_Λ)
    dz1Λ_2 = derivative_fit(Λgrid[cutoff_Λ:end], z1Λ[cutoff_Λ:end], Λgrid[cutoff_Λ:end], n_fit_Λ)
    dz1Λ = cat(dz1Λ_1, dz1Λ_2[2:end]; dims=1)
    dcΛ = dcΛs + dcΛu
    duΛ = 2.0*dz1Λ.*uΛ + dcΛ .* uΛ.^2

    for (idx, kΛ) in enumerate(Λgrid)
        deltagamma[idx] = aΛ[idx] + bΛ[idx]*uΛ[idx] + 2.0*z1Λ[idx]*uΛ[idx]
        # deltagamma[idx] = aΛ[idx] + bΛ[idx]*uΛ[idx] + (cΛs[idx] + cΛu[idx]) * (uΛ[idx])^2
    end

    uΛ_derive = zeros(Float64, NΛ)
    uΛ_derive_1 = zeros(Float64, NΛ)
    uΛ_derive_2 = zeros(Float64, NΛ)
    uΛ_derive[NΛ] = uΛ[NΛ]
    for jdx = 1:NΛ-1
        Λgrid_integrand = SimpleG.Arbitrary(Λgrid[jdx:NΛ])
        uΛ_derive[jdx] = deltagamma[jdx] - deltagamma[NΛ] - Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand) + uΛ_derive[NΛ]
        uΛ_derive_1[jdx] = deltagamma[jdx] - deltagamma[NΛ]
        uΛ_derive_2[jdx] = Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand)
        # if jdx == NΛ-1
        #     println("u_1=$(aΛ[jdx] - aΛ[NΛ]), u_2=$(Interp.integrate1D(duΛ[jdx:NΛ], Λgrid_integrand)), du=$(duΛ[jdx]),$(duΛ[NΛ])")
        # end
    end
    println("u_last = $(uΛ[1])\t u_derive = $(uΛ_derive[1])")
    uΛ_derive = mix * uΛ + (1 - mix) *uΛ_derive
    return dbΛ, dz1Λ, uΛ_derive_1, uΛ_derive_2, uΛ_derive
end 

function derive_f_Gamma3(para, Λgrid, Fs, uΛ, coeff_Lver3_direct, coeff_z1, Fs_grid, cutoff_F, n_fit; mix=0.8)
    Fs_derive = zeros(Float64, length(Fs))
    for (idx, F) in enumerate(Fs)
        para0 = ElectronLiquidRG.get_para(para, F)
        kamp = Λgrid[idx]
        # kamp = Λgrid[1]
        ct = true
        θgrid = ElectronLiquidRG.CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 32)
        φgrid = SimpleG.Uniform([0.1π, π], 100)
        avgPi_theta = zeros(Float64, length(θgrid))
        avgPi_phi = zeros(Float64, length(φgrid))
        avgRPi_theta = zeros(Float64, length(θgrid))
        avgRPi_phi = zeros(Float64, length(φgrid))
        for (idx, _theta) in enumerate(θgrid)
            for (jdx, _phi) in enumerate(φgrid)
                extK1 = [sin(0.5_theta), cos(0.5_theta), 0] * kamp
                extK2 = [-sin(0.5_theta), cos(0.5_theta), 0] * kamp
                extK3 = [sin(0.5_theta) * cos(_phi), cos(0.5_theta), sin(0.5_theta) * sin(_phi)] * kamp
                # qabs = sqrt(dot(extK1-extK2, extK1-extK2))
                qabs = sqrt(dot(extK1 - extK3, extK1 - extK3))
                avgPi_phi[jdx] = UEG.polarKW(qabs, 0, para0)
                avgRPi_phi[jdx] = UEG.KOstatic(qabs, para0; ct=ct) * avgPi_phi[jdx]

            end
            avgRPi_theta[idx] = Interp.integrate1D(avgRPi_phi, φgrid) / (π)
            avgPi_theta[idx] = Interp.integrate1D(avgPi_phi, φgrid) / (π)
        end
        avgPi = Ver4.Legrendre(0, avgPi_theta, θgrid)
        avgRPi = Ver4.Legrendre(0, avgRPi_theta, θgrid)

        bLambda_Lver3d = - getdata_fit_piecewise(coeff_Lver3_direct, idx, F, Fs_grid, cutoff_F, n_fit)
        z1 = getdata_fit_piecewise(coeff_z1, idx, F, Fs_grid, cutoff_F, n_fit)

        Fs_derive[idx] =  (bLambda_Lver3d + z1 + 0.5*uΛ[idx] * avgPi / para.NFstar ) / (avgPi / para.NFstar)
    end
    # println(Fs_derive)
    Fs_derive = mix*Fs + (1.0 - mix) * Fs_derive
    return Fs_derive
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
    cutoff_Λ = 1
    # println(label_Λ[cutoff_Λ])
    
    if _rs == 1.0
        Fs = SimpleG.Uniform([-0.99, 0.0], 100)
        Fs1 = SimpleG.Uniform([-1.5, -1.01], 50).grid
        Fs = cat(Fs1, Fs; dims=1)
        cutoff_F = 60
    elseif _rs == 2.0
        Fs = SimpleG.Uniform([-0.99, 0.0], 100)
        Fs1 = SimpleG.Uniform([-1.5, -1.01], 50).grid
        Fs = cat(Fs1, Fs; dims=1)
        cutoff_F = 60
    elseif _rs == 3.0  
        Fs = SimpleG.Uniform([-0.99, 0.0], 100)
        Fs1 = SimpleG.Uniform([-1.6, -1.01], 60).grid
        Fs = cat(Fs1, Fs; dims=1)
        cutoff_F = 80
    elseif _rs == 5.0
        Fs = SimpleG.Uniform([-0.99, 0.0], 100).grid
        Fs1 = SimpleG.Uniform([-2.4, -1.01], 140).grid
        Fs = cat(Fs1, Fs; dims=1)
        cutoff_F = 100
    end

    # zfactor 
    data_z1 = getdata("data_csv/zfactor_RG_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    # println(data_z1)
    data_z1 = data_z1[cutoff_Λ : length(Λgrid)]
    # println(data_z1)
    # aLambda 
    data_aΛ_PP_PHE = getdata("data_csv/aLambda_PP_PHE_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    data_aΛ_PP_PHE = data_aΛ_PP_PHE[cutoff_Λ : length(Λgrid)]
    data_aΛ_PH= getdata("data_csv/aLambda_PH_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    data_aΛ_PH = data_aΛ_PH[cutoff_Λ : length(Λgrid)]
    # bLambda 
    data_bΛ_PP = getdata("data_csv/bLambda_PP_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    data_bΛ_PP = data_bΛ_PP[cutoff_Λ : length(Λgrid)]
    data_bΛ_Lver3_direct = getdata("data_csv/bLambda_Lver3_direct_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    data_bΛ_Lver3_direct = data_bΛ_Lver3_direct[cutoff_Λ : length(Λgrid)]
    data_bΛ_Lver3_exchange = getdata("data_csv/bLambda_Lver3_exchange_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    data_bΛ_Lver3_exchange = data_bΛ_Lver3_exchange[cutoff_Λ : length(Λgrid)]
    data_bΛ_Lver3_bubble = getdata("data_csv/bLambda_Lver3_bubble_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    data_bΛ_Lver3_bubble = data_bΛ_Lver3_bubble[cutoff_Λ : length(Λgrid)]
    
    data_avgRPi = getdata("data_csv/avgRPi_static_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv")
    
    # fit
    n_fit_F = 4
    coeff_z1 = data_fit_F_piecewise(data_z1, Fs, cutoff_F, n_fit_F)
    coeff_aΛ_PP_PHE = data_fit_F_piecewise(data_aΛ_PP_PHE, Fs, cutoff_F, n_fit_F)
    coeff_aΛ_PH = data_fit_F_piecewise(data_aΛ_PH, Fs, cutoff_F, n_fit_F)
    coeff_bΛ_PP = data_fit_F_piecewise(data_bΛ_PP, Fs, cutoff_F, n_fit_F)
    coeff_bΛ_Lver3_direct = data_fit_F_piecewise(data_bΛ_Lver3_direct, Fs, cutoff_F, n_fit_F)
    coeff_bΛ_Lver3_exchange = data_fit_F_piecewise(data_bΛ_Lver3_exchange, Fs, cutoff_F, n_fit_F)
    coeff_bΛ_Lver3_bubble = data_fit_F_piecewise(data_bΛ_Lver3_bubble, Fs, cutoff_F, n_fit_F)


    filename_dcLambda = "data_csv/dcLambda_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"
    df = DataFrame(CSV.File(filename_dcLambda))
    dcΛs = zeros(Float64,length(Λgrid))
    dcΛu = zeros(Float64,length(Λgrid))
    for (idx, idx_Λ) in enumerate(label_Λ)
        dcΛs[idx] = df[idx_Λ, 1]
        dcΛu[idx] = df[idx_Λ, 2]
    end
    dcΛs = dcΛs[cutoff_Λ : length(Λgrid)]
    dcΛu = dcΛu[cutoff_Λ : length(Λgrid)]
    # dcΛs = -dcΛs
    # dcΛu = -dcΛu

    Λgrid = Λgrid[cutoff_Λ : length(Λgrid)]
    cΛs = zeros(Float64, length(Λgrid))
    cΛu = zeros(Float64, length(Λgrid))
    for (idx, kΛ) in enumerate(Λgrid)
        data_u, data_s = ElectronLiquidRG.c_coeff_pp(_para, kΛ, kΛ)
        data_s = data_s/_para.NF
        cΛs[idx] = data_s
        cΛu[idx] = ElectronLiquidRG.cLambda_Pi(_para; kΛ) / 2.0
    end
    # cΛs = -cΛs
    # cΛu = -cΛu
    
    # Large Lambda tail 
    filename_aLambda_Large = "data_csv/aLambda_PP_PHE_Large.jld2"
    Λgrid_measurement_Large = SimpleG.Log{Float64}([Λgrid[length(Λgrid)] , 100.0*kF], 101, 0.5*kF, true).grid 
    # Λgrid_measurement_Large = SimpleG.Log{Float64}([10.0*kF , 100.0*kF], 201, 0.5*kF, true).grid 
    
    B = _para.me * (_para.e0)^2 * (π - 2.0)
    # println("B=$(B/kF)kF")

    ac = [0.0]
    label_ac = 1
    a_c = ac[label_ac]
    # Λgrid_measurement_Large_Large = SimpleG.Log{Float64}([Λgrid_measurement_Large[end] , 10000.0*kF], 101, 0.5*kF, true).grid
    # uΛ_Large_Large = derive_equ_u_Large_analytic(_para, Λgrid_measurement_Large_Large/kF, B, a_c)
    # uΛ_Large_Large = derive_equ_u_Large_analytic(_para, Λgrid_measurement_Large/kF, B, a_c)
    # scatter(Λgrid_measurement_Large/kF, uΛ_Large_Large)
    # xscale("log")
    # show()
    # println(Λgrid_measurement_Large[1]/kF,",",uΛ_Large_Large[1])


    # boundary condition
    # ugrid = uΛ_Large_Large 
    # ugrid = ones(Float64, length(Λgrid_measurement_Large)) * uΛ_Large_Large[end]
    ugrid = ones(Float64, length(Λgrid_measurement_Large)) * 0.0

    FsΛ = zeros(Float64, length(Λgrid_measurement_Large))
    a_fit, b_fit = fit_LargeLambda(filename_aLambda_Large, _para)
    aΛ = derive_aLambda_Large(_para, Λgrid_measurement_Large, a_fit, b_fit)
    aΛ = - aΛ
    # aΛ = zeros(Float64, length(Λgrid_measurement_Large))
    println("aΛ=$(aΛ[1])")
    bΛ, dbΛ = derive_dbLambda_Large(_para, Λgrid_measurement_Large, B)
    dbΛ = 2.0*dbΛ
    dbΛ = - dbΛ
    bΛ = 2.0*bΛ
    bΛ = - bΛ
    # println(dbΛ)
    cΛ = _para.me / (6π *_para.NF) * Λgrid_measurement_Large
    dcΛ = ones(Float64, length(Λgrid_measurement_Large)) * _para.me/( 6π *_para.NF) *kF
    # dcΛ = - dcΛ
    uΛ = ugrid
    # println(uΛ)
    # uΛ = ugrid + aΛ
    # uΛ[length(Λgrid_measurement_Large)] = 0.0
    data_uΛ_Large_iter = DataFrame()
    uΛ_Large_final = []
    iter_Large = 100
    for iter = 1:iter_Large
        uΛ = derive_equ_u_Large(_para, Λgrid_measurement_Large/kF, aΛ, dbΛ, dcΛ, uΛ)
        # uΛ = derive_equ_u_Large_v2(_para, Λgrid_measurement_Large/kF, aΛ, bΛ, dcΛ, uΛ)
        # uΛ = derive_equ_u_Large_v3(_para, Λgrid_measurement_Large/kF, aΛ, bΛ, cΛ, uΛ)
        # uΛ = derive_equ_u_Large_v4(_para, Λgrid_measurement_Large/kF, aΛ, bΛ, dbΛ, cΛ, dcΛ, uΛ)

        # colname = "iter_$(iter)"
        # data_uΛ_Large_iter[!, colname] = uΛ
        println("iter=$(iter), u=$(uΛ[1])")
        if iter == iter_Large
            uΛ_Large_final = uΛ
        end
    end
    colname = "uLambda_Large"
    data_uΛ_Large_iter[!, colname] = uΛ_Large_final
    println(uΛ_Large_final[1],",", uΛ_Large_final[end])
    # du = test_equ_u_Large(_para, Λgrid_measurement_Large/kF, aΛ, dbΛ, dcΛ, uΛ_Large_Large)
    # du = test_equ_u_Large(_para, Λgrid_measurement_Large/kF, aΛ, dbΛ, dcΛ, uΛ)
    # plot(Λgrid_measurement_Large/kF,aΛ)
    # show()
    # close()

    # plot(Λgrid_measurement_Large/kF, uΛ_Large_final, label="iter")
    # # plot(Λgrid_measurement_Large/kF, uΛ_Large_Large, label="analytic")
    # xscale("log")
    # legend()
    # show()

    # df_u_inital = DataFrame(CSV.File("csv_final_u/uLambda_inital_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # df_u = DataFrame(CSV.File("csv_final_u/uLambda_final_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # df_f = DataFrame(CSV.File("csv_final_u/fLambda_final_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    df_f = DataFrame()
    df_u = DataFrame()
    df_ac = DataFrame()
    # df_f = DataFrame(CSV.File("csv_final_u/fLambda_final_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # df_u = DataFrame(CSV.File("csv_final_u/uLambda_final_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))
    # df_ac = DataFrame(CSV.File("csv_final_u/ac_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"))


    # u_inital = -0.5
    # uΛ = ones(Float64, length(Λgrid)) * u_inital
    uΛ = ones(Float64, length(Λgrid)) * uΛ_Large_final[1]
    ugrid = uΛ
    Ugrid = ugrid / _para.NF
    asgrid = _para.me/ 4π * uΛ
    data_aΛ_iter = DataFrame()
    data_bΛ_iter = DataFrame()
    data_dbΛ_iter = DataFrame()
    data_z1Λ_iter = DataFrame()
    data_dz1Λ_iter = DataFrame()
    data_uΛ_iter = DataFrame()
    data_uΛ_1_iter = DataFrame()
    data_uΛ_2_iter = DataFrame()
    data_fsΛ_iter = DataFrame()

    Λgrid_measurement = Λgrid/_para.kF
    n_fit_Λ = 4
    cutoff_Λ = 16
    println(Λgrid_measurement[cutoff_Λ])
    FsΛ = ElectronLiquidRG.treelevel_RG(_para, Λgrid, asgrid)
    
    # FsΛ = -0.8*ones(Float64, length(Λgrid))
    mix = 0.8
    # iteration
    for iter = 1:100
        # println("iter=$(iter)")
        # println("a=$(length(asgrid)), u=$(length(ugrid)) ")
        # FsΛ = ElectronLiquidRG.treelevel_RG(_para, Λgrid, asgrid)
        # FsΛ = ElectronLiquidRG.treelevel_RG_iter(_para, FsΛ, Λgrid, asgrid; mix=mix)
        FsΛ = derive_f_Gamma3(_para, Λgrid, FsΛ, uΛ, coeff_bΛ_Lver3_direct, coeff_z1, Fs, cutoff_F, n_fit_F; mix=mix)
        
        z1Λ = derive_z1Λ(_para, Λgrid, FsΛ, coeff_z1, Fs, cutoff_F, n_fit_F)
        # bΛ = derive_bLambda(_para, Λgrid, FsΛ, coeff_bΛ_PP, coeff_bΛ_Lver3_direct, coeff_bΛ_Lver3_exchange, coeff_bΛ_Lver3_bubble, n_fit_F)
        # aΛ = derive_aLambda(_para, Λgrid, FsΛ, Ugrid, coeff_aΛ_PP_PHE, coeff_aΛ_PH, coeff_z1, coeff_bΛ_Lver3_bubble, n_fit_F)
        bΛ = derive_bLambda(_para, Λgrid, FsΛ, coeff_bΛ_PP, coeff_bΛ_Lver3_direct, coeff_bΛ_Lver3_exchange, coeff_bΛ_Lver3_bubble, Fs, cutoff_F, n_fit_F)
        aΛ = derive_aLambda(_para, Λgrid, FsΛ, coeff_aΛ_PP_PHE, coeff_aΛ_PH, coeff_z1, coeff_bΛ_Lver3_bubble, Fs, cutoff_F, n_fit_F)
        aΛ = - aΛ
        bΛ = - bΛ
        # println("a=$aΛ")
        # println("b=$bΛ")
        # dbΛ, dz1Λ, uΛ_1, uΛ_2, uΛ = derive_equ_u(_para, Λgrid_measurement, aΛ, bΛ, dcΛs, dcΛu, z1Λ, ugrid, n_fit_Λ, cutoff_Λ; mix=mix)
        dbΛ, dz1Λ, uΛ_1, uΛ_2, uΛ = derive_equ_u_v2(_para, Λgrid_measurement, aΛ, bΛ, cΛs, cΛu, z1Λ, ugrid, n_fit_Λ, cutoff_Λ; mix=mix)
        # dbΛ, dz1Λ, uΛ_1, uΛ_2, uΛ = derive_equ_u_v3(_para, Λgrid_measurement, aΛ, bΛ, dcΛs, dcΛu, z1Λ, ugrid, n_fit_Λ, cutoff_Λ; mix=mix)
        ugrid = uΛ
        Ugrid = ugrid / _para.NF
        asgrid = _para.me/ 4π * uΛ
        colname = "iter_$(iter)"
        data_aΛ_iter[!, colname] = aΛ # [cutoff_Λ:end]
        data_bΛ_iter[!, colname] = bΛ # [cutoff_Λ:end]
        data_dbΛ_iter[!, colname] = dbΛ # [cutoff_Λ:end]
        data_z1Λ_iter[!, colname] = z1Λ # [cutoff_Λ:end]
        data_dz1Λ_iter[!, colname] = dz1Λ # [cutoff_Λ:end]
        # data_dcΛ_iter[!, colname] = dcΛ
        data_uΛ_iter[!, colname] = uΛ # [cutoff_Λ:end]
        data_uΛ_1_iter[!, colname] = uΛ_1
        data_uΛ_2_iter[!, colname] = uΛ_2
        data_fsΛ_iter[!, colname] = FsΛ # [cutoff_Λ:end]
        println("iter=$(iter), Fs=$(FsΛ[1]), a=$(aΛ[1]), uΛ=$(uΛ[1])")
        # plot(Λgrid_measurement, uΛ)
        # show()
        # println("iter=$(iter), Fs=$(FsΛ[cutoff_Λ]), a=$(aΛ[cutoff_Λ]), uΛ=$(uΛ[cutoff_Λ])")
    end

    delta_u_iter = zeros(Float64, size(data_uΛ_iter,2)-1)
    for jdx in eachindex(delta_u_iter)
        delta_u_iter[jdx] = data_uΛ_iter[1, jdx+1] - data_uΛ_iter[1, jdx]
    end
    figure(1)
    plot(1:size(data_uΛ_iter,2)-1, delta_u_iter)
    xlabel("iter")
    show()
    close()
    fileData = open("csv_final_u/uLambda_final_ac=$(a_c)_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).dat","w")
    writedlm(fileData, delta_u_iter)

    plot(Λgrid_measurement, data_fsΛ_iter[:,end])
    ylabel("f")
    show()
    close()

    plot(Λgrid_measurement, data_uΛ_iter[:,end])
    ylabel("u")
    show()
    close()

    
    colname = "ac=$(a_c)"
    df_f[!, colname] = data_fsΛ_iter[:, end]
    df_u[!, colname] = data_uΛ_iter[:, end]
    df_ac[!, colname] = [a_c]
    
    CSV.write("csv_final_u/uLambda_final_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", df_u)
    CSV.write("csv_final_u/fLambda_final_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", df_f)
    CSV.write("csv_final_u/ac_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", df_ac)

    # CSV.write("csv_test/uLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_uΛ_iter)
    # CSV.write("csv_test/uLambda_1_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_uΛ_1_iter)
    # CSV.write("csv_test/uLambda_2_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_uΛ_2_iter)
    # CSV.write("csv_test/uLambda_Large_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_uΛ_Large_iter)
    # CSV.write("csv_test/aLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_aΛ_iter)
    # CSV.write("csv_test/bLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_bΛ_iter)
    # CSV.write("csv_test/dbLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_dbΛ_iter)
    # CSV.write("csv_test/z1Lambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_z1Λ_iter)
    # CSV.write("csv_test/dz1Lambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_dz1Λ_iter)
    # CSV.write("csv_test/fsLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_fsΛ_iter)

    # CSV.write("csv_up/uLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_uΛ_iter)
    # CSV.write("csv_up/uLambda_1_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_uΛ_1_iter)
    # CSV.write("csv_up/uLambda_2_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_uΛ_2_iter)
    # CSV.write("csv_up/uLambda_Large_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_uΛ_Large_iter)
    # CSV.write("csv_up/aLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_aΛ_iter)
    # CSV.write("csv_up/bLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_bΛ_iter)
    # CSV.write("csv_up/dbLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_dbΛ_iter)
    # CSV.write("csv_up/z1Lambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_z1Λ_iter)
    # CSV.write("csv_up/dz1Lambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_dz1Λ_iter)
    # CSV.write("csv_up/fsLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_fsΛ_iter)

    # CSV.write("csv_temp/uLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_uΛ_iter)
    # CSV.write("csv_temp/uLambda_1_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_uΛ_1_iter)
    # CSV.write("csv_temp/uLambda_2_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_uΛ_2_iter)
    # CSV.write("csv_temp/uLambda_Large_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_uΛ_Large_iter)
    # CSV.write("csv_temp/aLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_aΛ_iter)
    # CSV.write("csv_temp/bLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_bΛ_iter)
    # CSV.write("csv_temp/dbLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_dbΛ_iter)
    # CSV.write("csv_temp/z1Lambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_z1Λ_iter)
    # CSV.write("csv_temp/dz1Lambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_dz1Λ_iter)
    # CSV.write("csv_temp/fsLambda_iter_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv", data_fsΛ_iter)

end