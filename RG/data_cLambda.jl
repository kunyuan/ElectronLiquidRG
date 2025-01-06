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
rs = [2.0,]
mass2 = [1e-2,]
beta = [40.0]
order = [1,]
isDynamic = true

for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)
    _para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
    kF = _para.kF
    # Λgrid = [1.0, 1.1, 1.2] *_para.kF
    Λgrid = SimpleG.Uniform([1.0*_para.kF, 8.0*_para.kF],501)
    

    cΛs = zeros(Float64, length(Λgrid))
    cΛu = zeros(Float64, length(Λgrid))
    for (idx, kΛ) in enumerate(Λgrid)
        data_u, data_s = ElectronLiquidRG.c_coeff_pp(_para, kΛ, kΛ)
        data_s = data_s/_para.NF
        cΛs[idx] = data_s
        cΛu[idx] = ElectronLiquidRG.cLambda_Pi(_para; kΛ) / 2.0
    end
    Λgrid_measurement = ElectronLiquidRG.CompositeGrid.LogDensedGrid(:cheb, [1.0 * _para.kF, 8.0 * _para.kF], [_para.kF,], 6, 0.02 * _para.kF, 8)
    cutoff_Λ = 1
    Λgrid_measurement = Λgrid_measurement[cutoff_Λ:end]
    # println(length(Λgrid_measurement))
    dcΛs = zeros(Float64, length(Λgrid_measurement))
    dcΛu = zeros(Float64, length(Λgrid_measurement))
    for (idx, kΛ) in enumerate(Λgrid_measurement)
        dcΛs[idx] = Interp.differentiate1D(cΛs, Λgrid, kΛ)
        dcΛu[idx] = Interp.differentiate1D(cΛu, Λgrid, kΛ)
    end
    filename_data_dcLambda = "data_csv/dcLambda_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).csv"
    dcLambda = DataFrame()
    colname = "dcLambda_s"
    dcLambda[!, colname] = dcΛs * kF
    colname = "dcLambda_u"
    dcLambda[!, colname] = dcΛu * kF
    CSV.write(filename_data_dcLambda, dcLambda)
end