using ElectronLiquidRG
using ElectronLiquid
using CompositeGrids

neval = 1e7
dim = 3
rs = [5.0,]
mass2 = [1e-2,]
beta = [40.0]
order = [1,]
# Fs = SimpleG.Uniform([-0.99,0.0],100).grid 
Fs = [-0.0]
isDynamic = true

println(Fs)
for (_rs, _mass2, _beta, _order, _Fs) in Iterators.product(rs, mass2, beta, order, Fs)
    _para = UEG.ParaMC(rs=_rs, beta=_beta, Fs = _Fs, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
    kF = _para.kF
    # println(UEG.short(_para))
    Λgrid = [1.0, 1.1, 1.2, 1.4, 1.5] *_para.kF
    # Λgrid = CompositeGrid.LogDensedGrid(:cheb, [1.0 * _para.kF, 8.0 * _para.kF], [_para.kF,], 6, 0.02 * _para.kF, 8)
    asgrid = zeros(Float64, length(Λgrid))
    FsΛ = ElectronLiquidRG.treelevel_RG(_para, Λgrid, asgrid)
    println(FsΛ)

    # θgrid = SimpleG.Uniform([0.0,π],21).grid
    # φgrid = SimpleG.Uniform([0.0,π],21).grid


    # filename_aPP  = "aLambda_PP_PHE.jld2"
    # filename_aPH  = "aLambda_PH.jld2"

    # filename_aPP_test  = "aLambda_PP_PHE_test.jld2"
    # filename_aPH_test  = "aLambda_PH_test.jld2"
    
    # filename_bPP  = "bLambda_PP.jld2"
    # filename_bLver3_direct  = "bLambda_Lver3_direct.jld2"
    # filename_bLver3_exchange  = "bLambda_Lver3_exchange.jld2"
    # filename_bLver3_bubble  = "bLambda_Lver3_bubble.jld2"
    # filename_sigma = "sigma_RG.jld2"
    # # filename1 = "sigma1_test_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).jld2"
    # # filename2 = "sigma2_test_rs=$(_rs)_mass2=$(_mass2)_beta=$(_beta).jld2"

    # KinL, KoutL, KinR = zeros(16), zeros(16), zeros(16)
    # KinL[1], KoutL[2], KinR[3] = 1.0, 1.0, 1.0
    


    #     # println("PH channel:")
    #     transferLoop = KinL - KoutL
    #     # ElectronLiquidRG.aLambda_avg(_para; kamp=Λgrid, n=[0, 0, 0, 0], theta=θgrid, phi=φgrid, neval = neval, filename=filename_aPH, filter=[ElectronLiquidRG.Parquet.NoHartree, ElectronLiquidRG.Parquet.NoBubble], channels=[ElectronLiquidRG.Parquet.PHr], transferLoop=transferLoop)
    #     ElectronLiquidRG.aLambda_thetaphi_avg(_para; kamp=Λgrid, n=[0, 0, 0, 0], theta=θgrid, phi=φgrid, neval = neval, filename=filename_aPH_test, filter=[ElectronLiquidRG.Parquet.NoHartree, ElectronLiquidRG.Parquet.NoBubble], channels=[ElectronLiquidRG.Parquet.PHr], transferLoop=transferLoop)
    #     # println("PHE+PP channel:")
    #     transferLoop = -KinR + KoutL
    #     # ElectronLiquidRG.aLambda_avg(_para; kamp=Λgrid, n=[0, 0, 0, 0], theta=θgrid, phi=φgrid, neval = neval, filename=filename_aPP, filter=[ElectronLiquidRG.Parquet.NoHartree, ElectronLiquidRG.Parquet.NoBubble, ElectronLiquidRG.Parquet.Proper], channels=[ElectronLiquidRG.Parquet.PPr, ElectronLiquidRG.Parquet.PHEr], transferLoop=transferLoop)  
    #     ElectronLiquidRG.aLambda_thetaphi_avg(_para; kamp=Λgrid, n=[0, 0, 0, 0], theta=θgrid, phi=φgrid, neval = neval, filename=filename_aPP_test, filter=[ElectronLiquidRG.Parquet.NoHartree, ElectronLiquidRG.Parquet.NoBubble, ElectronLiquidRG.Parquet.Proper], channels=[ElectronLiquidRG.Parquet.PPr, ElectronLiquidRG.Parquet.PHEr], transferLoop=transferLoop)  
        
        
    #     # ElectronLiquidRG.bLambda_MCMC(_para; kamp=Λgrid, n=[0, 0, 0], neval=neval, integrand=ElectronLiquidRG._PP, filename=filename_bPP)
    #     # ElectronLiquidRG.bLambda_MCMC(_para; kamp=Λgrid, n=[0,0,0], neval=neval, integrand=ElectronLiquidRG._Lver3_direct, filename=filename_bLver3_direct)
    #     # ElectronLiquidRG.bLambda_MCMC(_para; kamp=Λgrid, n=[0,0,0], neval=neval, integrand=ElectronLiquidRG._Lver3_exchange, filename=filename_bLver3_exchange)
    #     # ElectronLiquidRG.bLambda_MCMC(_para; kamp=Λgrid, n=[0,0,0], neval=neval, integrand=ElectronLiquidRG._Lver3_bubble, filename=filename_bLver3_bubble)
        
    #     # ElectronLiquidRG.sigma(_para; neval=neval, Λgrid=Λgrid, filename=filename_sigma)
    # GC.gc()
end