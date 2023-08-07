
fdict = Dict()
fdict[1.0] = collect(LinRange(-1.0, 0.0, 16))

Λgrid(kF) = CompositeGrid.LogDensedGrid(:gauss, [1.0 * kF, 100 * kF], [kF,], 8, 0.01 * kF, 8)

get_para(para, Fs) = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=Fs, Fa=-0.0, order=para.order,
    mass2=para.mass2, isDynamic=true, isFock=false)