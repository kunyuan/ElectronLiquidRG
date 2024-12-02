module ElectronLiquidRG

using JLD2
using CompositeGrids
using MCIntegration
using ElectronLiquid
using Measurements
using GreenFunc
using FeynmanDiagram
using FiniteDifferences
using ElectronGas

using LinearAlgebra
using Lehmann

include("para_tabel.jl")
include("R_derivative.jl")
# include("sigma.jl")
# include("vertex4.jl")
# include("vertex3.jl")

include("ver3_oneloop.jl")
include("sigma_oneloop.jl")
include("ver4_oneloop.jl")
end
