# using CUDA
using StaticArrays
using Interpolations
using LinearAlgebra
using DifferentialEquations
using Cubature
using DelimitedFiles
# using Interpolations
# using Plots
using Dierckx
using Distributed
using IntelVectorMath
using LoopVectorization
using BenchmarkTools
# using GR
# import Pkg; Pkg.add("SymPy")
# using SymPy
# const sympy_parsing_mathematica =
#     SymPy.PyCall.pyimport("sympy.parsing.mathematica")
# mathematica2julia(
#     s::AbstractString,
#     substitutions::Pair{<:AbstractString,<:AbstractString}...,
# ) = SymPy.walk_expression(sympy_parsing_mathematica."mathematica"(
#     s,
#     Dict(substitutions...),
# ))
#
# Pkg.add("SymPy")
#
# Pkg.build("PyCall")
#
# ENV["PYTHON"]
#
# ENV["PYTHON"]="python3"
#
# using Conda
# rm(Conda.ROOTENV, recursive=true)
# Conda.add("numpy")
# Conda.runconda(`search numpy=1.15.4`)
# Pkg.add("Conda")
#
# rm(raw"C:\Users\tyy\.juliapro\JuliaPro_v1.4.2-1\conda", recursive=true)
# Pkg.build("PyCall")
# #
# #
# #
# #
# # Pkg.add("Conda")
# #
# # using Conda
# #
# # Conda.add("numpy")
