################## PCD (Probability Constrained Dispatch) ################
# Formulation of probability constrained optimal power flow (OPF) problems.
#
# test_jac_mpb.jl is an initial script to compute constraint Jacobians using
# the MathProgBase routines.
###########################################################################
# 02/04/19, J.B.

# Declarations from 'runtests.jl'
using Test
using MPCCases, Printf
using JuMP, JuMPUtil, Ipopt, MathProgBase
using SparseArrays, LinearAlgebra
include("../src/OPF.jl")

const path = "/Users/johannesbrust/Dropbox/ANL/projects/LOAD_DYNAMIC/code/OPF/test/cases"
const case = "case9"

opfdata = load_case(case, path, other=false);

# Compute solution for stochastic problem formulation
sm = sacopf_model(opfdata)
sm = acopf_solve(sm, opfdata)

# Problem parameters
nbus    = length(opfdata.buses)
ngen    = length(opfdata.generators)
nvars   = 2*ngen+4*nbus
ncons   = 2*nbus

# Jump model and computation of Jacobian
# Uses MathProgBase
jmdl    = sm.m
d       = JuMP.NLPEvaluator(jmdl)

MathProgBase.initialize(d,[:Jac])

# Variables
# Pg, Qg, Vm, Va, Pd, Qd
xopt              = zeros(nvars)
PG                = getvalue(getindex(jmdl,:Pg))
QG                = getvalue(getindex(jmdl,:Qg))
VM                = getvalue(getindex(jmdl,:Vm))
VA                = getvalue(getindex(jmdl,:Va))
PD                = getvalue(getindex(jmdl,:Pd))
QD                = getvalue(getindex(jmdl,:Qd))

idxs              = 1
idxe              = ngen
xopt[idxs:idxe]   = PG; idxs = idxe + 1; idxe = idxs + ngen - 1
xopt[idxs:idxe]   = QG; idxs = idxe + 1; idxe = idxs + nbus - 1
xopt[idxs:idxe]   = VM; idxs = idxe + 1; idxe = idxs + nbus - 1
xopt[idxs:idxe]   = VA; idxs = idxe + 1; idxe = idxs + nbus - 1
xopt[idxs:idxe]   = PD; idxs = idxe + 1; idxe = idxs + nbus - 1
xopt[idxs:idxe]   = QD

# Jacobian computations
#full_jac          = zeros(2*ncons,nvars)
(Idxi,Idxj)       = MathProgBase.jac_structure(d)
jac_v             = zeros(length(Idxi))

MathProgBase.eval_jac_g(d,jac_v,xopt)

sjac              = sparse(Idxi,Idxj,jac_v,2*ncons,nvars) # sparse Jacobian
jac               = Array(sjac) # dense Jacobian

# Computing indices, and selecting columns and rows of
# the Jacobian.
# Columns PG      -- VA(L+G) define dFdy,
# Columns PD(L+G) -- QD(L+G) define dFdx
#
#         PG  QG(ref) QG VM(G) VM(L) VA(ref) VA(L+G) PD(ref) PD(L+G) QD(ref) QD(L+G)
# refIdx
# PIdx                X         X               X               X               X
# refIdx
# QIdx                X         X               X               X               X
#
refIdx            = opfdata.bus_ref # Reference bus. Assume at 1
loadIdx           = 1:nbus
loadIdx           = loadIdx[length.(opfdata.BusGenerators[1:nbus]).==0]

QGIdx             = (ngen + refIdx + 1):2*ngen
VMIdx             = QGIdx[end] .+ loadIdx
VAIdx             = (VMIdx[end] + refIdx + 1):VMIdx[end] + nbus
PDIdx             = (VAIdx[end] + refIdx + 1):VAIdx[end] + nbus
QDIdx             = (PDIdx[end] + refIdx + 1):PDIdx[end] + nbus

yIdx              = [QGIdx;VMIdx;VAIdx]
xIdx              = [PDIdx;QDIdx]

rowIdx_           = (refIdx+1):nbus
rowIdx            = [rowIdx_; (rowIdx_[end] + refIdx + 1):2*nbus]

dFdy              = jac[rowIdx,yIdx]
#dFdx              = jac[rowIdx,xIdx]
