using Test
using MPCCases, Printf
using JuMP, JuMPUtil, Ipopt, MathProgBase
#using JuMP, Ipopt, MathProgBase
using SparseArrays, LinearAlgebra
include("../src/OPF.jl")

#const path = "/Users/jakeroth/Desktop/tdoptimizationuncertainty/CODE/OPF/test/cases/"
const path = "/Users/johannesbrust/Dropbox/ANL/projects/LOAD_DYNAMIC/code/OPF/test/cases"
const case = "case9"
# const case = "case5"
const tol = 1e-9
opfdata = load_case(case, path, other=false);

## deterministic acopf
dm = acopf_model(opfdata)
dm = acopf_solve(dm, opfdata)
dm_eval = setup(dm.m);
dxbar = copy(dm_eval.last_x);
acopf_outputAll(dm, opfdata)

## stochastic acopf
sm = sacopf_model(opfdata)
sm = acopf_solve(sm, opfdata)
sm_eval = setup(sm.m);
sxbar = copy(sm_eval.last_x);
acopf_outputAll(sm, opfdata)

@testset "indexing" begin
@test norm([ getvalue(dm.m[:Pg]);
             getvalue(dm.m[:Qg]);
             getvalue(dm.m[:Vm]);
             getvalue(dm.m[:Va]) ] - dxbar) <= tol
@test norm([ getvalue(sm.m[:Pg]);
             getvalue(sm.m[:Qg]);
             getvalue(sm.m[:Vm]);
             getvalue(sm.m[:Va]);
             getvalue(sm.m[:Pd]);
             getvalue(sm.m[:Qd]) ] - sxbar) <= tol
@test norm([ getvalue(getindex(dm.m, :Pg));
             getvalue(getindex(dm.m, :Qg));
             getvalue(getindex(dm.m, :Vm));
             getvalue(getindex(dm.m, :Va)) ] - dxbar) <= tol
@test norm([ getvalue(getindex(sm.m, :Pg));
             getvalue(getindex(sm.m, :Qg));
             getvalue(getindex(sm.m, :Vm));
             getvalue(getindex(sm.m, :Va));
             getvalue(getindex(sm.m, :Pd));
             getvalue(getindex(sm.m, :Qd)) ] - sxbar) <= tol
end

dFdy, dFdx = dFdy_dFdx(sm_eval, sxbar, opfdata)
dydx = Matrix(dFdy) \ Matrix(dFdx)
