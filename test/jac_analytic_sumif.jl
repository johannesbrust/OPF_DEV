################## PCD (Probability Constrained Dispatch) ################
# Formulation of probability constrained optimal power flow (OPF) problems.
#
# jac_analytic_sumif.jl is an initial implementation to compute constraint Jacobians using
# derivatives of the power flow equations.
#
# INPUTS:
# opf_data  := Structure with raw problem information
# (Pg)      := Array, ngx1 (Power generation)
# (Qg)      := Array, ngx1 (Imaginary power generation)
# Vm        := Array, nbx1 (Voltages)
# Va        := Array, nbx1 (Angles)
# Pd        := Array, nbx1 (Loads)
# Qd        := Array, nbx1 (Imaginary loads)
#
#
# OUTPUTS:
# jac       := Array, 2nbx(2ng+4nb)
#
# This function tests the use of 'if' conditions in computing the Jacobian in
# order to use 'sum if' constructions
###########################################################################
# 04/19/19, J.B.

function jac_analytic_sumif(opf_data,Pg,Qg,Vm,Va,Pd,Qd)

    # Initializations
    lines       = opf_data.lines;
    buses       = opf_data.buses;
    generators  = opf_data.generators;
    baseMVA     = opf_data.baseMVA;
    busIdx      = opf_data.BusIdx;
    FromLines   = opf_data.FromLines;
    ToLines     = opf_data.ToLines;
    BusGeners   = opf_data.BusGenerators;

    nbus        = length(buses);
    nline       = length(lines);
    ngen        = length(generators);

    jac         = zeros(2*nbus,2*ngen+4*nbus);

    YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI = computeAdmitances(lines, buses, baseMVA)

    # Derivatives w.r.t. Vm, Va. Jacobian formed rowwise jac[b,j]

    idxoffVm    = 2*ngen; # Offset Vm
    idxoffVa    = idxoffVm + nbus; # Offset Va

     for b in 1:nbus
         for j in 1:nbus

             # fl     = FromLines[b];
             # tl     = ToLines[b];
             # nfl    = length(fl);
             # ntl    = length(tl);

             row    = nbus + b;
             colm   = idxoffVm + j;
             cola   = idxoffVa + j;

            # Derivatives Vm, Va
            # Function modifies the location of 'if' conditions

            if (j == b)

                jac[b,colm]      = 2*(YshR[b])*Vm[b]; # Real
                jac[row,colm]    = 2*(-YshI[b])*Vm[b]; # Imaginary

            end

            # From lines
            #for l in 1:nfl
            for l in FromLines[b]

                if (j==b)

                    # Vm derivatives
                    jac[b,colm]      += 2*(YffR[l]*Vm[b]);
                    # + sum(YttR[tl[:]]) + YshR[b])*Vm[b]; # Real
                    # Vm derivatives
                    jac[row,colm]    += 2*(-YffI[l]*Vm[b]); # Imaginary

                    # Vm derivatives
                    jac[b,colm]      += Vm[busIdx[lines[l].to]] *
                     ( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]] ) +
                        YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )); # Real

                    jac[row,colm]    += Vm[busIdx[lines[l].to]] *
                         ( -YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]] ) +
                            YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )); # Imaginary

                    # Va derivatives
                    jac[b,cola]      += Vm[b]*Vm[busIdx[lines[l].to]] *
                     ( -YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]] ) +
                        YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  )); # Real

                    jac[row,cola]    += Vm[b]*Vm[busIdx[lines[l].to]] *
                         ( YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]] ) +
                            YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  )); # Imaginary

                elseif (j == busIdx[lines[l].to]) # Condition for off diagonal elements

                    # Vm derivatives
                    jac[b,colm]  = Vm[b] *
                     ( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]] ) +
                        YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )); # Real

                    jac[row,colm]  = Vm[b] *
                     ( -YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]] ) +
                        YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )); # Imaginary

                    # Va derivatives
                    jac[b,cola]  = Vm[b] * Vm[busIdx[lines[l].to]]*
                         ( YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]] ) +
                            -YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  )); # Real

                    jac[row,cola]  = Vm[b] * Vm[busIdx[lines[l].to]]*
                        ( -YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]] ) +
                            -YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  )); # Imaginary
                end

            end
            # To lines
            for l in ToLines[b]

                if (j == b)

                    jac[b,colm]      += 2*(YttR[l])*Vm[b]; # Real
                    jac[row,colm]    += 2*(-YttI[l])*Vm[b]; # Imaginary

                    # Vm derivatives
                    jac[b,colm] += Vm[busIdx[lines[l].from]] *
                     ( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]] ) +
                        YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]]  )); # Real

                    jac[row,colm] += Vm[busIdx[lines[l].from]] *
                     ( -YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]] ) +
                        YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]]  )); # Imaginary

                    # Va derivatives
                    jac[b,cola] += Vm[b]*Vm[busIdx[lines[l].from]] *
                         ( -YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]] ) +
                            YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]  )); # Real

                    jac[row,cola] += Vm[b]*Vm[busIdx[lines[l].from]] *
                         ( YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]] ) +
                            YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]  )); # Imaginary

                elseif (j == busIdx[lines[l].from])

                    # Vm derivatives
                    jac[b,colm]  = Vm[b] *
                     ( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]] ) +
                        YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]]  )); # Real

                    jac[row,colm]  = Vm[b] *
                     ( -YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]] ) +
                        YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]]  )); # Imaginary

                    # Va derivatives
                    jac[b,cola]  = Vm[b] * Vm[busIdx[lines[l].from]] *
                     ( YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]] ) +
                        -YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]  )); # Real

                    jac[row,cola]  = Vm[b] * Vm[busIdx[lines[l].from]] *
                     ( -YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]] ) +
                        -YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]  )); # Imaginary

                end
            end
         end
     end

     return jac
end
