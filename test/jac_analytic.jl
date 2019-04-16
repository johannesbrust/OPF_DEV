################## PCD (Probability Constrained Dispatch) ################
# Formulation of probability constrained optimal power flow (OPF) problems.
#
# jac_analytic.jl is an initial implementation to compute constraint Jacobians using
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
# NOTE: This function assumes one reference bus at index 1, at most one generator
# at a bus, or one load at a bus.
###########################################################################
# 12/04/19, J.B.

function jac_analytic(opf_data,Pg,Qg,Vm,Va,Pd,Qd)

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

             fl     = FromLines[b];
             tl     = ToLines[b];
             nfl    = length(fl);
             ntl    = length(tl);

             row    = nbus + b;
             colm   = idxoffVm + j;
             cola   = idxoffVa + j;

                # Derivatives Vm, Va
                if (j == b)

                    jac[b,colm]      = 2*(sum(YffR[fl[:]]) + sum(YttR[tl[:]]) + YshR[b])*Vm[b]; # Real
                    jac[row,colm]    = 2*(sum(-YffI[fl[:]]) + sum(-YttI[tl[:]]) - YshI[b])*Vm[b]; # Imaginary

                    # From lines
                    for l in 1:nfl

                        # Vm derivatives
                        jac[b,colm]      += Vm[busIdx[lines[fl[l]].to]] *
                         ( YftR[fl[l]]*cos(Va[b]-Va[busIdx[lines[fl[l]].to]] ) +
                            YftI[fl[l]]*sin(Va[b]-Va[busIdx[lines[fl[l]].to]]  )); # Real

                        jac[row,colm]    += Vm[busIdx[lines[fl[l]].to]] *
                             ( -YftI[fl[l]]*cos(Va[b]-Va[busIdx[lines[fl[l]].to]] ) +
                                YftR[fl[l]]*sin(Va[b]-Va[busIdx[lines[fl[l]].to]]  )); # Imaginary

                        # Va derivatives
                        jac[b,cola]      += Vm[b]*Vm[busIdx[lines[fl[l]].to]] *
                         ( -YftR[fl[l]]*sin(Va[b]-Va[busIdx[lines[fl[l]].to]] ) +
                            YftI[fl[l]]*cos(Va[b]-Va[busIdx[lines[fl[l]].to]]  )); # Real

                        jac[row,cola]    += Vm[b]*Vm[busIdx[lines[fl[l]].to]] *
                             ( YftI[fl[l]]*sin(Va[b]-Va[busIdx[lines[fl[l]].to]] ) +
                                YftR[fl[l]]*cos(Va[b]-Va[busIdx[lines[fl[l]].to]]  )); # Imaginary

                    end
                    # To lines
                    for l in 1:ntl

                        # Vm derivatives
                        jac[b,colm] += Vm[busIdx[lines[tl[l]].from]] *
                         ( YtfR[tl[l]]*cos(Va[b]-Va[busIdx[lines[tl[l]].from]] ) +
                            YtfI[tl[l]]*sin(Va[b]-Va[busIdx[lines[tl[l]].from]]  )); # Real

                        jac[row,colm] += Vm[busIdx[lines[tl[l]].from]] *
                         ( -YtfI[tl[l]]*cos(Va[b]-Va[busIdx[lines[tl[l]].from]] ) +
                            YtfR[tl[l]]*sin(Va[b]-Va[busIdx[lines[tl[l]].from]]  )); # Imaginary

                        # Va derivatives
                        jac[b,cola] += Vm[b]*Vm[busIdx[lines[tl[l]].from]] *
                             ( -YtfR[tl[l]]*sin(Va[b]-Va[busIdx[lines[tl[l]].from]] ) +
                                YtfI[tl[l]]*cos(Va[b]-Va[busIdx[lines[tl[l]].from]]  )); # Real

                        jac[row,cola] += Vm[b]*Vm[busIdx[lines[tl[l]].from]] *
                             ( YtfI[tl[l]]*sin(Va[b]-Va[busIdx[lines[tl[l]].from]] ) +
                                YtfR[tl[l]]*cos(Va[b]-Va[busIdx[lines[tl[l]].from]]  )); # Imaginary

                    end

                else


                    # Check if variables are in from busses or to busses

                    hasfrom  = 0;
                    hasto    = 0;
                    l        = 0;

                    for ll in 1:nfl

                        if j == busIdx[lines[fl[ll]].to]

                            l       = ll;
                            hasfrom = 1;
                            break;

                        end
                    end

                    for ll in 1:ntl

                        if j == busIdx[lines[tl[ll]].from]

                            l       = ll;
                            hasto   = 1;
                            break;

                        end
                    end




                    #if (sum(j .== fl[:])==1)

                    #l           = j;

                    if hasfrom == 1

                        #l = j;

                        # Vm derivatives
                        jac[b,colm]  = Vm[b] *
                         ( YftR[fl[l]]*cos(Va[b]-Va[busIdx[lines[fl[l]].to]] ) +
                            YftI[fl[l]]*sin(Va[b]-Va[busIdx[lines[fl[l]].to]]  )); # Real

                        jac[row,colm]  = Vm[b] *
                         ( -YftI[fl[l]]*cos(Va[b]-Va[busIdx[lines[fl[l]].to]] ) +
                            YftR[fl[l]]*sin(Va[b]-Va[busIdx[lines[fl[l]].to]]  )); # Imaginary

                        # Va derivatives
                        jac[b,cola]  = Vm[b] * Vm[busIdx[lines[fl[l]].to]]*
                             ( YftR[fl[l]]*sin(Va[b]-Va[busIdx[lines[fl[l]].to]] ) +
                                -YftI[fl[l]]*cos(Va[b]-Va[busIdx[lines[fl[l]].to]]  )); # Real

                        jac[row,cola]  = Vm[b] * Vm[busIdx[lines[fl[l]].to]]*
                            ( -YftI[fl[l]]*sin(Va[b]-Va[busIdx[lines[fl[l]].to]] ) +
                                -YftR[fl[l]]*cos(Va[b]-Va[busIdx[lines[fl[l]].to]]  )); # Imaginary

                    end

                #elseif (sum(j .== tl[:])==1)

                #    l           = j;

                    if hasto == 1

                        #l = j;

                        # Vm derivatives
                        jac[b,colm]  = Vm[b] *
                         ( YtfR[tl[l]]*cos(Va[b]-Va[busIdx[lines[tl[l]].from]] ) +
                            YtfI[tl[l]]*sin(Va[b]-Va[busIdx[lines[tl[l]].from]]  )); # Real

                        jac[row,colm]  = Vm[b] *
                         ( -YtfI[tl[l]]*cos(Va[b]-Va[busIdx[lines[tl[l]].from]] ) +
                            YtfR[tl[l]]*sin(Va[b]-Va[busIdx[lines[tl[l]].from]]  )); # Imaginary

                        # Va derivatives
                        jac[b,cola]  = Vm[b] * Vm[busIdx[lines[tl[l]].from]] *
                         ( YtfR[tl[l]]*sin(Va[b]-Va[busIdx[lines[tl[l]].from]] ) +
                            -YtfI[tl[l]]*cos(Va[b]-Va[busIdx[lines[tl[l]].from]]  )); # Real

                        jac[row,cola]  = Vm[b] * Vm[busIdx[lines[tl[l]].from]] *
                         ( -YtfI[tl[l]]*sin(Va[b]-Va[busIdx[lines[tl[l]].from]] ) +
                            -YtfR[tl[l]]*cos(Va[b]-Va[busIdx[lines[tl[l]].from]]  )); # Imaginary

                    end
                end
         end
     end

     return jac
end
