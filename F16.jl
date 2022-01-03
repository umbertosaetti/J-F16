# DESCRIPTION
# Nonlinear simulation model of the rigid-body dynamics of an F-16
# aircraft.
#
# INPUT
# - x: state vector
# - p: control input vector
#
# OUTPUT
# - XD: state derivative
# - Y: output vector 
#
# ------------------------------------------------------------------------------

function F16!(x,p)
    # unpack state
    VT=x[1]
    ALPHA=x[2]*R2D
    BETA=x[3]*R2D
    PHI=x[4]
    THETA=x[5]
    PSI=x[6]
    P=x[7]
    Q=x[8]
    R=x[9]
    NORTH=x[10]
    EAST=x[11]
    ALT=x[12]
    if ALT<0.0
        ALT=0.0
    end
    POW=x[13]
    # unpack controls
    THTL=p[1]
    AIL=p[2]
    EL=p[3]
    RDR=p[4]
    # define state derivative vector
    XD=zeros(13,1)
    # convert attitude and angular rates from [rad] to [deg]
    PHID = PHI*R2D
    THETAD = THETA*R2D
    PD = P*R2D
    QD = Q*R2D
    RD = R*R2D
    # air data computer and engine model
    MACH, QBAR = ADC(VT,ALT)
    CPOW = TGEAR(THTL)
    XD[13] = PDOT(POW,CPOW)
    T = THRUST(POW,ALT,MACH)
    # component buildup
    CXT = CX(ALPHA,EL)
    CYT = CY(BETA,AIL,RDR)
    CZT = CZ(ALPHA,BETA,EL)
    #
    DAIL= AIL/20.0
    DRDR= RDR/30.0
    #
    CLT = CL(ALPHA,BETA) + DLDA(ALPHA,BETA)*DAIL + DLDR(ALPHA,BETA)*DRDR
    CMT = CM(ALPHA,EL)
    CNT = CN(ALPHA,BETA) + DNDA(ALPHA,BETA)*DAIL + DNDR(ALPHA,BETA)*DRDR
    # add damping derivatives
    TVT = 0.5/VT
    B2V = B*TVT
    CQ = CBAR*Q*TVT
    D = DAMP(ALPHA)
    #
    CXT = CXT + CQ*D[1]
    CYT = CYT + B2V*(D[2]*R+D[3]*P)
    CZT = CZT + CQ*D[4]
    #
    CLT = CLT + B2V*(D[5]*R+D[6]*P)
    CMT = CMT + CQ*D[7] + CZT*(XCGR-XCG)
    CNT = CNT + B2V*(D[8]*R+D[9]*P) - CYT*(XCGR-XCG)*CBAR/B
    # prepare state equations
    CBTA  = cos(x[3])
    U = VT*cos(x[2])*CBTA
    V = VT*sin(x[3])
    W = VT*sin(x[2])*CBTA
    #
    STH = sin(THETA)
    CTH = cos(THETA)
    SPH = sin(PHI)
    CPH = cos(PHI)
    SPSI = sin(PSI)
    CPSI = cos(PSI)
    #
    RM = G/WEIGHT
    QS = QBAR*S
    QSB = QS*B
    RMQS = RM*QS
    GCTH = G*CTH
    QSPH = Q*SPH
    AX = RM*(QS*CXT+T)
    AY = RMQS*CYT
    AZ = RMQS*CZT
    # equations of motion
    # longitudinal dynamics
    UDOT = R*V - Q*W - G*STH + AX
    VDOT = P*W - R*U + GCTH*SPH + AY
    WDOT = Q*U - P*V + GCTH*CPH + AZ
    DUM = (U*U + W*W)
    XD[1] = (U*UDOT+V*VDOT+W*WDOT)/VT
    XD[2] = (U*WDOT-W*UDOT)/DUM
    XD[3] = (VT*VDOT-V*XD[1])*CBTA/DUM
    # rotational kinematics
    XD[4] = P + (STH/CTH)*(QSPH + R*CPH)
    XD[5] = Q*CPH - R*SPH
    XD[6] = (QSPH + R*CPH)/CTH
    # rotational dynamics
    XD[7] = (C2*P + C1*R + C4*HE)*Q + QSB*(C3*CLT + C4*CNT)
    XD[8] = (C5*P - C7*HE)*R + C6*(R*R-P*P) +QS*CBAR*C7*CMT
    XD[9] = (C8*P-C2*R+C9*HE)*Q + QSB*(C4*CLT + C9*CNT)
    # navigation
    T1 = SPH*CPSI
    T2 = CPH*STH
    T3 = SPH*SPSI
    S1 = CTH*CPSI
    S2 = CTH*SPSI
    S3 = T1*STH-CPH*SPSI
    S4 = T3*STH+CPH*CPSI
    S5 = SPH*CTH
    S6 = T2*CPSI+T3
    S7 = T2*SPSI-T1
    S8 = CPH*CTH
    XD[10] = U*S1 +V*S3 + W*S6
    XD[11] = U*S2 + V*S4 + W*S7
    XD[12] = -(-U*STH + V*S5 + W*S8)
    # outputs
    # measured accelerometer at pilot location
    AYN = AY/G
    AXN = AX/G
    AZN = (-AZ+XA*XD[8])/G
    # cstar parameter
    CSTAR = AZN + 12.4*Q
    Y = [AZN; AYN; AXN; QBAR; MACH; CSTAR]
    # return state derivative and vectors
    return XD, Y
end
