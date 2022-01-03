# DESCRIPTION
# Closed-loop aircraft dynamics.
#
# INPUT
# - x: current state vector
# - u_curr: current control vector
#
# OUTPUT
# - xdot: system dynamics
# - out: output vector 
#
# ------------------------------------------------------------------------------

function F16_DI!(x,u_curr)
    # unpack states
    VT = x[1]
    phi = x[4]
    p = x[7]
    q = x[8]
    r = x[9]
    # unpack command filter states
    p_cmd = x[NSTATES+1]
    intperr_cmd = x[NSTATES+2]
    q_cmd = x[NSTATES+3]
    intqerr_cmd = x[NSTATES+4]
    r_cmd = x[NSTATES+5]
    intrerr_cmd = x[NSTATES+6]
    δ_a = x[NSTATES+7]
    δ_e = x[NSTATES+8]
    δ_r = x[NSTATES+9]
    # control perturbations
    Δδ_a = u_curr[2]-δ_a
    Δδ_e = u_curr[3]-δ_e
    Δδ_r = u_curr[4]-δ_r
    # feedforward
    pdot_cmd = -1/τ_p*p_cmd + k_a/τ_p*Δδ_a
    qdot_cmd = -1/τ_q*q_cmd + k_e/τ_q*Δδ_e
    rdot_cmd = -1/τ_r*r_cmd + k_r/τ_r*Δδ_r
    # feedback compensation
    νp = pdot_cmd + kp_roll*(p_cmd-p) + ki_roll*(intperr_cmd)
    νq = qdot_cmd + kp_pitch*(q_cmd-q) + ki_pitch*(intqerr_cmd)
    νr = rdot_cmd + kp_yaw*(r_cmd-r) + ki_yaw*(intrerr_cmd)
    # pseudo control
    ν=[νp, νq, νr]
    # total velocity [kts]
    V=VT/1.688
    # coefficient matrix interpolation
    idx1=Int(min(max(floor(V/100)+1,1),7))
    idx2=idx1+1;
    intscale=min((V-(idx1-1)*100)/100,1.)
    CA=CATAB[:,:,idx1]+intscale*(CATAB[:,:,idx2]-CATAB[:,:,idx1])
    CBinv=CBinvTAB[:,:,idx1]+intscale*(CBinvTAB[:,:,idx2]-CBinvTAB[:,:,idx1])
    # dynamic inversion
    Δu=CBinv*(ν-CA*[p, q, r])
    # control inputs
    u=[u_curr[1]; Δu[1]+δ_a; Δu[2]+δ_e; Δu[3]+δ_r]
    # state derivative
    der, out = F16!(x[1:NSTATES],u)
    # rigid-body state derivative
    xdot=zeros(NSTATES+NDISTATES,1)
    xdot[1:NSTATES]=der
    # command filters state derivative
    xdot[NSTATES+1] = pdot_cmd
    xdot[NSTATES+2] = p_cmd-p
    xdot[NSTATES+3] = qdot_cmd+r*sin(phi)
    xdot[NSTATES+4] = q_cmd-q
    xdot[NSTATES+5] = rdot_cmd+G/VT*phi
    xdot[NSTATES+6] = r_cmd-r
    xdot[NSTATES+7] = (u_curr[2]-δ_a)/τ_u
    xdot[NSTATES+8] = (u_curr[3]-δ_e)/τ_u
    xdot[NSTATES+9] = (u_curr[4]-δ_r)/τ_u
    #
    return xdot, out
end
