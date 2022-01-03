# DESCRIPTION
# Computes the linearized dynamics of a nonlinear system (acfun), given the
# trim state and control vectors (x0 and u0).
#
# INPUT
# - acfun: name of system dynamics function
# - x0: trim state vector
# - u0: trim control vector
#
# OUTPUT
# - A, B, C, D: coefficient matrices of linearized dynamics
#
# ------------------------------------------------------------------------------

function linearize(acfun,x0,u0)
    #xdot0=H60!(xdot_in,x0,u0,t)
    #xdot0=acfun(xdot_in,x0,u0,t)
    A=zeros(NSTATES,NSTATES)
    B=zeros(NSTATES,NCTRLS)
    C=zeros(NOUT,NSTATES)
    D=zeros(NOUT,NCTRLS)
    # A and B matrices
    x_p=zeros(NSTATES)
    u_p=zeros(NCTRLS)
    for k=1:NSTATES
        x_p[:]=x0[:]
        x_p[k]=x_p[k]+DELXLIN[k]
        xdot_p1, y_p1=acfun(x_p,u0)
        x_p[k]=x_p[k]-2*DELXLIN[k]
        xdot_p2, y_p2=acfun(x_p,u0)
        #
        A[:,k]=(xdot_p1-xdot_p2)/(2*DELXLIN[k])
        C[:,k]=(y_p1-y_p2)/(2*DELXLIN[k])
    end
    # C and D matrices
    for k=1:NCTRLS
        u_p[:]=u0[:]
        u_p[k]=u_p[k]+DELCLIN[k]
        xdot_p1, y_p1=acfun(x0,u_p)
        u_p[k]=u_p[k]-2*DELCLIN[k]
        xdot_p2, y_p2=acfun(x0,u_p)
        B[:,k]=(xdot_p1-xdot_p2)/(2*DELCLIN[k])
        D[:,k]=(y_p1-y_p2)/(2*DELCLIN[k])
    end
    #
    return A, B, C, D
end
