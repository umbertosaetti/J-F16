# DESCRIPTION
# Runge-Kutta 4 integration scheme.
#
# INPUT
# - fdyn: name of system dynamics function
# - finp: name of control input function
# - t: time
# - dt: time step
# - x_curr: current state vector
# - u0: initial control vector
#
# OUTPUT
# - x_new: new state vector after integration
#
# ------------------------------------------------------------------------------

function rk4(fdyn,finp,t,dt,x_curr,u0)
    #
    u=finp(t)+u0
    xd, y = fdyn(x_curr,u)
    k1=dt*xd
    #
    u=finp(t+0.5*dt)+u0
    xd, y = fdyn(x_curr+0.5*k1,u)
    k2=dt*xd
    #
    xd, y = fdyn(x_curr+0.5*k2,u)
    k3=dt*xd
    #
    u=finp(t+dt)+u0
    xd, y = fdyn(x_curr+k3,u)
    k4=dt*xd
    # weighted average
    x_new=x_curr+k1/6+k2/3+k3/3+k4/6
    #
    return x_new
end
