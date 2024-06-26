using PeriodicOrbits
using PeriodicOrbits: lambdaperms

thomas_cyclical(u0 = [1.0, 0, 0]; b = 0.2) = CoupledODEs(thomas_rule, u0, [b])

function thomas_rule(u, p, t)
    x,y,z = u
    b = p[1]
    xdot = sin(y) - b*x
    ydot = sin(z) - b*y
    zdot = sin(x) - b*z
    return SVector{3}(xdot, ydot, zdot)
end
function thomas_jacob(u, p, t)
    x,y,z = u
    b = p[1]
    return SMatrix{3,3}(-b, 0, cos(x), cos(y), -b, 0, 0, cos(z), -b)
end

begin
    ds = thomas_cyclical(b = 0.1665);
    plane = (3, 0.0)
    pmap = PoincareMap(ds, plane)
    xg = yg = range(-6.0, 6.0; length = 5)
    igs = InitialGuess[InitialGuess([x, y, z], 0.0) for x in xg for y in yg for z in [0.0]]
    indss, signss = lambdaperms(dimension(ds))
    alg = SchmelcherDiakonos(o=1, indss=indss, signss=signss, λs=[0.1])
    po = periodic_orbits(pmap, alg, igs)
    po2 = find_minimal_period(pmap, po, 1e-3)
    display(po.points)
    display(po2.points)
end

begin
    alg = DampedNewtonRaphsonMees(δ=2^(-1), J=thomas_jacob, maxiter=5000, disttol=1e-6)
    traj, t = trajectory(ds, 10; Dt=0.1)
    # igs = InitialGuess[InitialGuess(x, 5*rand()) for x in traj]
    step!(ds, 100*rand())
    igs = InitialGuess[InitialGuess(SVector{3}(current_state(ds)), 10*rand())]
    @time pos = periodic_orbits(ds, alg, igs)

    if length(pos.points) > 0
        reinit!(ds, pos.points[1].u)
        step!(ds, 1.0)
        if DynamicalSystemsBase.norm(pos.points[1].u - current_state(ds)) < 1e-3
            print("FP")
        else
            print("UPO")
        end
    end
end


reinit!(ds, pos.points[1].u)