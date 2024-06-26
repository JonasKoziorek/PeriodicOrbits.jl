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
    igs = InitialGuess[InitialGuess([0.0, 0.0, 0.0], 0.0)]
    indss, signss = lambdaperms(dimension(ds))
    alg = DavidchackLai(n=3, m=2)
    periodic_orbits(pmap, alg, igs)
end