using PeriodicOrbits: lambdamatrix, lambdaperms
using PeriodicOrbits
using DynamicalSystems

begin
    ds = Systems.henon()
    traj, t = trajectory(ds, 100)
    igs = InitialGuess[InitialGuess(x, 0.0) for x in traj]
end

begin
    indss, signss = lambdaperms(dimension(ds))
    alg2 = SchmelcherDiakonos(o=5, indss=indss, signss=signss, λs=[0.1])
    alg3 = SchmelcherDiakonos(5, dimension(ds), 0.1)
    alg4 = SchmelcherDiakonos(5, [0.1], indss, signss)
    @time orbits2 = periodic_orbits(ds, alg2, igs)
    @time orbits3 = periodic_orbits(ds, alg3, igs)
    @time orbits4 = periodic_orbits(ds, alg4, igs)
end;