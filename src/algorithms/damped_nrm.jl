export periodic_orbits, DampedNewtonRaphsonMees

@kwdef struct DampedNewtonRaphsonMees <: PeriodicOrbitFinder
    δ = 2^(-3)
    J = nothing
    maxiter = 10000
    disttol = 1e-6
end

function periodic_orbits(ds::CoupledODEs, alg::PeriodicOrbitFinder, igs::Vector{InitialGuess})
    tands = TangentDynamicalSystem(ds; J = alg.J, k=dimension(ds))
    pos = PeriodicOrbit{typeof(current_state(ds)), typeof(current_time(ds))}([])

    j = 0
    for ig in igs
        j+=1
        println(j)
        # try
            prev_step = SVector{length(ig.u0)+1}(vcat(ig.u0, ig.T))

            i = 1
            reinit!(tands, prev_step[1:end-1])
            step!(tands, prev_step[end])
            # while DynamicalSystemsBase.norm(current_state(tands) - prev_step[1:end-1]) > alg.disttol && i < alg.maxiter
            while true
                if DynamicalSystemsBase.norm(current_state(tands) - prev_step[1:end-1]) < alg.disttol || i > alg.maxiter
                    push!(pos, PeriodicPoint{typeof(current_state(ds)), typeof(current_time(ds))}(prev_step[1:end-1], prev_step[end]))
                    break
                end
                # println(i)
                prev_step = next_step(prev_step, tands, alg.δ)
                i+=1
            end
            # println(i)
            # println(prev_step)
        # catch e
            # @warn e.msg
        # end
    end
    return pos
end

function next_step(prev_step, tands, δ)
    X = prev_step[1:end-1]
    T = prev_step[end]
    reinit!(tands, X)
    step!(tands, T)
    Phi = current_deviations(tands)
    phi = current_state(tands)

    F = tands.original_f
    p = current_parameters(tands)
    t = current_time(tands)

    I = DynamicalSystemsBase.I
    A = hcat(vcat(Phi-I, F(X, p, t)'), vcat(F(phi, p, t), 0))
    b = vcat(X-phi, 0)

    Δ = A\b
    next_step = prev_step + δ*Δ

    return next_step
end