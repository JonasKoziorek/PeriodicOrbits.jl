export InitialGuess, PeriodicPoint, PeriodicOrbit, PeriodicOrbitFinder, find_minimal_period

struct InitialGuess{U<:AbstractArray{<:Real}, R<:Union{Real, Nothing}}
    u0::U
    T::R
end
# InitialGuess(u0) = InitialGuess(u0, nothing) #  what does it mean to have period equal to nothing?
# if you have several guesses:
# guesses = [InitialGuess(u, t) for (u, t) in zip(guesses)]

struct PeriodicPoint{U<:AbstractArray{<:Real}, R<:Real}
    u::U # point of the periodic orbit
    T::R # period
end

# struct PeriodicOrbit{U<:AbstractArray{<:Real}, R}
#     points::Vector{U}
#     T::R
# end

# -- change since we also want to store POs with different periods
struct PeriodicOrbit{U<:AbstractArray{<:Real}, R}
    points::Vector{PeriodicPoint{U, R}}
end

Base.show(io::IO, ::MIME"text/plain", po::PeriodicOrbit) = print(io, po.points)

function Base.push!(po::PeriodicOrbit, pp::PeriodicPoint)
    push!(po.points, pp)
end


function Base.:∈(u0::PeriodicPoint, POs::PeriodicOrbit)
    # custom search
    # (discrete) - linear search through the set
    # (continuous) - distinguish identical periodic orbits
end

function stable(ds, po::PeriodicPoint; jac=autodiff_jac(ds))::Bool
end

function find_minimal_period(ds::DynamicalSystem, pp::PeriodicPoint, disttol::Real)
    result :: PeriodicPoint # with minimal period
    return result
end

function find_minimal_period(ds::DiscreteTimeDynamicalSystem, pp::PeriodicPoint, disttol::Real=1e-6)
    period = pp.T
    for i in 1:period
        reinit!(ds, pp.u)
        step!(ds, i)
        if norm(pp.u - current_state(ds)) < disttol
            return PeriodicPoint(pp.u, i)
        end
    end
    return pp
end

function find_minimal_period(ds, po::PeriodicOrbit, disttol::Real=1e-6)
    return PeriodicOrbit(find_minimal_period.(ds, po.points, disttol))
end

function complete_orbit(ds, po::PeriodicPoint)
    # compute trajectory for period po.T
    result :: PeriodicOrbit
    return result
end

abstract type PeriodicOrbitFinder end

@kwdef struct Algorithm1 <: PeriodicOrbitFinder
    param1 = 1
    param2 = 2
    param3 = 3
end

function check_parameters(alg::A) where {A <: PeriodicOrbitFinder}
    # check if parameters are valid
end

function periodic_orbit(ds::DynamicalSystem, alg::PeriodicOrbitFinder, igs::Vector{InitialGuess} = [InitialGuess(ds)])
    result::PeriodicOrbit
    return result
end

function periodic_orbits(ds::DynamicalSystem, alg::PeriodicOrbitFinder, igs::Vector{InitialGuess} = [InitialGuess(ds)])
    result::PeriodicOrbit
    return result
end

periodic_orbits(ds::DynamicalSystem, alg::PeriodicOrbitFinder, igs::Vector{AbstractVector{Real}}) = periodic_orbits(ds, alg, InitialGuess[InitialGuess(u, 0.0) for u in igs])

# InitialGuess(ds::DynamicalSystem, T = nothing) = InitialGuess(current_state(ds), T) # what does it mean to have period equal to nothing?