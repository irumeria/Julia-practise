using GSL

include("constants.jl")

struct Orbital
    n
    l
    m
    rho_over_r::Float64
    root_term::Float64
end

function newOrbital(n, l, m)
    root_term = sqrt(
        8.0 * factorial(n-l-1) /
        (
            (n * n * n * n * 2.0 * factorial(n+l))
            * REDUCED_BOHR_RADIUS
            * REDUCED_BOHR_RADIUS
            * REDUCED_BOHR_RADIUS
        )
    )
    rho_over_r = 2.0 / (n * REDUCED_BOHR_RADIUS)

    Orbital(
        n,
        l,
        m,
        rho_over_r,
        root_term,
    )
end

function probability_with_phase(
    orbit::Orbital, 
    coord::Vector,  # x y z
    volume::Float64)

    (ψ, phase) = ψ_with_phase(orbit, coord)

    prob = get_sqr_complex(ψ) * volume

    (prob, phase)
end

function get_sqr_complex(complex_number)
    real_part = real(complex_number)
    img_part = real(im * complex_number)
    real_part^2 + img_part^2
end

function ψ_with_phase(orbit::Orbital, coord::Vector)
    rho = radius(coord) * orbit.rho_over_r
    radial = orbit.root_term * exp((-rho / 2.0)) *
                rho^orbit.l * 
                GSL.sf_laguerre_n(
                        orbit.n - orbit.l - 1, 
                        2 * orbit.l + 1, 
                        rho
                    )

    sph_harmonics = harmonics(orbit.l, orbit.m, coord)

    ψ = radial * sph_harmonics

    phase = getPhase(orbit, sph_harmonics, radial, coord)

    (ψ, phase)
end


function getPhase(orbit, sph_harmonics::ComplexF64, radial::Float64, coord::Vector)
    # // Phase calculation
    # // This is the sign of R(r)Y_(m, l)(theta, phi), but Y_(m, l) is in its real form
    # // (since we can't take the sign of a the complex number psi)

    # // We use the Condon-Shortley phase convention for our definition of Y
    condon_shortley_sign = abs(orbit.m) % 2 == 0 ? 1. : -1. 
    r_sph_harm = 0
    if orbit.m < 0 
        # // Since we need to calculate Y_(-m,l) anyway, it's cheaper to calculate the real value directly
        r_sph_harm = condon_shortley_sign * real(harmonics(orbit.l, orbit.m, coord))
    elseif orbit.m == 0 
        r_sph_harm = real(sph_harmonics)
    else 
        r_sph_harm = condon_shortley_sign * SQRT_2 * real(sph_harmonics)
    end
    res = 0
    if radial * r_sph_harm > 0.0 
        res = 1
    elseif radial * r_sph_harm == 0.0 
        res = 0
    else 
        res = -1
    end
    res
end

function radius(coord)
    return sqrt(coord[1]^2 + coord[2]^2 + coord[3]^2)
end

function theta(coord)
    acos(coord[3] / radius(coord))
end

function phi(coord)
    atan(coord[2] / coord[1])
end

function harmonics(l, m, coord::Vector)
    (-1)^(m) *
    GSL.sf_legendre_sphPlm(l, m, cos(theta(coord))) *
    exp(im * m * phi(coord))
end

# orbit = newOrbital(2, 1, 0)
# r_bound = 0.000000001175595458901099
# grid_size = 1024
# volume = (r_bound / (grid_size / 2.0 - 1.0))^3

# coord = [r_bound/4,0.0,r_bound/8]
# rho = radius(coord) * orbit.rho_over_r
# radial = orbit.root_term * exp((-rho / 2.0)) *
#          rho^orbit.l * GSL.sf_laguerre_n(orbit.n - orbit.l - 1, 2 * orbit.l + 1, rho)

# sph_harmonics = harmonics(orbit.l, orbit.m, coord)

# ψ = radial * sph_harmonics
# prob = get_sqr_complex(ψ) * volume

# using Polynomials
# print(
#     # harmonics(1, 0, [0.00000001, 0.00000001, 1.0])
#     # GSL.sf_laguerre_n(1,3,rho)
#     # GSL.sf_laguerre_n(1, 3, rho)
#     # convert(Polynomial, orbit.laguerre)
#     # Laguerre{3}(3)(rho)
#     sph_harmonics
# )

# @time harmonics(1, 1, [11.66, -35, -35])
