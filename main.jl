using Plots, LinearAlgebra

mutable struct Particle
    m::Real
    x
    v
    F
    F_old
    Particle(m, x, v, F, F_old) = new(m, x, v, F, F_old)
end



function timeIntegration_basis(t::Real, dt::Real, t_end::Real, p, N::Int64)
    result_x = []
    result_y = []
    result_v = []
    compF_basis(p, N)
    while (t < t_end)
        t += dt
        compX_basis(p, N, dt)
        compF_basis(p, N)
        compV_basis(p, N, dt)
        compoutStatistic_basis(p, N, t)
        for i = 1:N
            result_x = push!(result_x, p[i].x[1])
            result_y = push!(result_y, p[i].x[2])
        end
    end
    return result_x, result_y, result_v
end

function compoutStatistic_basis(p, N::Int64, t::Real)
    e = 0
    for i=1:N
        v = 0
        for j = 1:2
            v += (p[i].v[j])^2
        end
        e += 0.5*p[i].m*v
    end
    return e
end

function compX_basis(p, N::Int64, dt::Real)
    for i = 1:N
        updateX(p[i], dt)
    end
end

function updateX(p, dt)
    a = dt*0.5/p.m
    p.x = p.x + dt*(p.v + a*p.F)
    p.F_old = p.F
end

function compV_basis(p, N::Int64, dt::Real)
    for i = 1:N
        updateV(p[i], dt)
    end
end

function updateV(p, dt)
    a = dt*0.5/p.m
    p.v = p.v + a*(p.F + p.F_old)
end

function compF_basis(p, N::Int64)
    for i = 1:N
        p[i].F = [0,0]
    end
    for i = 1:N
        for j = 1:N
            if i != j
                p[i].F = force(p[i], p[j])
            end
        end
    end
end

function force(p1, p2)
    r_vector = p2.x - p1.x
    rmag = norm(p1.x - p2.x)
    f = ((p1.m*p2.m)/rmag^3)*r_vector
    print(f)
    p1.F = p1.F + f
end

function main(particles, dt, t_end)
    N = length(particles)
    timeIntegration_basis(0, dt, t_end, particles, N)
end

#mass, x0, v0, F
p_sun = Particle(1, [0,0], [0,0], [0,0], [0,0])
p_earth = Particle(3e-6, [0, 1], [-1, 0], [0,0], [0,0])
p_jupiter = Particle(9.55e-4, [0,5.36], [-0.425, 0], [0,0], [0,0])
p_halley = Particle(1e-14, [34.75, 0], [0,0.0296], [0,0], [0,0])

x, y, v = main([p_sun, p_earth, p_jupiter, p_halley], .1, 95)
plot([0],[0], label="Sun", marker=:star, markersize=20)
plot!(x[2:4:end], y[2:4:end], label="Jupiter", lw=2)
plot!(x[3:4:end], y[3:4:end], label="Earth")
savefig("Orbit.pdf")
#plot!(x[4:3:end], y[4:3:end], label="Halley") ???

