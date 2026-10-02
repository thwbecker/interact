# bp1_pt.jl
#
# Explicit (accelerated pseudo-transient) quasi-dynamic SEAS BP1-type solver,
# written with ParallelStencil in the style of the JustRelax miniapps.
#
# STATUS: this file is a translation of the tested NumPy prototype bp1_pt.py.
# It has not been run (no Julia available where it was written). Revised 2026-09-29:
# Data.Array instead of PTArray, index arrays instead of Bool masks, explicit
# kernel ranges, check_every and wall_limit keywords. Expect small
# syntax or indexing errors; the algorithm and all constants match bp1_pt.py,
# which passed a manufactured-solution convergence test and a fault-solve
# consistency test. Compare the first steps against run80.csv from the Python
# version before trusting anything.
#
# Discretization and algorithm: see the docstring of bp1_pt.py.
# Arrays: u[i, j], i = 1..nx+1 along x (fault at i = 1), j = 1..nz+1 along z
# (free surface at j = 1). Fault RSF nodes: j with z[j] <= Wf.
#
# Backend: set USE_GPU = true and load CUDA (or AMDGPU) to run on a GPU.

const USE_GPU = false

using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    using CUDA
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
using Printf

# ---------------------------------------------------------------------------
# BP1 parameters
# ---------------------------------------------------------------------------
Base.@kwdef struct BP1Params
    rho::Float64 = 2670.0
    cs::Float64 = 3464.0
    mu::Float64 = rho * cs^2
    sigma_n::Float64 = 50.0e6
    a0::Float64 = 0.010
    amax::Float64 = 0.025
    b::Float64 = 0.015
    Dc::Float64 = 0.008
    V0::Float64 = 1.0e-6
    f0::Float64 = 0.6
    Vp::Float64 = 1.0e-9
    Vinit::Float64 = 1.0e-9
    H::Float64 = 15.0e3
    h::Float64 = 3.0e3
    Wf::Float64 = 40.0e3
    eta::Float64 = mu / (2.0 * cs)
end

function a_profile(z, p::BP1Params)
    a = fill(p.amax, length(z))
    for (k, zk) in enumerate(z)
        if zk < p.H
            a[k] = p.a0
        elseif zk < p.H + p.h
            a[k] = p.a0 + (p.amax - p.a0) * (zk - p.H) / p.h
        end
    end
    return a
end

# scalar friction functions, used inside kernels
@inline function friction(V, theta, a, sigma_n, f0, b, V0, Dc)
    arg = (f0 + b * log(V0 * theta / Dc)) / a
    return sigma_n * a * asinh(V / (2.0 * V0) * exp(arg))
end

@inline function friction_dV(V, theta, a, sigma_n, f0, b, V0, Dc)
    arg = (f0 + b * log(V0 * theta / Dc)) / a
    q = exp(arg) / (2.0 * V0)
    return sigma_n * a * q / sqrt(1.0 + (V * q)^2)
end

@inline function theta_aging_exact(theta_old, V, dt, Dc)
    x = V * dt / Dc
    ex = x < 1e-8 ? 1.0 - x : exp(-x)
    return Dc / V + (theta_old - Dc / V) * ex
end

# ---------------------------------------------------------------------------
# Elastic residual and PT update
# ---------------------------------------------------------------------------
# Interior nodes i = 2..nx, j = 1..nz (j = 1 is the free surface, treated by mirror).
# Dirichlet at i = 1 (fault), i = nx+1, j = nz+1.
@parallel_indices (i, j) function compute_residual!(R, u, mu, _dx2, _dz2, nx, nz)
    if 2 <= i <= nx && j <= nz
        uxx = (u[i + 1, j] - 2.0 * u[i, j] + u[i - 1, j]) * _dx2
        if j == 1
            uzz = 2.0 * (u[i, 2] - u[i, 1]) * _dz2
        else
            uzz = (u[i, j + 1] - 2.0 * u[i, j] + u[i, j - 1]) * _dz2
        end
        R[i, j] = mu * (uxx + uzz)
    end
    return nothing
end

@parallel_indices (i, j) function update_u!(u, dudtau, R, damp, dtau, _mu, nx, nz)
    if 2 <= i <= nx && j <= nz
        dudtau[i, j] = damp * dudtau[i, j] + dtau * R[i, j] * _mu
        u[i, j] += dtau * dudtau[i, j]
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Fault: nonlinear Robin condition solved per node by safeguarded Newton in log V
# ---------------------------------------------------------------------------
# One thread per RSF fault node. rsf_idx holds the z indices of the RSF nodes
# (1-based, stored as Float64), so the kernel runs over k = 1:length(rsf_idx).
@parallel_indices (k) function solve_fault!(V, u, u0_old, theta_old, a_f, tau0, rsf_idx,
        dt, mu, dx, sigma_n, f0, b, V0, Dc, eta)
    j = Int(rsf_idx[k])
    c0 = tau0[j] + mu * (4.0 * u[2, j] - u[3, j]) / (2.0 * dx)
    cu = -3.0 * mu / (2.0 * dx)
    lo = log(1e-25)
    hi = log(1e3)
    s = log(clamp(V[j], 1e-25, 1e3))
    aj = a_f[j]
    thj = theta_old[j]
    for _ in 1:100
        Vs = exp(s)
        gv = c0 + cu * (u0_old[j] + 0.5 * Vs * dt) -
             friction(Vs, thj, aj, sigma_n, f0, b, V0, Dc) - eta * Vs
        dg = 0.5 * cu * dt - friction_dV(Vs, thj, aj, sigma_n, f0, b, V0, Dc) - eta
        if gv > 0.0
            lo = s
        else
            hi = s
        end
        ds = clamp(-gv / (dg * Vs), -5.0, 5.0)
        s_raw = s + ds
        s_cl = clamp(s_raw, lo, hi)
        s_new = (s_cl != s_raw && s_cl == s) ? 0.5 * (lo + hi) : s_cl
        done = abs(s_new - s) < 1e-12 || abs(gv) < 1e-4 || (hi - lo) < 1e-12
        s = s_new
        done && break
    end
    V[j] = exp(s)
    u[1, j] = u0_old[j] + 0.5 * V[j] * dt
    return nothing
end

@parallel_indices (k) function update_theta!(theta, theta_old, V, rsf_idx, dt, Dc)
    j = Int(rsf_idx[k])
    theta[j] = theta_aging_exact(theta_old[j], V[j], dt, Dc)
    return nothing
end

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
function run_bp1(; nx = 80, nz = 80, Lx = 80e3, Lz = 80e3, t_end_yr = 400.0,
        tol_pa = 1.0, xi = 0.1, nu = 4.0, cfl = 0.9, out_every = 25, logfile = "bp1_log_jl.csv",
        check_every = 10, wall_limit = Inf)
    p = BP1Params()
    yr = 365.25 * 86400.0
    dx, dz = Lx / nx, Lz / nz
    z = collect(range(0.0, Lz, length = nz + 1))
    rsf_cpu = z .<= p.Wf
    rsf_idx_cpu = findall(rsf_cpu)            # RSF node indices along z
    creep_idx_cpu = findall(.!rsf_cpu)        # nodes below Wf, creeping at Vp
    nrsf = length(rsf_idx_cpu)
    a_cpu = a_profile(z, p)
    theta0 = p.Dc / p.Vinit
    tau0_cpu = [rsf_cpu[j] ? friction(p.Vinit, theta0, a_cpu[j], p.sigma_n, p.f0, p.b, p.V0, p.Dc) + p.eta * p.Vinit : 0.0
                for j in 1:nz+1]

    # device arrays
    u = @zeros(nx + 1, nz + 1)
    u_prev = @zeros(nx + 1, nz + 1)
    u_guess = @zeros(nx + 1, nz + 1)
    R = @zeros(nx + 1, nz + 1)
    dudtau = @zeros(nx + 1, nz + 1)

    V = @zeros(nz + 1)
    V[rsf_idx_cpu] .= p.Vinit
    # V = Data.Array(fill(p.Vinit, nz + 1))

theta = Data.Array(fill(theta0, nz + 1))
    theta_old = copy(theta)
    u0_old = @zeros(nz + 1)
    a_f = Data.Array(a_cpu)
    tau0 = Data.Array(tau0_cpu)
    # Data.Array is fixed to the number type given in @init_parallel_stencil, so the
    # RSF node indices are stored as Float64 on the device and converted in the kernels
    rsf_idx = Data.Array(Float64.(rsf_idx_cpu))

    dtau = cfl * min(dx, dz) / sqrt(2.0)
    damp = 1.0 - nu / max(nx, nz)
    _dx2, _dz2, _mu = 1.0 / dx^2, 1.0 / dz^2, 1.0 / p.mu
    max_it = 50 * max(nx, nz)
    dt_max = 0.1 * yr

    t = 0.0
    dt_prev = -1.0
    step = 0
    io = open(logfile, "w")
    println(io, "step,t_yr,dt_s,it,resid_pa,Vmax")
    t_wall0 = time()

    while t < t_end_yr * yr
        Vmax = maximum(Array(V))
        dt = min(dt_max, xi * p.Dc / Vmax)
        t_new = t + dt
        u_far = 0.5 * p.Vp * t_new

        # warm start by linear extrapolation
        if dt_prev > 0
            @. u_guess = u + (u - u_prev) * (dt / dt_prev)
        else
            u_guess .= u
        end
        u_prev .= u
        u0_old .= view(u, 1, :)
        theta_old .= theta

        u .= u_guess
        u[nx + 1, :] .= u_far
        u[:, nz + 1] .= u_far
        # creeping fault nodes below Wf
        u[1, creep_idx_cpu] .= u_far
        dudtau .= 0.0

        it = 0
        err = Inf
        while true
            @parallel (1:nrsf) solve_fault!(V, u, u0_old, theta_old, a_f, tau0, rsf_idx,
                dt, p.mu, dx, p.sigma_n, p.f0, p.b, p.V0, p.Dc, p.eta)
            @parallel (2:nx, 1:nz) compute_residual!(R, u, p.mu, _dx2, _dz2, nx, nz)
            @parallel (2:nx, 1:nz) update_u!(u, dudtau, R, damp, dtau, _mu, nx, nz)
            it += 1
            # check_every = 1 reproduces the Python iteration counts exactly;
            # a larger value saves reductions at the cost of a few extra iterations
            if it % check_every == 0
                err = maximum(abs.(view(R, 2:nx, 1:nz))) * dx
                (err < tol_pa || it >= max_it) && break
            end
        end

        @parallel (1:nrsf) update_theta!(theta, theta_old, V, rsf_idx, dt, p.Dc)
        t = t_new
        dt_prev = dt
        step += 1

        if step % out_every == 0
            Vm = maximum(Array(V))
            @printf(io, "%d,%.10e,%.10e,%d,%.10e,%.16e\n", step, t / yr, dt, it, err, Vm)
            @printf("step %6d t=%9.4f yr dt=%9.3e s it=%5d res=%8.2e Pa Vmax=%9.3e\n",
                step, t / yr, dt, it, err, Vm)
            flush(io)
        end
        if time() - t_wall0 > wall_limit
            println("wall limit reached")
            break
        end
    end
    close(io)
    return nothing
end

# run_bp1()
