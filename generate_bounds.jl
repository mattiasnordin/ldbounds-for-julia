# This code builds on the ldbounds package for R,
# by Charlie Casper, Thomas Cook and Oscar A. Perez
# see https://CRAN.R-project.org/package=ldbounds.
# That package, in turn, builds on  FORTRAN program ld98
# implemented by Reboussin, et al. (2000), doi:10.1016/s0197-2456(00)00057-x

using Distributions

function qp(xq, nlast, nints, yam1, ybm1, stdv)
    hlast = (ybm1 - yam1) / nints
    grid = LinRange(yam1, ybm1, nints + 1)
    fun = nlast .* cdf.(Normal(xq, stdv), grid)
    qp = 0.5 * hlast * (2 * sum(fun) - fun[1] - fun[length(fun)])
    return qp
end

function num_find(nlast, nints, target, stdv, ya, yb)
    tol = 1e-20
    l_est = BigFloat(yb - 50)
    u_est = BigFloat(yb + 50)
    pl = qp(l_est, nlast, nints, ya, yb, stdv)
    pu = qp(u_est, nlast, nints, ya, yb, stdv)
    while (abs(pu - target) > tol) & (abs(pl - target) > tol)
        println(u_est, " ", l_est)
        t = BigFloat((u_est + l_est) / 2)
        println(t)
        np = qp(t, nlast, nints, ya, yb, stdv)
        if np > target + tol
            l_est = t
            pl = np
        elseif np < target - tol
            u_est = t
            pu = np
        end
    end
    if abs(pl - target) < abs(pu - target)
        return Float64(l_est)
    else
        return Float64(u_est)
    end
end

function yu_search(nlast, nints, pd, stdv, ya, yb)
    tol = 10^(-15)
    del = 10
    max_search = 500
    uppr = yb
    q = qp(uppr, nlast, nints, ya, yb, stdv)
    while (abs(q - pd) > tol)
        del /= 10
        incr = 2 * (q > pd + tol) - 1
        j = 1
        while j <= max_search
            uppr += incr * del
            q = qp(uppr, nlast, nints, ya, yb, stdv)
            if (abs(q - pd) > tol) & (j == max_search)
                error("Convergence not achieved")
            elseif ((incr == 1) & (q <= pd + tol)) |
                   ((incr == -1) & (q >= pd - tol))
                j = max_search
            end
            j += 1
        end
    end
    ybval = uppr
    return ybval
end

function get_bounds(t_list, p_list)
    zl = -8
    h = .05
    stdv = sqrt(t_list[1])

    zu = quantile(Normal(), 1 - p_list[1])
    if zu > 8
        zu = 8.0
    end
    z_list = [zu]
    yl = zl * stdv
    yu = zu * stdv
    nh = Int.(ceil((yu - yl) / (h * stdv)))

    grid = LinRange(yl, yu, nh + 1)
    flast = pdf.(Normal(0, stdv), grid)

    for i = 2:length(t_list)
        yup = deepcopy(yu)
        ylp = deepcopy(yl)
        yu = yu_search(flast, nh, p_list[i],
                       sqrt(t_list[i] - t_list[i-1]), yl, yu)
        yl = zl * sqrt(t_list[i])
        zu = yu / sqrt(t_list[i])
        push!(z_list, zu)

        nh_last = deepcopy(nh)
        nh = Int.(ceil((yu - yl) / (h * sqrt(t_list[i] - t_list[i-1]))))
        grid = LinRange(yl, yu, nh + 1)

        flast = get_new_f(flast, grid, ylp, yup,
                          sqrt(t_list[i] - t_list[i-1]), nh_last)
    end
    return z_list
end

function get_bounds3(t_list, p_list)
    zl = -8
    h = .05
    stdv = sqrt(t_list[1])

    zu = quantile(Normal(), 1 - p_list[1])
    if zu > 8
        zu = 8.0
    end
    z_list = [zu]
    yl = zl * stdv
    yu = zu * stdv
    nh = Int.(ceil((yu - yl) / (h * stdv)))

    grid = LinRange(yl, yu, nh + 1)
    flast = pdf.(Normal(0, stdv), grid)

    for i = 2:length(t_list)
        if p_list[i] > 1e-14
            yup = deepcopy(yu)
            ylp = deepcopy(yl)
            yu = yu_search(flast, nh, p_list[i],
                           sqrt(t_list[i] - t_list[i-1]), yl, yu)
            yl = zl * sqrt(t_list[i])
            zu = yu / sqrt(t_list[i])

            nh_last = deepcopy(nh)
            nh = Int.(ceil((yu - yl) / (h * sqrt(t_list[i] - t_list[i-1]))))
            grid = LinRange(yl, yu, nh + 1)

            flast = get_new_f(flast, grid, ylp, yup,
                              sqrt(t_list[i] - t_list[i-1]), nh_last)
        else
            zu = z_list[i-1] +
            (z_list[i-1] - z_list[i-2])^2 / (z_list[i-2] - z_list[i-3])
        end
        push!(z_list, zu)
    end
    return z_list
end

function get_bounds2(t_list, p_list)
    zl = -8
    h = .05
    stdv = sqrt(t_list[1])

    zu = quantile(Normal(), 1 - p_list[1])
    if zu > 8
        zu = 8.0
    end
    z_list = [zu]
    yl = zl * stdv
    yu = zu * stdv
    nh = Int.(ceil((yu - yl) / (h * stdv)))

    grid = LinRange(yl, yu, nh + 1)
    flast = pdf.(Normal(0, stdv), grid)

    for i = 2:length(t_list)
        yup = deepcopy(yu)
        ylp = deepcopy(yl)
        yu = num_find(BigFloat.(flast), nh, p_list[i],
                      sqrt(t_list[i] - t_list[i-1]), yl, yu)
        yl = zl * sqrt(t_list[i])
        zu = yu / sqrt(t_list[i])
        push!(z_list, zu)
        println(zu)

        nh_last = deepcopy(nh)
        nh = Int.(ceil((yu - yl) / (h * sqrt(t_list[i] - t_list[i-1]))))
        grid = LinRange(yl, yu, nh + 1)

        flast = get_new_f(flast, grid, ylp, yup,
                          sqrt(t_list[i] - t_list[i-1]), nh_last)
    end
    return z_list
end

function get_new_f(fold, grid, yl_prev, yu_prev, stdv, nints)
    hlast = (yu_prev - yl_prev) / nints
    v = hlast .* collect(0:1:nints) .+ yl_prev
    out = zeros(length(grid))
    for j = 1:nints+1
        a = fold[j] .* pdf.(Normal.(grid, stdv), v[j])
        if (j == 1) | (j == nints+1)
            out .+= 0.5 .* a
        else
            out .+= a
        end
    end
    return out .* hlast
end

function gen_alpha_list(alpha, n_list, spending_function)
    t_list = accumulate(+, n_list) ./ sum(n_list)
    if spending_function == "O'Brien-Fleming"
        pe_list = 2 .* (1 .- cdf.(Normal(), quantile(Normal(),
                  (1-alpha/2)) ./ sqrt.(t_list)))
    elseif spending_function == "Pocock"
        pe_list = alpha .* log.(1 .+ (exp(1)-1) .* t_list)
    elseif spending_function == "Uniform"
        pe_list = alpha .* deepcopy(t_list)
    end
    pd_list = pe_list .- vcat([0], pe_list[1:length(pe_list)-1])
    if minimum(pd_list) == 0.0
        error("Floating-point error")
    end
    return pd_list
end

function gen_n_lists(n_start, n_end, n_inc)
    n_list = vcat((n_start), fill(n_inc, Int((n_end - n_start) / n_inc)))
    n_accum_list = accumulate(+, n_list)
    t_list = n_accum_list / sum(n_list)
    return n_list, n_accum_list, t_list
end
