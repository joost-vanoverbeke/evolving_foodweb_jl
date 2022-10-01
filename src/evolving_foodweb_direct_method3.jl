
# module Evolving_foodweb

using Random
using StatsBase
using LinearAlgebra
using Distributions
# using Profile

# using ExportAll

struct Init_values
    init_resource::Float64
    init_N::Int
        function Init_values(; init_resource = 200, init_N = 10000)
            return new(init_resource, init_N)
        end
end


struct Ecol_parameters
    # patches
    grid::NamedTuple{(:X, :Y),Tuple{Int,Int}}
    torus::NamedTuple{(:X, :Y),Tuple{Symbol,Symbol}}
    env_range::NamedTuple{(:X, :Y),Tuple{Tuple{Float64, Float64},Tuple{Float64, Float64}}}
    env_step::Float64
    dt_env::Float64
    patches::Int
    # dispersal
    m::Float64
    rho::Float64
    m_tl::Symbol
    m_power::Float64
    # resource
    in_rate::Float64
    out_rate::Float64
    # mortality
    d::Float64
    d_power::Float64
    # species
    species::Int
    trophic_levels::Int
    # bodymass classes
    bodymass_tl::Array{Float64, 1}
    tl_species::Vector{UInt}
    # feeding
    uptake_pars::Vector{Float64}
    i_power::Float64
    uptake_tl::Array{Float64, 2}
    resource_conversion::Float64
    resource_assimilation::Float64
    assimilation_eff::Float64
    conversion_tl::Vector{Float64}
    scale_uptake::Float64
    scale_assim::Float64

    function Ecol_parameters(;
        # patches
        grid = (X = 5, Y = 2), torus = (X = :NO, Y = :NO), env_range = (X = (-1., 1.), Y = (0., 0.)), env_step = 0., dt_env = 0.1,
        # dispersal
        m = 0.1, rho = 2., m_tl = :EQUAL,
        # resource
        in_rate = 200., out_rate = 0.1,
        # mortality
        d = 0.1, d_power = -0.25,
        # species
        species = 3, trophic_levels = 3,
        # feeding
        uptake_pars = [0.0001, 0.4, 1.], i_power = 0.75, resource_conversion = 1., resource_assimilation = 10., assimilation_eff = 0.7, scale_uptake = 1., scale_assim = 0.)

        patches = grid.X * grid.Y
        m_power = m_tl == :DECR ? -0.25 : (m_tl == :INCR ? 0.25 : 0.)
        in_rate = in_rate*scale_uptake
        uptake_pars[1] /= scale_uptake
        species = patches*trophic_levels
        tl_species = [(s - 1) % trophic_levels + 1 for s = 1:species]
        bodymass_tl = [10^(l - 1) for l = 1:trophic_levels]
        uptake_tl = zeros(trophic_levels, 2)
        uptake_tl[:, 1] .= @. uptake_pars[1] * bodymass_tl^i_power
        for tl in 1:trophic_levels
            if tl == 1
                uptake_tl[tl, 2] = uptake_pars[2] * bodymass_tl[tl]^(-i_power) * resource_conversion
            else
                uptake_tl[tl, 2] = uptake_pars[2] * bodymass_tl[tl]^(-i_power) * bodymass_tl[tl-1]
            end
        end
        conversion_tl = zeros(trophic_levels)
        for tl in 1:trophic_levels
            if tl == 1
                conversion_tl[tl] = resource_assimilation/bodymass_tl[tl]
            else
                conversion_tl[tl] = assimilation_eff*(1-scale_assim*bodymass_tl[tl]^(-i_power)) * bodymass_tl[tl-1]/bodymass_tl[tl]
            end
        end
        return new(grid, torus, env_range, env_step, dt_env, patches, m, rho, m_tl, m_power, in_rate, out_rate, d, d_power, species, trophic_levels, bodymass_tl, tl_species, uptake_pars, i_power, uptake_tl, resource_conversion, resource_assimilation, assimilation_eff, conversion_tl, scale_uptake, scale_assim)
    end
end


struct Evol_parameters
    omega_e::Float64
    trait_loci::Int
    tot_genes::Int
    mu::Float64
    sigma_z::Float64
    div_f::Float64
    all_mother::Array{Int, 1}
    all_father::Array{Int, 1}
    all_genes::Array{Int, 1}

    function Evol_parameters(; omega_e = 4.0, trait_loci = 20, mu = 1e-4, sigma_z = 0.1)
        div_f = 2 * omega_e^2
        tot_genes = 2 * trait_loci
        all_mother = 1:trait_loci
        all_father = all_mother .+ trait_loci
        all_genes = [all_mother; all_father]
        return new(omega_e, trait_loci, tot_genes, mu, sigma_z, div_f, all_mother, all_father, all_genes)
    end
end


mutable struct Direct_method
    c_b_resource::Vector{Float64}
    c_d_resource::Vector{Float64}
    c_b_ind::Array{Float64, 2}
    c_d_ind::Array{Float64, 2}
    c_vec::Vector{Float64}
    cs_vec::Vector{Float64}
    c_total::Float64
    c_b_res_pos::Vector{Int}
    c_d_res_pos::Vector{Int}
    c_b_ind_pos::Array{Int, 2}
    c_d_ind_pos::Array{Int, 2}
    c_cart_inds::CartesianIndices{2,Tuple{Base.OneTo{Int64},Base.OneTo{Int64}}}

    function Direct_method(ecol::Ecol_parameters; c_b_res = 200., c_d_res = 500., c_b = 1., c_d = 1.)
        c_b_resource = fill(c_b_res, ecol.patches)
        c_d_resource = fill(c_d_res, ecol.patches)
        c_b_ind = [c_b/1.3^ecol.tl_species[s] for s in 1:ecol.species, p in 1:ecol.patches]
        c_d_ind = [c_d/1.3^ecol.tl_species[s] for s in 1:ecol.species, p in 1:ecol.patches]
        c_vec = [c_b_resource; c_d_resource; c_b_ind[:]; c_d_ind[:]]
        cs_vec = cumsum(c_vec)
        c_total = cs_vec[end]
        c_b_res_pos = collect(1:ecol.patches)
        c_d_res_pos = ecol.patches .+ collect(1:ecol.patches)
        c_b_ind_pos = [(p-1)*ecol.species .+ s .+ 2*ecol.patches for s in 1:ecol.species, p in 1:ecol.patches]
        c_d_ind_pos = c_b_ind_pos .+ ecol.patches*ecol.species
        c_cart_inds = CartesianIndices(c_b_ind)
        return new(c_b_resource, c_d_resource, c_b_ind, c_d_ind, c_vec, cs_vec, c_total, c_b_res_pos, c_d_res_pos, c_b_ind_pos, c_d_ind_pos, c_cart_inds)
    end
end

function update_dm_patch(dm::Direct_method, ecol::Ecol_parameters, patch)
    p = patch.patch_ID
    dm.c_b_resource[p] = patch.resource.gain
    dm.c_d_resource[p] = patch.resource.loss
    # dm.c_total -= dm.c_vec[dm.c_b_res_pos[p]]
    # dm.c_total -= dm.c_vec[dm.c_d_res_pos[p]]
    dm.c_vec[dm.c_b_res_pos[p]] = dm.c_b_resource[p]
    dm.c_vec[dm.c_d_res_pos[p]] = dm.c_d_resource[p]
    # dm.c_total += dm.c_vec[dm.c_b_res_pos[p]]
    # dm.c_total += dm.c_vec[dm.c_d_res_pos[p]]
    for s in 1:ecol.species
        if patch.N_s[s] > 0
            tl = ecol.tl_species[s]
            dm.c_b_ind[s,p] = patch.gain_tl[tl]
            dm.c_d_ind[s,p] = patch.loss_tl[tl] + patch.mort_s[s]
        else
            dm.c_b_ind[s,p] = 0
            dm.c_d_ind[s,p] = 0
        end
        # dm.c_total -= dm.c_vec[dm.c_b_ind_pos[s,p]]
        # dm.c_total -= dm.c_vec[dm.c_d_ind_pos[s,p]]
        dm.c_vec[dm.c_b_ind_pos[s,p]] = dm.c_b_ind[s,p].*patch.N_s[s]
        dm.c_vec[dm.c_d_ind_pos[s,p]] = dm.c_d_ind[s,p].*patch.N_s[s]
        # dm.c_total += dm.c_vec[dm.c_b_ind_pos[s,p]]
        # dm.c_total += dm.c_vec[dm.c_d_ind_pos[s,p]]
    end
    # cumsum!(dm.cs_vec, dm.c_vec)
    # dm.c_total = dm.cs_vec[end]
    dm.c_total = sum(dm.c_vec)
end

function next_time(dm::Direct_method)
    t_next = rand(Exponential(1. / dm.c_total))
    return t_next
end

function sample_c(dm::Direct_method, ecol::Ecol_parameters)
    # c_sel = binary_search((@view dm.cs_vec[:]), rand()*dm.c_total)
    # c_sel = wsample((@view dm.c_vec[:]))
    # c_sel = my_sample(dm.c_vec, dm.c_total)
    c_sel = my_sample3(dm.c_vec, dm.c_total)
    event = :nothing
    s = 0
    p = 0
    if 1 <= c_sel <= ecol.patches
        event = :resource_gain
        s = 0
        p = c_sel
    elseif (ecol.patches + 1) <= c_sel <= (2*ecol.patches)
        event = :resource_loss
        s = 0
        p = c_sel - ecol.patches
    elseif (2*ecol.patches + 1) <= c_sel <= (2*ecol.patches + ecol.patches*ecol.species)
        event = :reproduction
        c_sel -= 2*ecol.patches
        ci = dm.c_cart_inds[c_sel]
        s = ci[1]
        p = ci[2]
    else
        event = :death
        c_sel -= 2*ecol.patches + ecol.patches*ecol.species
        ci = dm.c_cart_inds[c_sel]
        s = ci[1]
        p = ci[2]
    end
    return event, s, p
end

# struct Run
#     dt::Float64
#     runs::Int
#     time_steps::Int
#     pre_change::Int
#     post_change::Int
#     print_steps::Int
#     save_steps::Int
#     file::String
#
#     function Run(;
#         dt = 0.1,
#         runs = 1,
#         time_steps = 100,
#         pre_change = 10,
#         post_change = 10,
#         print_steps = 10,
#         save_steps = 100,
#         file = "output_evolving_foodweb.csv",
#     )
#         return new(dt, runs, time_steps, pre_change, post_change, print_steps, save_steps, file)
#     end
# end


mutable struct Resource
    in_rate::Float64
    out_rate::Float64
    level::Float64
    gain::Float64
    loss::Float64

    function Resource(ecol::Ecol_parameters, level = 200.)
        return new(ecol.in_rate, ecol.out_rate, level, 0., 0.)
    end
end

function update_gain(resource::Resource)
    resource.gain = resource.in_rate
end

function update_loss(resource::Resource, uptake)
    resource.loss = resource.out_rate*resource.level + uptake
end

function increase_level(resource::Resource, p, dm::Direct_method)
    event = rand() <= (resource.gain/dm.c_b_resource[p])
    if event
        resource.level += 1
    end
    return event
end

function decrease_level(resource::Resource, p, dm::Direct_method)
    event = rand() <= (resource.loss/dm.c_d_resource[p])
    if event
        resource.level -= 1
    end
    return event
end


mutable struct Individual
    species_ID::Int
    trophic_level::UInt
    bodymass::Float64
    genotype::Vector{Int8}
    phenotype::Float64
    fitness::Float64
    mortality_rate::Float64

    function Individual(ecol::Ecol_parameters, evol::Evol_parameters, e::Float64, s::Int, gtp::Vector{Int8})
        species_ID = s
        trophic_level = ecol.tl_species[s]
        bodymass = 10.0^(trophic_level - 1)
        genotype = gtp
        phenotype = mean(genotype) + randn() * evol.sigma_z
        fitness = exp(-(phenotype - e)^2 / evol.div_f)
        mortality_rate = 1 - fitness * (1 -  ecol.d * bodymass^ecol.d_power)
        return new(species_ID, trophic_level, bodymass, genotype, phenotype, fitness, mortality_rate)
    end

    function Individual(ecol::Ecol_parameters, evol::Evol_parameters, e = 0.; s = 1)
        genotype = Int8.(round.(rand(evol.tot_genes) .* 0.5 .* rand([-1.0, 1.0], evol.tot_genes) .+ e))
        return Individual(ecol, evol, e, s, genotype)
    end

    function Individual(ecol::Ecol_parameters, evol::Evol_parameters, mother::Individual, e = 0.)
        species_ID = mother.species_ID
        genotype = [i for i in mother.genotype]
        muts = rand(evol.tot_genes) .<= evol.mu
        muts_sum = sum(muts)
        if muts_sum > 0
            genotype[muts] .+= rand([-1,1] , muts_sum)
        end
        return Individual(ecol, evol, e, species_ID, genotype)
    end

    function Individual(ecol::Ecol_parameters, evol::Evol_parameters, mother::Individual, father::Individual, e = 0.)
        species_ID = mother.species_ID
        genotype = zeros(Int8, evol.tot_genes)
        genes = rand(evol.trait_loci) .< 0.5
        genotype[evol.all_mother[genes]] .= mother.genotype[evol.all_mother[genes]]
        genotype[evol.all_mother[.!genes]] .= mother.genotype[evol.all_father[.!genes]]
        genes = rand(evol.trait_loci) .< 0.5
        genotype[evol.all_father[genes]] .= father.genotype[evol.all_mother[genes]]
        genotype[evol.all_father[.!genes]] .= father.genotype[evol.all_father[.!genes]]
        muts = rand(evol.tot_genes) .<= evol.mu
        muts_sum = sum(muts)
        if muts_sum > 0
            genotype[muts] .+= rand([-1,1] , muts_sum)
        end
        return Individual(ecol, evol, e, species_ID, genotype)
    end
end

function update_fitness(ind::Individual, env::Float64, ecol::Ecol_parameters, evol::Evol_parameters)
    ind.fitness =
        exp(-(ind.phenotype-env)^2 / evol.div_f)
    ind.mortality_rate = 1 - ind.fitness * (1 -  ecol.d * ind.bodymass^ecol.d_power)
    return ind.mortality_rate
end


mutable struct Patch
    patch_ID::Int
    environment::Float64
    resource::Resource
    individuals::Vector{Vector{Individual}}
    morts::Vector{Vector{Float64}}
    N_s::Vector{Int}
    N_tl::Vector{Int}
    total_N::Int
    total_uptake::Vector{Float64}
    gain_tl::Vector{Float64}
    loss_tl::Vector{Float64}
    mort_s::Vector{Float64}

    function Patch(ecol::Ecol_parameters, evol::Evol_parameters, init::Init_values, id = 1, e = 0.)
        patch_ID = id
        environment = e
        resource = Resource(ecol, init.init_resource)
        individuals = [Individual[] for i in 1:ecol.species]
        morts = [Float64[] for i in 1:ecol.species]
        N_s = zeros(Int64, ecol.species)
        N_tl = zeros(Int64, ecol.trophic_levels)
        total_uptake = zeros(ecol.trophic_levels)
        prob_species = [((s - 1) ÷ (ecol.species / ecol.patches) == (patch_ID - 1) ? 1 : 0) /
            10^(1.5*(ecol.tl_species[s]-1)) for s in 1:ecol.species]
        # prob_species = [1,0,0]
        s_id = wsample(1:ecol.species, prob_species, init.init_N)
        for i in 1:init.init_N
            push!(individuals[s_id[i]], Individual(ecol, evol, environment; s = s_id[i]))
            push!(morts[s_id[i]], individuals[s_id[i]][end].mortality_rate)
            N_s[s_id[i]] += 1
            N_tl[ecol.tl_species[s_id[i]]] += 1
        end
        total_N = sum(N_tl)
        gain_tl = zeros(ecol.species)
        loss_tl = zeros(ecol.species)
        mort_s = [N_s[s] > 0 ? maximum(morts[s]) : 0 for s in 1:ecol.species]
        return new(patch_ID, environment, resource, individuals, morts, N_s, N_tl, total_N, total_uptake, gain_tl, loss_tl, mort_s)
    end
end

function update_uptake(patch::Patch, ecol::Ecol_parameters, dm::Direct_method)
    for tl in 1:ecol.trophic_levels
        if tl == 1
            upt = ecol.uptake_tl[tl, 1]*patch.resource.level^ecol.uptake_pars[3]
        else
            upt = ecol.uptake_tl[tl, 1]*patch.N_tl[tl-1]^ecol.uptake_pars[3]
        end
        patch.total_uptake[tl] = upt/(1. + upt*ecol.uptake_tl[tl, 2])*patch.N_tl[tl]
    end
    for tl in 1:ecol.trophic_levels
        patch.gain_tl[tl] = patch.total_uptake[tl]/patch.N_tl[tl] * ecol.conversion_tl[tl]
        if tl < ecol.trophic_levels
            patch.loss_tl[tl] = patch.total_uptake[tl+1]/patch.N_tl[tl]
        end
    end
    update_gain(patch.resource)
    update_loss(patch.resource, patch.total_uptake[1])
    update_dm_patch(dm, ecol, patch)
end

function change_environment(patch::Patch, ecol::Ecol_parameters, evol::Evol_parameters)
    patch.environment += ecol.env_step * ecol.dt_env
    for s in 1:ecol.species
        if patch.N_s[s] > 0
            # patch.morts[s] = map(i -> update_fitness(i, patch.environment, ecol, evol), patch.individuals[s])
            for i in 1:patch.N_s[s]
                patch.morts[s][i] = update_fitness(patch.individuals[s][i], patch.environment, ecol, evol)
            end
            patch.mort_s[s] = maximum(patch.morts[s])
        end
    end
end


mutable struct World
    # patches
    nbr_patches::Int
    patch_XY::Array{Int, 2}
    patches::Vector{Patch}
    # dispersal
    neighbours::Array{Float64, 2}
    m_neighbours::Array{Float64, 3}
    cs_m_neighbours::Array{Float64, 3}

    function World(ecol::Ecol_parameters, evol::Evol_parameters, init::Init_values)
        # patches
        grid = ecol.grid
        nbr_patches = ecol.patches
        patch_XY = zeros(Int, 2, nbr_patches)
        for x in (1:grid.X) .- 1, y in (1:grid.Y) .- 1
            patch_XY[1, y*grid.X+x+1] = x + 1
            patch_XY[2, y*grid.X+x+1] = y + 1
        end
        # dispersal
        neighbours = zeros(nbr_patches, nbr_patches)
        for x in (1:grid.X), y in (1:grid.Y), x2 in (1:grid.X), y2 in (1:grid.Y)
            d_X = abs(x - x2)
            d_X = (ecol.torus.X == :YES ? min(d_X, grid.X - d_X) : d_X)^2
            d_Y = abs(y - y2)
            d_Y = (ecol.torus.Y == :YES ? min(d_Y, grid.Y - d_Y) : d_Y)^2
            neighbours[(y-1)*grid.X+x, (y2-1)*grid.X+x2] = sqrt(d_X + d_Y)
        end
        temp_m_neighbours = @. ecol.rho * exp(-ecol.rho * neighbours)
        temp_m_neighbours[diagind(temp_m_neighbours)] .= 0.
        temp_m_neighbours ./= sum(temp_m_neighbours; dims = 1)
        m_neighbours = repeat(temp_m_neighbours, outer = [1, 1, ecol.trophic_levels])
        disp_tl = ecol.m .* ecol.bodymass_tl .^ ecol.m_power
        m_neighbours .= [
            i == j ? 1 - disp_tl[c] : disp_tl[c] * m_neighbours[i, j, c]
            for i = 1:nbr_patches, j = 1:nbr_patches, c = 1:ecol.trophic_levels
        ]
        cs_m_neighbours = cumsum(m_neighbours; dims = 1)
        # environment
        step_X =
            grid.X == 1 ? 0.0 : (ecol.env_range.X[2] - ecol.env_range.X[1]) / (grid.X - 1)
        step_Y =
            grid.Y == 1 ? 0.0 : (ecol.env_range.Y[2] - ecol.env_range.Y[1]) / (grid.Y - 1)
        environment = [
            ecol.env_range.X[1] +
            step_X * (patch_XY[1, p] - 1) +
            ecol.env_range.Y[1] +
            step_Y * (patch_XY[2, p] - 1) for p = 1:nbr_patches
        ]

        patches = [Patch(ecol, evol, init, p, environment[p]) for p in 1:nbr_patches]
        return new(nbr_patches, patch_XY, patches, neighbours, m_neighbours, cs_m_neighbours)
    end
end


function next_event(world::World, ecol::Ecol_parameters, evol::Evol_parameters, dm::Direct_method)
    dt = 0.
    smpl_event = :nothing
    event = false
    p = 0
    new_p = 0
    while !event
        dt += next_time(dm)
        smpl_event, s, p = sample_c(dm, ecol)
        patch = world.patches[p]
        if smpl_event == :resource_gain
            event = increase_level(patch.resource, p, dm)
        elseif smpl_event == :resource_loss
            event = decrease_level(patch.resource, p, dm)
        elseif smpl_event == :reproduction
            tl = ecol.tl_species[s]
            event = rand() <= (patch.gain_tl[tl]/dm.c_b_ind[s, p])
            if event
                i = rand(1:patch.N_s[s])
                ind = patch.individuals[s][i]
                new_p = binary_search((@view world.cs_m_neighbours[:, p, tl]), rand())
                new_patch = world.patches[new_p]
                new_patch.N_s[s] += 1
                new_patch.total_N += 1
                new_patch.N_tl[tl] += 1
                push!(new_patch.individuals[s], Individual(ecol, evol, ind, new_patch.environment))
                d = new_patch.individuals[s][end].mortality_rate
                push!(new_patch.morts[s], d)
                if d > new_patch.mort_s[s]
                    new_patch.mort_s[s] = d
                end
            end
        else
            tl = ecol.tl_species[s]
            i = rand(1:patch.N_s[s])
            ind = patch.individuals[s][i]
            event = rand() <= ((ind.mortality_rate + patch.loss_tl[tl])/dm.c_d_ind[s, p])
            if event
                patch.N_tl[tl] -= 1
                d = patch.individuals[s][i].mortality_rate
                ind2 = pop!(patch.individuals[s])
                d2 = pop!(patch.morts[s])
                if i < patch.N_s[s]
                    patch.individuals[s][i] = ind2
                    patch.morts[s][i] = d2
                end
                patch.N_s[s] -= 1
                patch.total_N -= 1
                if patch.N_s[s] == 0
                    patch.mort_s[s] = 0
                else
                    if d == patch.mort_s[s]
                        patch.mort_s[s] = maximum(patch.morts[s])
                    end
                end
            end
        end
    end
    update_uptake(world.patches[p], ecol, dm)
    if new_p > 0 && new_p != p
        update_uptake(world.patches[new_p], ecol, dm)
    end
    return dt, smpl_event
end


function evolving_foodweb_dm(; time_steps = 1)
    init = Init_values();
    ecol = Ecol_parameters(; grid = (X = 5, Y = 2), env_step = 0.001, in_rate = 200., scale_uptake = 1.);
    evol = Evol_parameters();
    dm = Direct_method(ecol; c_b_res = 200., c_d_res = 500., c_b = 1., c_d = 1.);
    world = World(ecol, evol, init);
    for p in 1:ecol.patches
        update_uptake(world.patches[p], ecol, dm)
    end
    time = 0.
    time_env = 0.
    events = zeros(Int, 4)
    while time < time_steps
        dt, event = next_event(world, ecol, evol, dm)
        time += dt
        if ecol.env_step > 0. && time > (time_env + ecol.dt_env)
            for p in 1:ecol.patches
                change_environment(world.patches[p], ecol, evol)
            end
            time_env += ecol.dt_env
        end
        if event == :resource_gain
            events[1] += 1
        elseif event == :resource_loss
            events[2] += 1
        elseif event == :reproduction
            events[3] += 1
        else
            events[4] += 1
        end
    end
    return time, events, world, ecol, dm
end


function binary_search(lst, val)
    low = 1
    high = length(lst)
    while low ≤ high
        mid = (low + high) ÷ 2
        if lst[mid] >= val
            high = mid - 1
        else
            low = mid + 1
        end
    end
    return low
end

function my_sample(wv, tot)
    # t = rand() * sum(wv)
    t = rand() * tot
    n = length(wv)
    i = 1
    @inbounds cw = wv[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += wv[i]
    end
    return i
end

function my_sample2(wv, tot)
    cs = cumsum(wv)
    # t = rand() * cs[end]
    t = rand() * tot
    i = binary_search(cs, t)
    return i
end

function my_sample3(wv, tot)
    # t = rand() * sum((@view wv[:]))
    t = rand() * tot
    low = 1
    high = length(wv)
    midold = 0
    smid = 0
    while low ≤ high
        mid = (low + high) ÷ 2
        if mid > midold
            smid = smid + sum((@view wv[(midold+1):mid]))
        else
            smid = smid - sum((@view wv[(mid+1):midold]))
        end
        midold = mid
        if smid >= t
            high = mid - 1
        else
            low = mid + 1
        end
    end
    return low
end



# @exportAll()
# end
# using .Evolving_foodweb

# @time time, events, maxs, world, ecol = evolving_foodweb_dm(; time_steps = 10);
@time time, events, world, ecol, dm = evolving_foodweb_dm(; time_steps = 100);
time
events
fails
@. events/(fails+events)
print(maxs[4])

# @profiler time, events, maxs, world, ecol = evolving_foodweb_dm(; time_steps = 50);
@profiler time, events, world, ecol, dm = evolving_foodweb_dm(; time_steps = 100);
# @profview time, events, fails, maxs, patch, ecol = evolving_foodweb_dm(; time_steps = 100);
time
events
fails
@. events/(fails+events)
maxs

1. ./ 1.3.^[1,2,3]
