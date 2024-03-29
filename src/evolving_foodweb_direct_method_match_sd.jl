
# module Evolving_foodweb

using Random
using StatsBase
using LinearAlgebra
using Distributions
using Statistics

# using ExportAll

mutable struct Init_values
    ## ecological values
    grid::NamedTuple{(:X, :Y),Tuple{Int,Int}}
    torus::NamedTuple{(:X, :Y),Tuple{Symbol,Symbol}}
    env_range::NamedTuple{(:X, :Y),Tuple{Tuple{Float64, Float64},Tuple{Float64, Float64}}}
    CC_vec::Vector{Float64}
    time_CC_vec::Vector{Int}
    CC::Float64
    time_CC::Int
    env_step_CC::Float64
    env_step_local::Float64
    dt_env::Float64
    # dispersal
    m_vec::Vector{Float64}
    m::Float64
    rho::Float64
    m_tl::Symbol
    # resource
    resource::Float64
    in_rate::Float64
    out_rate::Float64
    # species
    N::Int
    spec_dist::Symbol
    spec_dens::Int
    rep_type::Symbol
    # species matching
    match_sd_vec::Vector{Float64}
    match_sd::Float64
    # trophic levels
    trophic_levels::Int
    bm_offset::Float64
    bm_power::Float64
    # mortality
    d::Float64
    d_power::Float64
    # feeding
    uptake_pars::Vector{Float64}
    i_power::Float64
    resource_conversion::Float64
    resource_assimilation::Float64
    assimilation_eff::Float64
    scale_uptake::Float64
    scale_assim::Float64

    ## evolutionary values
    omega_e::Float64
    trait_loci::Int
    mu::Float64
    sigma_z::Float64

    ## run values
    runs::Int
    pre_post_change::Int
    print_steps::Int
    log_steps::Int
    output_file::String

    function Init_values(;
        ## ecological values
        # patches
        grid = (X = 5, Y = 2), torus = (X = :NO, Y = :NO), env_range = (X = (-1., 1.), Y = (0., 0.)), 
        CC_vec = [0.], time_CC_vec = 0, CC = CC_vec[1], time_CC = time_CC_vec[1], env_step_CC = CC/time_CC, env_step_local = 0.01, dt_env = 0.1,
        # dispersal
        m_vec = [0.1], m = m_vec[1], rho = 2., m_tl = :EQUAL,
        # resource
        resource = 200., in_rate = 200., out_rate = 0.1,
        # species
        N = 10000, spec_dist = :X, spec_dens = 1, rep_type = :ASEXUAL,
        # species matching
        match_sd_vec = [1.], match_sd = match_sd_vec[1], 
        # trophic levels
        trophic_levels = 3, bm_offset = 1., bm_power = 1.,
        # mortality
        d = 0.1, d_power = -0.25,
        # feeding
        uptake_pars = [0.0001, 0.4, 1.], i_power = 0.75, resource_conversion = 1., resource_assimilation = 10., assimilation_eff = 0.7, scale_uptake = 1., scale_assim = 0.,

        ## evolutionary values
        omega_e = 4.0, trait_loci = 20, mu = 1e-4, sigma_z = 0.1,

        ## run values
        runs = 1,
        pre_post_change = 10,
        print_steps = 10,
        log_steps = 10,
        output_file = "output_evolving_foodweb.csv"
        )
        return new(grid, torus, env_range, CC_vec, time_CC_vec, CC, time_CC, env_step_CC, env_step_local, dt_env, m_vec, m, rho, m_tl, resource, in_rate, out_rate, N, spec_dist, spec_dens, rep_type, match_sd_vec, match_sd, trophic_levels, bm_offset, bm_power, d, d_power, uptake_pars, i_power, resource_conversion, resource_assimilation, assimilation_eff, scale_uptake, scale_assim, omega_e, trait_loci, mu, sigma_z, runs, pre_post_change, print_steps, log_steps, output_file)
    end
end

# struct Ecol_parameters
mutable struct Ecol_parameters
    # patches
    grid::NamedTuple{(:X, :Y),Tuple{Int,Int}}
    torus::NamedTuple{(:X, :Y),Tuple{Symbol,Symbol}}
    env_range::NamedTuple{(:X, :Y),Tuple{Tuple{Float64, Float64},Tuple{Float64, Float64}}}
    CC::Float64
    env_step_CC::Float64
    time_CC::Int
    env_step_local::Float64
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
    # species
    species::Int
    rep_type::Symbol
    # trophic levels
    trophic_levels::Int
    bodymass_tl::Array{Float64, 1}
    tl_species::Vector{Int}
    species_per_tl::Int

    match_resource::Vector{Int}
    species_match::Vector{Float64}
    match_sd::Float64
    match_matr::Array{Float64, 2}
    match_rel::Array{Float64, 2}

    bm_offset::Float64
    bm_power::Float64
    # mortality
    d::Float64
    d_power::Float64
    d_tl::Vector{Float64}
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

    function Ecol_parameters(init::Init_values)
        # patches
        grid = init.grid
        torus = init.torus
        env_range = init.env_range
        CC = init.CC
        env_step_CC = init.env_step_CC
        time_CC = init.time_CC
        env_step_local = init.env_step_local
        dt_env = init.dt_env
        # dispersal
        m = init.m
        rho = init.rho
        m_tl = init.m_tl
        # resource
        in_rate = init.in_rate
        out_rate = init.out_rate
        # species
        rep_type = init.rep_type
        # species matching
        match_sd = init.match_sd
        # trophic levels
        trophic_levels = init.trophic_levels
        bm_offset = init.bm_offset
        bm_power = init.bm_power
        # mortality
        d = init.d
        d_power = init.d_power
        # feeding
        uptake_pars = copy(init.uptake_pars)
        i_power = init.i_power
        resource_conversion = init.resource_conversion
        resource_assimilation = init.resource_assimilation
        assimilation_eff = init.assimilation_eff
        scale_uptake = init.scale_uptake
        scale_assim = init.scale_assim

        patches = grid.X * grid.Y
        if patches == 1
            m = 0.
        end
        m_power = m_tl == :DECR ? -0.25 : (m_tl == :INCR ? 0.25 : 0.)
        in_rate = in_rate*scale_uptake
        uptake_pars[1] /= scale_uptake
        if init.spec_dist == :X
            species = grid.X*trophic_levels*init.spec_dens # different species for different X
        else
            species = patches*trophic_levels*init.spec_dens # different species in each patch
            # species = 1
        end
        tl_species = [(s - 1) % trophic_levels + 1 for s = 1:species]
        species_per_tl = species/trophic_levels

        match_resource = [tl_species[s] == 1 ? 1 : 0 for s in 1:species]
        # species_match = [(s - 1) ÷ trophic_levels + 1 for s = 1:species]
        # match_shuffle = shuffle(1:(species/trophic_levels))
        match_shuffle = round.(rand(convert(Int, species/trophic_levels)), digits=3)
        species_match = [match_shuffle[(s - 1) ÷ trophic_levels + 1] for s = 1:species]
        # species_match = rand(1:10, species/trophic_levels)
        match_matr = [tl_species[s1] == tl_species[s2] + 1 ? exp(-(min(abs(species_match[s1] - species_match[s2]), 1 - abs(species_match[s1] - species_match[s2])))^2 / (2*match_sd^2)) : 0 for s1 in 1:species, s2 in 1:species]
        match_rel = match_matr./[x == 0 ? 1 : x for x in sum(match_matr; dims=2)]
        
        bodymass_tl = [bm_offset*10^((l - 1)*bm_power) for l = 1:trophic_levels]
        d_tl = d .* bodymass_tl.^d_power
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
                # # conversion_tl[tl] = resource_assimilation/bodymass_tl[tl]
                conversion_tl[tl] = resource_assimilation/(bodymass_tl[tl]^(1/scale_assim))
            else
                # # conversion_tl[tl] = assimilation_eff*(1-scale_assim*bodymass_tl[tl]^(-i_power)) * bodymass_tl[tl-1]/bodymass_tl[tl]
                conversion_tl[tl] = assimilation_eff * bodymass_tl[tl-1]/(bodymass_tl[tl]^(1/scale_assim))
            end
        end
        return new(grid, torus, env_range, CC, env_step_CC, time_CC, env_step_local, dt_env, patches, m, rho, m_tl, m_power, in_rate, out_rate, species, rep_type, trophic_levels, bodymass_tl, tl_species, species_per_tl, match_resource, species_match, match_sd, match_matr, match_rel, bm_offset, bm_power, d, d_power, d_tl, uptake_pars, i_power, uptake_tl, resource_conversion, resource_assimilation, assimilation_eff, conversion_tl, scale_uptake, scale_assim)
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

    function Evol_parameters(init::Init_values)
        omega_e = init.omega_e
        trait_loci = init.trait_loci
        mu = init.mu
        sigma_z = init.sigma_z
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
        c_total = sum(c_vec)
        c_b_res_pos = collect(1:ecol.patches)
        c_d_res_pos = ecol.patches .+ collect(1:ecol.patches)
        c_b_ind_pos = [(p-1)*ecol.species .+ s .+ 2*ecol.patches for s in 1:ecol.species, p in 1:ecol.patches]
        c_d_ind_pos = c_b_ind_pos .+ ecol.patches*ecol.species
        c_cart_inds = CartesianIndices(c_b_ind)
        return new(c_b_resource, c_d_resource, c_b_ind, c_d_ind, c_vec, c_total, c_b_res_pos, c_d_res_pos, c_b_ind_pos, c_d_ind_pos, c_cart_inds)
    end
end

function update_dm_patch(dm::Direct_method, ecol::Ecol_parameters, patch)
    p = patch.patch_ID
    gain = 0.
    loss = 0.
    # dm.c_total -= dm.c_vec[dm.c_b_res_pos[p]]
    # dm.c_total -= dm.c_vec[dm.c_d_res_pos[p]]
    dm.c_vec[dm.c_b_res_pos[p]] = dm.c_b_resource[p] = patch.resource_gain
    dm.c_vec[dm.c_d_res_pos[p]] = dm.c_d_resource[p] = patch.resource_loss
    # dm.c_total += dm.c_vec[dm.c_b_res_pos[p]]
    # dm.c_total += dm.c_vec[dm.c_d_res_pos[p]]
    for s in 1:ecol.species
        # dm.c_total -= dm.c_vec[dm.c_b_ind_pos[s,p]]
        # dm.c_total -= dm.c_vec[dm.c_d_ind_pos[s,p]]
        N_s = patch.N_s[s]
        if N_s > 0
            # tl = ecol.tl_species[s]
            gain = patch.gain_s[s]
            loss = patch.loss_s[s] + patch.mort_max_s[s]
        else
            gain = 0.
            loss = 0.
        end
        dm.c_b_ind[s,p] = gain
        dm.c_d_ind[s,p] = loss
        dm.c_vec[dm.c_b_ind_pos[s,p]] = gain*N_s
        dm.c_vec[dm.c_d_ind_pos[s,p]] = loss*N_s
        # dm.c_total += dm.c_vec[dm.c_b_ind_pos[s,p]]
        # dm.c_total += dm.c_vec[dm.c_d_ind_pos[s,p]]
    end
    dm.c_total = sum(dm.c_vec)
end

function update_dm_loss(dm::Direct_method, ecol::Ecol_parameters, patch)
    p = patch.patch_ID
    loss = 0.
    for s in 1:ecol.species
        N_s = patch.N_s[s]
        if N_s > 0
            # tl = ecol.tl_species[s]
            loss = patch.loss_s[s] + patch.mort_max_s[s]
        else
            loss = 0.
        end
        dm.c_d_ind[s,p] = loss
        dm.c_vec[dm.c_d_ind_pos[s,p]] = loss*N_s
    end
    dm.c_total = sum(dm.c_vec)
end

function next_time(dm::Direct_method)

    # println("next_time: c_total = $(dm.c_total)")

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

struct Run
    runs::Int
    time_steps::Int
    pre_change::Int
    post_change::Int
    print_steps::Int
    log_steps::Int
    output_file::String

    function Run(init::Init_values)
        pre_change = init.pre_post_change
        post_change = init.pre_post_change
        time_steps = pre_change + init.time_CC + post_change
        return new(init.runs, time_steps, pre_change, post_change, init.print_steps, init.log_steps, init.output_file)
    end
end


mutable struct Patch
    patch_ID::Int
    X::Int
    Y::Int
    environment::Float64
    resource::Float64
    species_ID::Vector{Vector{Int}}
    genotype::Vector{Vector{Vector{Int8}}}
    phenotype::Vector{Vector{Float64}}
    fitness::Vector{Vector{Float64}}
    mortality::Vector{Vector{Float64}}

    N_s::Vector{Int}
    N_tl::Vector{Int}
    total_N::Int
    N_prey::Vector{Float64}
    tl_ratio::Vector{Float64}
    uptake::Vector{Float64}
    total_uptake::Vector{Float64}
    resource_gain::Float64
    resource_loss::Float64
    # gain_tl::Vector{Float64}
    # loss_tl::Vector{Float64}
    gain_s::Vector{Float64}
    loss_s::Vector{Float64}
    mort_max_s::Vector{Float64}

    function Patch(ecol::Ecol_parameters, id = 1, X = 1, Y = 1, r = 0., e = 0.)
        patch_ID = id
        environment = e
        resource = r
        species_ID = [Int[] for i in 1:ecol.species]
        genotype = [Array{Int8,1}[] for i in 1:ecol.species]
        phenotype = [Float64[] for i in 1:ecol.species]
        fitness = [Float64[] for i in 1:ecol.species]
        mortality = [Float64[] for i in 1:ecol.species]
        N_s = zeros(Int64, ecol.species)
        N_tl = zeros(Int64, ecol.trophic_levels)
        total_N = sum(N_tl)
        N_prey = zeros(Float64, ecol.species)
        tl_ratio = zeros(Float64, ecol.trophic_levels)
        uptake = zeros(Float64, ecol.species)
        total_uptake = zeros(ecol.trophic_levels)
        resource_gain = 0.
        resource_loss = 0.
        # gain_tl = zeros(ecol.trophic_levels)
        # loss_tl = zeros(ecol.trophic_levels)
        gain_s = zeros(ecol.species)
        loss_s = zeros(ecol.species)
        mort_max_s = [N_s[s] > 0 ? maximum(mortality[s]) : 0 for s in 1:ecol.species]
        return new(patch_ID, X, Y, environment, resource, species_ID, genotype, phenotype, fitness, mortality, N_s, N_tl, total_N, N_prey, tl_ratio, uptake, total_uptake, resource_gain, resource_loss, gain_s, loss_s, mort_max_s)
    end
end

function increase_resource(patch::Patch, dm::Direct_method)
    event = rand() <= (patch.resource_gain/dm.c_b_resource[patch.patch_ID])
    if event
        patch.resource += 1
    end
    return event
end

function decrease_resource(patch::Patch, dm::Direct_method)
    event = rand() <= (patch.resource_loss/dm.c_d_resource[patch.patch_ID])
    if event
        patch.resource -= 1
    end
    return event
end

function push_individual!(patch::Patch, ecol::Ecol_parameters, evol::Evol_parameters, s::Int, genotype::Vector{Int8})
    tl = ecol.tl_species[s]
    phenotype = mean(genotype) + randn() * evol.sigma_z
    fitness = exp(-(phenotype - patch.environment)^2 / evol.div_f)
    mortality = 1 - fitness * (1 -  ecol.d_tl[tl])

    push!(patch.species_ID[s], s)
    push!(patch.genotype[s], genotype)
    push!(patch.phenotype[s], phenotype)
    push!(patch.fitness[s], fitness)
    push!(patch.mortality[s], mortality)

    patch.N_s[s] += 1
    patch.N_tl[tl] += 1
    patch.total_N += 1
    if mortality > patch.mort_max_s[s]
        patch.mort_max_s[s] = mortality
    end
end

function push_individual!(patch::Patch, ecol::Ecol_parameters, evol::Evol_parameters, s::Int)
    g = evol.tot_genes
    e = patch.environment
    v = sign(e)*cld(round(abs(e)*g),g)
    n = abs(round(e*g/v))
    genotype = Int8.([(i > cld(n,2) && i <= g/2) || (i > g/2 + fld(n,2)) ? 0 : v for i in 1:g])
    # genotype = Int8.(round.(rand(evol.tot_genes) .* 0.5 .* rand([-1.0, 1.0], evol.tot_genes) .+ patch.environment))
    push_individual!(patch, ecol, evol, s, genotype)
end

function push_individual!(patch::Patch, parent_patch::Patch, ecol::Ecol_parameters, evol::Evol_parameters, s::Int, mother)
    # genotype = [i for i in parent_patch.genotype[s][mother]]
    genotype = copy(parent_patch.genotype[s][mother])

    # muts = rand(evol.tot_genes) .<= evol.mu
    # muts_sum = sum(muts)
    # if muts_sum > 0
    #     genotype[muts] .+= rand([-1,1] , muts_sum)
    # end

    k = rand(Binomial(evol.tot_genes, evol.mu))
    if k > 0
        muts = shuffle(1:evol.tot_genes)[1:k]
        genotype[muts] .+= rand([-1,1] , k)
    end

    push_individual!(patch, ecol, evol, s, genotype)
end

function push_individual!(patch::Patch, parent_patch::Patch, ecol::Ecol_parameters, evol::Evol_parameters, s::Int, mother, father)
    genotype = zeros(Int8, evol.tot_genes)

    k_genes_m = rand(Binomial(evol.trait_loci, 0.5))
    genes_sh = randperm(evol.trait_loci)
    # genes_m = genes_sh[1:k_genes_m]
    # genes_f = genes_sh[(k_genes_m+1):evol.trait_loci]
    for i in (@view genes_sh[1:k_genes_m])
        genotype[evol.all_mother[i]] = parent_patch.genotype[s][mother][evol.all_mother[i]]
    end
    for i in (@view genes_sh[(k_genes_m+1):evol.trait_loci])
        genotype[evol.all_mother[i]] = parent_patch.genotype[s][mother][evol.all_father[i]]
    end
    # genes = rand(evol.trait_loci) .< 0.5
    # genotype[evol.all_mother[genes]] .= parent_patch.genotype[s][mother][evol.all_mother[genes]]
    # genotype[evol.all_mother[.!genes]] .= parent_patch.genotype[s][mother][evol.all_father[.!genes]]
    k_genes_m = rand(Binomial(evol.trait_loci, 0.5))
    genes_sh = randperm(evol.trait_loci)
    # genes_m = genes_sh[1:k_genes_m]
    # genes_f = genes_sh[(k_genes_m+1):evol.trait_loci]
    for i in (@view genes_sh[1:k_genes_m])
        genotype[evol.all_father[i]] = parent_patch.genotype[s][father][evol.all_mother[i]]
    end
    for i in (@view genes_sh[(k_genes_m+1):evol.trait_loci])
        genotype[evol.all_father[i]] = parent_patch.genotype[s][father][evol.all_father[i]]
    end
    # genes = rand(evol.trait_loci) .< 0.5
    # genotype[evol.all_father[genes]] .= parent_patch.genotype[s][father][evol.all_mother[genes]]
    # genotype[evol.all_father[.!genes]] .= parent_patch.genotype[s][father][evol.all_father[.!genes]]

        # muts = rand(evol.tot_genes) .<= evol.mu
        # muts_sum = sum(muts)
        # if muts_sum > 0
        #     genotype[muts] .+= rand([-1,1] , muts_sum)
        # end

        k = rand(Binomial(evol.tot_genes, evol.mu))
        if k > 0
            muts = shuffle(1:evol.tot_genes)[1:k]
            genotype[muts] .+= rand([-1,1] , k)
        end

    push_individual!(patch, ecol, evol, s, genotype)
end

function populate_patch(patch::Patch, ecol::Ecol_parameters, evol::Evol_parameters, N, spec_dist::Symbol)
    # if patch.X == 6
    if spec_dist == :X
        prob_species = [(((s - 1) ÷ (ecol.species / ecol.grid.X) == (patch.X - 1) ? 1 : 0) /
            ecol.bodymass_tl[ecol.tl_species[s]]) ^ (1) for s in 1:ecol.species] # probabilities for different species per X
    else
        prob_species = [(((s - 1) ÷ (ecol.species / ecol.patches) == (patch.patch_ID - 1) ? 1 : 0) /
            ecol.bodymass_tl[ecol.tl_species[s]]) ^ (1) for s in 1:ecol.species] # probabilities for different species per patch
        # prob_species = [1 /
        #     ecol.bodymass_tl[ecol.tl_species[s]] for s in 1:ecol.species]
    end
    s_id = wsample(1:ecol.species, prob_species, N)
    for i in 1:N
        push_individual!(patch, ecol, evol, s_id[i])
    end
    # end
end

function remove_individual!(patch::Patch, ecol::Ecol_parameters, s::Int, i)
    tl = ecol.tl_species[s]
    patch.N_tl[tl] -= 1
    mort_remove = patch.mortality[s][i]
    species_ID = pop!(patch.species_ID[s])
    genotype = pop!(patch.genotype[s])
    phenotype = pop!(patch.phenotype[s])
    fitness = pop!(patch.fitness[s])
    mortality = pop!(patch.mortality[s])
    if i < patch.N_s[s]
        patch.species_ID[s][i] = species_ID
        patch.genotype[s][i] = genotype
        patch.phenotype[s][i] = phenotype
        patch.fitness[s][i] = fitness
        patch.mortality[s][i] = mortality
    end
    patch.N_s[s] -= 1
    patch.total_N -= 1
    if patch.N_s[s] == 0
        patch.mort_max_s[s] = 0
    else
        if mort_remove == patch.mort_max_s[s]
            patch.mort_max_s[s] = maximum(patch.mortality[s])
        end
    end
end

function update_uptake(patch::Patch, ecol::Ecol_parameters, dm::Direct_method)
    match_N = ecol.match_matr .* patch.N_s'
    patch.N_prey = vec(sum(match_N; dims=2))
    @. patch.tl_ratio = ecol.species_per_tl/patch.N_tl
    patch.uptake = zeros(ecol.species)
    for s in 1:ecol.species
        tl = ecol.tl_species[s]
        if tl == 1
            upt = ecol.uptake_tl[tl, 1] * patch.resource^ecol.uptake_pars[3]
        else
            upt = ecol.uptake_tl[tl, 1] * patch.N_prey[s]^ecol.uptake_pars[3]
        end
        patch.uptake[s] = upt/(1. + upt*ecol.uptake_tl[tl, 2])
        patch.gain_s[s] = patch.uptake[s] * ecol.conversion_tl[tl]
    end
    patch.total_uptake = patch.uptake .* patch.N_s
    match_loss = patch.total_uptake .* ecol.match_rel
    patch.loss_s .= vec(sum(match_loss; dims = 1)) .* patch.tl_ratio[ecol.tl_species]
    resource_uptake = sum(patch.total_uptake .* ecol.match_resource)
    patch.resource_gain = ecol.in_rate
    patch.resource_loss = ecol.out_rate*patch.resource + resource_uptake
    update_dm_patch(dm, ecol, patch)
end

function update_uptake(patch::Patch, ecol::Ecol_parameters, dm::Direct_method, s, dir=:UP)
    if s == 0
        tl = 1
        s_comp = (1:ecol.species)[ecol.tl_species .== tl]
        s_comp_nonz = s_comp[patch.N_s[s_comp] .> 0]
        resource_uptake = 0.
        upt = ecol.uptake_tl[tl, 1] * patch.resource^ecol.uptake_pars[3]
        upt = upt/(1. + upt*ecol.uptake_tl[tl, 2])
        for s2 in s_comp_nonz 
            patch.uptake[s2] = upt
            patch.gain_s[s2] = patch.uptake[s2] * ecol.conversion_tl[tl]
            patch.total_uptake[s2] = patch.uptake[s2] .* patch.N_s[s2]
            resource_uptake += patch.total_uptake[s2]
        end
        # patch.resource_gain = ecol.in_rate
        patch.resource_loss = ecol.out_rate*patch.resource + resource_uptake
    elseif s > 0
        tl_comp = ecol.tl_species[s]
        tl_prey = tl_comp - 1
        tl_pred = tl_comp + 1
        s_prey = (1:ecol.species)[ecol.tl_species .== tl_prey]
        s_prey_nonz = s_prey[patch.N_s[s_prey] .> 0]
        s_comp = (1:ecol.species)[ecol.tl_species .== tl_comp]
        s_comp_nonz = s_comp[patch.N_s[s_comp] .> 0]
        s_pred = (1:ecol.species)[ecol.tl_species .== tl_pred]
        s_pred_nonz = s_pred[patch.N_s[s_pred] .> 0]
        patch.tl_ratio[tl_comp] = ecol.species_per_tl/patch.N_tl[tl_comp]
        if tl_comp == 1
            if dir == :UP
                if patch.N_s[s] == 1
                    upt = ecol.uptake_tl[tl_comp, 1] * patch.resource^ecol.uptake_pars[3]
                    patch.uptake[s] = upt/(1. + upt*ecol.uptake_tl[tl_comp, 2])
                    patch.gain_s[s] = patch.uptake[s] * ecol.conversion_tl[tl_comp]
                end
                patch.resource_loss += patch.uptake[s]
                patch.N_prey[s_pred] .+= ecol.match_matr[s_pred, s]
                patch.total_uptake[s] += patch.uptake[s]
            else
                patch.resource_loss -= patch.uptake[s]
                patch.N_prey[s_pred] .-= ecol.match_matr[s_pred, s]
                patch.total_uptake[s] -= patch.uptake[s]
            end
            for s2 in s_pred_nonz 
                upt = ecol.uptake_tl[tl_pred, 1] * patch.N_prey[s2]^ecol.uptake_pars[3]
                patch.uptake[s2] = upt/(1. + upt*ecol.uptake_tl[tl_pred, 2])
                patch.gain_s[s2] = patch.uptake[s2] * ecol.conversion_tl[tl_pred]
                patch.total_uptake[s2] = patch.uptake[s2] * patch.N_s[s2]
            end
            match_loss = patch.total_uptake[s_pred_nonz] .* ecol.match_rel[s_pred_nonz, s_comp_nonz]
            patch.loss_s[s_comp_nonz] .= vec(sum(match_loss; dims = 1)) .* patch.tl_ratio[tl_comp]
        elseif tl_comp > 1 && tl_comp < ecol.trophic_levels
            if dir == :UP
                if patch.N_s[s] == 1
                    upt = ecol.uptake_tl[tl_comp, 1] * patch.N_prey[s]^ecol.uptake_pars[3]
                    patch.uptake[s] = upt/(1. + upt*ecol.uptake_tl[tl_comp, 2])
                    patch.gain_s[s] = patch.uptake[s] * ecol.conversion_tl[tl_comp]
                end
                patch.N_prey[s_pred] .+= ecol.match_matr[s_pred, s]
            else
                patch.N_prey[s_pred] .-= ecol.match_matr[s_pred, s]
            end
            for s2 in s_pred_nonz 
                upt = ecol.uptake_tl[tl_pred, 1] * patch.N_prey[s2]^ecol.uptake_pars[3]
                patch.uptake[s2] = upt/(1. + upt*ecol.uptake_tl[tl_pred, 2])
                patch.gain_s[s2] = patch.uptake[s2] * ecol.conversion_tl[tl_pred]
                patch.total_uptake[s2] = patch.uptake[s2] * patch.N_s[s2]
            end
            match_loss = patch.total_uptake[s_pred_nonz] .* ecol.match_rel[s_pred_nonz, s_comp_nonz]
            patch.loss_s[s_comp_nonz] .= vec(sum(match_loss; dims = 1)) .* patch.tl_ratio[tl_comp]
            @. patch.loss_s[s_prey_nonz] -= patch.total_uptake[s] * ecol.match_rel[s, s_prey_nonz] * patch.tl_ratio[tl_prey]
            patch.total_uptake[s] = patch.uptake[s] * patch.N_s[s]
            @. patch.loss_s[s_prey_nonz] += patch.total_uptake[s] * ecol.match_rel[s, s_prey_nonz] * patch.tl_ratio[tl_prey]
        else
            if dir == :UP && patch.N_s[s] == 1
                upt = ecol.uptake_tl[tl_comp, 1] * patch.N_prey[s]^ecol.uptake_pars[3]
                patch.uptake[s] = upt/(1. + upt*ecol.uptake_tl[tl_comp, 2])
                patch.gain_s[s] = patch.uptake[s] * ecol.conversion_tl[tl_comp]
            end
            patch.loss_s[s] = 0.
            @. patch.loss_s[s_prey_nonz] -= patch.total_uptake[s] * ecol.match_rel[s, s_prey_nonz] * patch.tl_ratio[tl_prey]
            patch.total_uptake[s] = patch.uptake[s] * patch.N_s[s]
            @. patch.loss_s[s_prey_nonz] += patch.total_uptake[s] * ecol.match_rel[s, s_prey_nonz] * patch.tl_ratio[tl_prey]
        end
        if dir == :DOWN && patch.N_s[s] == 0
            patch.uptake[s] = patch.gain_s[s] = patch.total_uptake[s] = patch.loss_s[s] = 0.
        end
    end
    update_dm_patch(dm, ecol, patch)
 end

function adjust_to_range(val::Float64, min::Float64, max::Float64)
    range = max - min
    quot = fld(val - max, range)
    rem = mod((val - max), range)
    minAdd = mod(quot, 2)
    maxAdd = 1 - minAdd
    return minAdd * (min + rem) + maxAdd * (max - rem)
end

function change_environment_local(patch::Patch, ecol::Ecol_parameters, evol::Evol_parameters, dm::Direct_method, env_X)
    new_env = patch.environment + ecol.env_step_local * ecol.dt_env * rand([-1.,1.])
    patch.environment = adjust_to_range(new_env, env_X + ecol.env_range.Y[1], env_X + ecol.env_range.Y[2])
    for s in 1:ecol.species
        if patch.N_s[s] > 0
            tl = ecol.tl_species[s]
            for i in 1:patch.N_s[s]
                patch.fitness[s][i] = exp(-(patch.phenotype[s][i]-patch.environment)^2 / evol.div_f)
                patch.mortality[s][i] = 1 - patch.fitness[s][i] * (1 -  ecol.d_tl[tl])
            end
            patch.mort_max_s[s] = maximum(patch.mortality[s])
        end
    end
    update_dm_loss(dm, ecol, patch)
end


mutable struct World
    # patches
    nbr_patches::Int
    patch_XY::Array{Int, 2}
    env_X::Vector{Float64}
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
        range_X = ecol.env_range.X[2] - ecol.env_range.X[1]
        range_Y = ecol.env_range.Y[2] - ecol.env_range.Y[1]
        step_X =
            grid.X == 1 ? 0.0 : range_X / (grid.X - 1)
        # step_Y =
        #     grid.Y == 1 ? 0.0 : (ecol.env_range.Y[2] - ecol.env_range.Y[1]) / (grid.Y - 1)
        # init_environment = [
        #     ecol.env_range.X[1] +
        #     step_X * (patch_XY[1, p] - 1) +
        #     ecol.env_range.Y[1] +
        #     step_Y * (patch_XY[2, p] - 1) for p = 1:nbr_patches
        # ]

        env_X = ecol.env_range.X[1] .+ step_X .* ((1:grid.X) .- 1)
        init_environment = [env_X[patch_XY[1, p]] + ecol.env_range.Y[1] + rand() * range_Y for p = 1:nbr_patches]

        patches = [Patch(ecol, p, patch_XY[1, p], patch_XY[2, p], init.resource, init_environment[p]) for p in 1:nbr_patches];
        for p in 1:nbr_patches
            populate_patch(patches[p], ecol, evol, Int(init.N*ecol.scale_uptake), init.spec_dist)
        end
        return new(nbr_patches, patch_XY, env_X, patches, neighbours, m_neighbours, cs_m_neighbours)
    end
end

function calc_m_neighbours(world::World, ecol::Ecol_parameters)
    temp_m_neighbours = @. ecol.rho * exp(-ecol.rho * world.neighbours)
    temp_m_neighbours[diagind(temp_m_neighbours)] .= 0.
    temp_m_neighbours ./= sum(temp_m_neighbours; dims = 1)
    world.m_neighbours = repeat(temp_m_neighbours, outer = [1, 1, ecol.trophic_levels])
    disp_tl = ecol.m .* ecol.bodymass_tl .^ ecol.m_power
    world.m_neighbours .= [
        i == j ? 1 - disp_tl[c] : disp_tl[c] * world.m_neighbours[i, j, c]
        for i = 1:world.nbr_patches, j = 1:world.nbr_patches, c = 1:ecol.trophic_levels
    ]
    world.cs_m_neighbours = cumsum(world.m_neighbours; dims = 1)
end

function get_patch_X(world::World, p::Int)
    return world.patch_XY[1, p]
end

function get_patch_Y(world::World, p::Int)
    return world.patch_XY[2, p]
end

function change_environment_CC(world::World, ecol::Ecol_parameters)
    world.env_X .+= ecol.env_step_CC * ecol.dt_env
end

function change_environment_local(world::World, ecol::Ecol_parameters, evol::Evol_parameters, dm::Direct_method)
    for p in 1:world.nbr_patches
        x = world.patch_XY[1, p]
        change_environment_local(world.patches[p], ecol, evol, dm, world.env_X[x])
    end
end

function next_event(world::World, ecol::Ecol_parameters, evol::Evol_parameters, dm::Direct_method)
    dt = 0.
    smpl_event = :nothing
    event = false
    p = 0
    new_p = 0
    s = 0
    tl = 0
    dir = :UP
    while !event
        dt += next_time(dm)
        smpl_event, s, p = sample_c(dm, ecol)

        # println("next_event1: smpl_event = $smpl_event; s = $s; p = $p")

        patch = world.patches[p]
        if smpl_event == :resource_gain
            event = increase_resource(patch, dm)
            dir = :UP
        elseif smpl_event == :resource_loss
            event = decrease_resource(patch, dm)
            dir = :DOWN
        elseif smpl_event == :reproduction
            tl = ecol.tl_species[s]
            # event = rand() <= (patch.gain_s[s]/dm.c_b_ind[s, p])
            event = true
            dir = :UP
            if event
                mother = rand(1:patch.N_s[s])
                new_p = binary_search((@view world.cs_m_neighbours[:, p, tl]), rand())
                new_patch = world.patches[new_p]
                if ecol.rep_type == :ASEXUAL
                    push_individual!(new_patch, patch, ecol, evol, s, mother)
                else
                    father = rand(1:patch.N_s[s])
                    push_individual!(new_patch, patch, ecol, evol, s, mother, father)
                end
            end
        else
            dir = :DOWN
            tl = ecol.tl_species[s]
            i = rand(1:patch.N_s[s])
            event = rand() <= ((patch.mortality[s][i] + patch.loss_s[s])/dm.c_d_ind[s, p])
            if event
                remove_individual!(patch, ecol, s, i)
            end
        end
    end
    if event

        # println("next_event2: smpl_event = $smpl_event; s = $s; p = $p; new_p = $new_p")

        if new_p > 0
            # update_uptake(world.patches[new_p], ecol, dm)
            update_uptake(world.patches[new_p], ecol, dm, s, dir)
        else
            # update_uptake(world.patches[p], ecol, dm)
            update_uptake(world.patches[p], ecol, dm, s, dir)
        end
    end
    return dt, smpl_event
end

function total_N(world::World)
    N = sum(map(p -> p.total_N, world.patches))
    return N
end

function genotype_mean(patch::Patch, s)
    mn = 0.
    if patch.N_s[s] > 0
        mn = mean(map(mean, patch.genotype[s]))
    else
        mn = NaN
    end
    return mn
end

function genotype_var(patch::Patch, s)
    vr = 0.
    if patch.N_s[s] > 1
        vr = var(map(mean, patch.genotype[s]))
    else
        vr = NaN
    end
    return vr
end

function phenotype_mean(patch::Patch, s)
    mn = 0.
    if patch.N_s[s] > 0
        mn = mean(patch.phenotype[s])
    else
        mn = NaN
    end
    return mn
end

function phenotype_var(patch::Patch, s)
    vr = 0.
    if patch.N_s[s] > 1
        vr = var(patch.phenotype[s])
    else
        vr = NaN
    end
    return vr
end

function fitness_mean(patch::Patch, s)
    mn = 0.
    if patch.N_s[s] > 0
        mn = mean(patch.fitness[s])
    else
        mn = NaN
    end
    return mn
end

function fitness_var(patch::Patch, s)
    vr = 0.
    if patch.N_s[s] > 1
        vr = var(patch.fitness[s])
    else
        vr = NaN
    end
    return vr
end

function log_titles(f)
    write(f, 
    "grid_X;grid_Y;torus_X;torus_Y;patches;nbr_species;match_sd;m;rho;" * 
    "CC;time_CC;pre_CC;post_CC;e_step_local;" * 
    "nbr_loci;sigma_z;mu;omega_e;d;rep_type;" * 
    "run;time;patch;X;Y;environment;resource;" * 
    "species;match;trophic_level;bodymass;mortality;N;biomass;" * 
    "genotype_mean;genotype_var;phenotype_mean;phenotype_var;fitness_mean;fitness_var\n")
end

function log_results(f, world::World, ecol::Ecol_parameters, evol::Evol_parameters, run::Run, r, t)
    for p in 1:world.nbr_patches, s in 1:ecol.species
        patch = world.patches[p]
        tl = ecol.tl_species[s]
        mch = ecol.species_match[s]
        write(f, 
        "$(ecol.grid.X);$(ecol.grid.Y);$(ecol.torus.X);$(ecol.torus.Y);$(ecol.patches);$(ecol.species);$(ecol.match_sd);$(ecol.m);$(ecol.rho);" *
        "$(ecol.CC);$(ecol.time_CC);$(run.pre_change);$(run.post_change);$(ecol.env_step_local);" * 
        "$(evol.trait_loci);$(evol.sigma_z);$(evol.mu);$(evol.omega_e);$(ecol.d);$(ecol.rep_type);" * 
        "$(r);$(t);$(p);$(patch.X);$(patch.Y);$(round(patch.environment, digits=2));$(patch.resource);" * 
        "$(s);$(mch);$(tl);$(ecol.bodymass_tl[tl]);$(round(ecol.d_tl[tl], digits=3));$(patch.N_s[s]);$(ecol.bodymass_tl[tl]*patch.N_s[s]);" * 
        "$(round(genotype_mean(patch, s), digits=2));$(round(genotype_var(patch, s), digits=5));$(round(phenotype_mean(patch, s), digits=2));$(round(phenotype_var(patch, s), digits=5));$(round(fitness_mean(patch, s), digits=2));$(round(fitness_var(patch, s), digits=5))\n")
    end
end

function print_results(world::World, r, t, m, cc, tcc, msd)
    println("run = $r;  time = $t; m = $m; CC = $cc; time_CC = $tcc; match_sd = $msd;  N = $(total_N(world))")
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


function evolving_foodweb_dm(init::Init_values)
    # run = Run(init);
    f = open(init.output_file,"w")
    log_titles(f)

    for r in 1:init.runs
        for m in init.m_vec, cc in init.CC_vec, tcc in init.time_CC_vec, msd in init.match_sd_vec
            init.m = m
            init.match_sd = msd
            init.CC = cc
            init.time_CC = tcc
            init.env_step_CC = tcc > 0 ? cc/tcc : 0.
            run = Run(init);
            ecol = Ecol_parameters(init);
            evol = Evol_parameters(init);
            dm = Direct_method(ecol; c_b_res = ecol.in_rate, c_d_res = ecol.in_rate*2, c_b = 1., c_d = 1.);
            world = World(ecol, evol, init);
            for p in 1:world.nbr_patches
                update_uptake(world.patches[p], ecol, dm)
            end
            time = 0.
            time_env = ecol.dt_env
            time_log = run.log_steps
            time_print = run.print_steps
            events = zeros(Int, 4)

            # m_min = 0.0
            # m_max = 0.2
            # m_step = 0.01
            # m_sign = 1.
            # m_timestep = 1000
            # ecol.m = m_min
            # time_m = m_timestep

            # println("init step_CC = $(init.env_step_CC)")
            # println("ecol step_CC = $(ecol.env_step_CC)")
            # println("dt_env = $(ecol.dt_env)")

            log_results(f, world, ecol, evol, run, r, round(time))
            println()
            print_results(world, r, round(time), m, cc, tcc, msd)

            while time < run.time_steps
                dt, event = next_event(world, ecol, evol, dm)
                time += dt
                if event == :resource_gain
                    events[1] += 1
                elseif event == :resource_loss
                    events[2] += 1
                elseif event == :reproduction
                    events[3] += 1
                else
                    events[4] += 1
                end
                if time > time_env
                    if ecol.env_step_CC > 0. && time > run.pre_change && time < (run.time_steps - run.post_change)
                        change_environment_CC(world, ecol)
                    end
                    change_environment_local(world, ecol, evol, dm)
                    time_env += ecol.dt_env
                end

                # if time > time_m
                #     if ecol.m < m_min || ecol.m > m_max
                #         m_sign = -m_sign
                #     end
                #     ecol.m += m_sign * m_step
                #     calc_m_neighbours(world, ecol)
                #     time_m += m_timestep
                # end

                if time > time_log
                    log_results(f, world, ecol, evol, run, r, round(time))
                    time_log += run.log_steps
                end
                if time > time_print
                    print_results(world, r, round(time), m, cc, tcc, msd)
                    time_print += run.print_steps
                end
            end
        end
    end
    close(f)
    # return time, events, world, ecol, dm
end

