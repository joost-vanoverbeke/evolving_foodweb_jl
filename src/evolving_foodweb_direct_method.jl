

using Random
using StatsBase
using LinearAlgebra
using Distributions

struct Ecol_parameters
    in_rate::Float64
    out_rate::Float64
    d::Float64
    d_power::Float64
    patches::Int
    species::Int
    trophic_levels::Int
    bodymass_tl::Array{Float64, 1}
    tl_species::Vector{UInt}
    uptake_pars::Vector{Float64}
    i_power::Float64
    uptake_tl::Array{Float64, 2}
    resource_conversion::Float64
    resource_assimilation::Float64
    assimilation_eff::Float64
    scale_uptake::Float64
    scale_assim::Float64

    function Ecol_parameters(; in_rate = 200., out_rate = 0.1, d = 0.1, d_power = -0.25,
        patches = 1, species = 3, trophic_levels = 3, uptake_pars = [0.0001, 0.4, 1.], i_power = 0.75,
        resource_conversion = 1., resource_assimilation = 10., assimilation_eff = 0.7,
        scale_uptake = 1., scale_assim = 0.)
        in_rate = in_rate*scale_uptake
        uptake_pars[1] /= scale_uptake
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
        return new(in_rate, out_rate, d, d_power,
            patches, species, trophic_levels, bodymass_tl, tl_species, uptake_pars, i_power, uptake_tl, resource_conversion, resource_assimilation, assimilation_eff)
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
        return new(omega_e,  trait_loci, tot_genes, mu, sigma_z, div_f, all_mother, all_father, all_genes)
    end
end


mutable struct Direct_method
    c_b_resource::Float64
    c_d_resource::Float64
    # c_b_ind::Vector{Float64}
    # c_d_ind::Vector{Float64}
    c_b_ind::Float64
    c_d_ind::Float64
    c_vec::Vector{Float64}
    c_total::Float64

    function Direct_method(ecol::Ecol_parameters; c_b_resource = 200., c_d_resource = 5000., c_b = 2., c_d = 2.)
        # c_b_ind = [c_b/i for i in 1:ecol.trophic_levels]
        # c_d_ind = [c_d for i in 1:ecol.trophic_levels]
        c_b_ind = c_b
        c_d_ind = c_d
        c_vec = [c_b_resource; c_d_resource; c_b_ind; c_d_ind]
        c_total = sum(c_vec)
        return new(c_b_resource, c_d_resource, c_b_ind, c_d_ind, c_vec, c_total)
    end
end

function update_dm(dm::Direct_method, resource, uptake, N)
    dm.c_b_resource = resource.in_rate
    dm.c_d_resource = resource.out_rate*resource.level + uptake[1]
    # dm.c_vec .= [dm.c_b_resource; dm.c_d_resource; dm.c_b_ind.*N_s; dm.c_d_ind.*N_s]
    dm.c_vec[1] = dm.c_b_resource
    dm.c_vec[2] = dm.c_d_resource
    dm.c_vec[3] = dm.c_b_ind*N
    dm.c_vec[4] = dm.c_d_ind*N
    dm.c_total = sum(dm.c_vec)
end

function next_time(dm::Direct_method)
    t_next = rand(Exponential(1. / dm.c_total))
    return t_next
end

function sample_c(dm::Direct_method)
    c_sel = wsample(dm.c_vec)
    return c_sel
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


# mutable struct Unit
#     gain::Float64
#     loss::Float64
# end
#
# function update_gain(unit::Unit, in_rate, level)
#     unit.gain = in_rate/level
# end
#
# function update_loss(unit::Unit, out_rate, level, uptake::Vector{Float64})
#     unit.loss = out_rate + uptake[1]/level
# end


mutable struct Resource
    in_rate::Float64
    out_rate::Float64
    # units::Vector{Unit}
    level::Float64
    gain::Float64
    loss::Float64

    function Resource(ecol::Ecol_parameters, level = 2.)
        # tot::Int = Int(ecol.in_rate/ecol.out_rate)
        # units = [Unit(ecol.in_rate/level, ecol.out_rate) for i=1:tot]
        # return new(ecol.in_rate, ecol.out_rate, units, level, 0., 0.)
        return new(ecol.in_rate, ecol.out_rate, level, 0., 0.)
    end
end

function update_gain(resource::Resource)
    # resource.gain = accumulate(x,y -> x+y.gain, resource.units, init = 0)
    resource.gain = resource.in_rate
end

function update_loss(resource::Resource, uptake::Vector{Float64})
    # resource.loss = accumulate(x,y -> x+y.loss, resource.units, init = 0)
    resource.loss = resource.out_rate*resource.level + uptake[1]
end

function increase_level(resource::Resource, dm::Direct_method)
    # i = rand(1:resource.level)
    # event = rand() <= (resource.units[i].gain/dm.c_b_resource)
    # if event
    #     resource.level += 1
    #     units[level] = Unit(in_rate/level, out_rate + uptake[1]/level)
    # end
    event = rand() <= (resource.gain/dm.c_b_resource)
    if event
        resource.level += 1
    end
    return event
end

function decrease_level(resource::Resource, dm::Direct_method)
    # i = rand(1:resource.level)
    # event = rand() <= (resource.units[i].loss/dm.c_d_resource)
    # if event
    #     units[i] = units[level]
    #     resource.level -= 1
    # end
    event = rand() <= (resource.loss/dm.c_d_resource)
    if event
        resource.level -= 1
    end
    return event
end


mutable struct Individual
    species_ID::Int
    trophic_level::Int
    bodymass::Float64
    genotype::Vector{Int8}
    phenotype::Float64
    fitness::Float64
    mortality_rate::Float64
    uptake_rate::Vector{Float64}
    conversion_rate::Float64
    uptake::Float64
    assimilation::Float64
    gain::Float64
    loss::Float64

    function Individual(ecol::Ecol_parameters, evol::Evol_parameters, e::Float64, s::Int, gtp::Vector{Int8})
        species_ID = s
        trophic_level = ecol.tl_species[s]
        bodymass = 10.0^(trophic_level - 1)
        genotype = gtp
        phenotype = mean(genotype) + randn() * evol.sigma_z
        fitness = exp(-(phenotype - e)^2 / evol.div_f)
        mortality_rate = 1 - fitness * (1 -  ecol.d * bodymass^ecol.d_power)
        if trophic_level == 1
            uptake_rate =
            [ecol.uptake_pars[1] * bodymass^ecol.i_power,
            ecol.uptake_pars[2] * bodymass^(-ecol.i_power) * ecol.resource_conversion,
            ecol.uptake_pars[3]]
            conversion_rate = ecol.resource_assimilation/bodymass
        else
            uptake_rate =
            [ecol.uptake_pars[1] * bodymass^ecol.i_power,
            ecol.uptake_pars[2] * bodymass^(-ecol.i_power) * 10.0^(trophic_level - 2),
            ecol.uptake_pars[3]]
            conversion_rate = ecol.assimilation_eff*(1-ecol.scale_assim*bodymass^(-ecol.i_power)) * 10.0^(trophic_level - 2)/bodymass
        end
        uptake = 0.
        assimilation = 0.
        gain = 0.
        loss = 0.
        return new(species_ID, trophic_level, bodymass, genotype, phenotype, fitness, mortality_rate, uptake_rate, conversion_rate, uptake, assimilation, gain, loss)
    end

    function Individual(ecol::Ecol_parameters, evol::Evol_parameters, e = 0.; s = 1)
        genotype = Int8.(round.(rand(evol.tot_genes) .* 0.5 .* rand([-1.0, 1.0], evol.tot_genes) .+ e))
        return Individual(ecol, evol, e, s, genotype)
    end

    function Individual(ecol::Ecol_parameters, evol::Evol_parameters, mother::Individual, e = 0.)
        species_ID = mother.species_ID
        genotype = deepcopy(mother.genotype)
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

function update_fitness(ind::Individual, env::Float64, evol::Evol_parameters)
    ind.fitness =
        exp(-(ind.phenotype-env)^2 / evol.div_f)
end

function update_uptake(ind::Individual, resource, N)
    if ind.trophic_level == 1
        upt = ind.uptake_rate[1]*resource^ind.uptake_rate[3]
    else
        upt = ind.uptake_rate[1]*N[ind.trophic_level-1]^ind.uptake_rate[3]
    end
    ind.uptake = upt/(1. + upt*ind.uptake_rate[2])
end

function update_assimilation(ind::Individual)
    ind.assimilation = ind.uptake*ind.conversion_rate
end

function update_gain(ind::Individual)
    ind.gain = ind.assimilation
end

function update_loss(ind::Individual, uptake, N)
    if ind.trophic_level < length(uptake)
        ind.loss = ind.mortality_rate + uptake[ind.trophic_level+1]/N[ind.trophic_level]
    else
        ind.loss = ind.mortality_rate
    end
end

function die(ind::Individual, dm::Direct_method)
    event = rand() <= (ind.loss/dm.c_d_ind)
    return event
end

function reproduce(ind::Individual, dm::Direct_method)
    event = rand() <= (ind.gain/dm.c_b_ind)
    return event
end


mutable struct Patch
    patch_ID::Int
    environment::Float64
    resource::Resource
    microsites::Int
    individuals::Vector{Individual}
    N::Vector{Int}
    total_N::Int
    # total_mortality::Vector{Float64}
    total_uptake::Vector{Float64}
    # total_assimilation::Vector{Float64}

    function Patch(ecol::Ecol_parameters, evol::Evol_parameters, id = 1, e = 0.; r::Float64, m::Int, init_N::Int)
        patch_ID = id
        environment = e
        resource = Resource(ecol, r)
        microsites = m
        individuals = [Individual(ecol, evol) for i=1:microsites]
        N = zeros(Int64, ecol.trophic_levels)
        # total_mortality = zeros(ecol.trophic_levels)
        total_uptake = zeros(ecol.trophic_levels)
        # total_assimilation = zeros(ecol.trophic_levels)
        prob_species = [((s - 1) ÷ (ecol.species / ecol.patches) == (patch_ID - 1) ? 1 : 0) /
            10^(1.5*(ecol.tl_species[s]-1)) for s = 1:ecol.species]
        # prob_species = [1,0,0]
            s_id = wsample(1:ecol.species, prob_species, init_N)
        for i in 1:init_N
            individuals[i] = Individual(ecol, evol, environment; s = s_id[i])
            N[individuals[i].trophic_level] += 1
            # total_mortality[individuals[i].trophic_level] += individuals[i].mortality_rate
        end
        total_N = sum(N)
        return new(patch_ID, environment, resource, microsites, individuals, N, total_N, total_uptake)
    end
end

function update_uptake(ecol::Ecol_parameters, patch::Patch, dm::Direct_method)
    for tl in 1:ecol.trophic_levels
        if tl == 1
            upt = ecol.uptake_tl[tl, 1]*patch.resource.level^ecol.uptake_pars[3]
        else
            upt = ecol.uptake_tl[tl, 1]*patch.N[tl-1]^ecol.uptake_pars[3]
        end
        patch.total_uptake[tl] = upt/(1. + upt*ecol.uptake_tl[tl, 2])*patch.N[tl]
    end
    update_dm(dm, patch.resource, patch.total_uptake, patch.total_N)
end

# function update_flows(patch::Patch, dm::Direct_method)
#     fill!(patch.total_uptake, 0.)
#     for i in 1:patch.total_N
#         update_uptake(patch.individuals[i], patch.resource.level, patch.N)
#         update_assimilation(patch.individuals[i])
#         patch.total_uptake[patch.individuals[i].trophic_level] += patch.individuals[i].uptake
#     end
#     update_gain(patch.resource)
#     update_loss(patch.resource, patch.total_uptake)
#     for i in 1:1:patch.total_N
#         update_gain(patch.individuals[i])
#         update_loss(patch.individuals[i], patch.total_uptake, patch.N)
#     end
#     # dm.c_b_resource = patch.resource.gain
#     # dm.c_d_resource = patch.resource.loss
# end

function update_fitness(patch::Patch, evol::Evol_parameters)
    for i in 1:sum(N)
        update_fitness(individuals[i], environment, evol)
    end
end

function change_environment(patch::Patch, ecol::Ecol_parameters, dt)
end

function next_event(patch::Patch, ecol::Ecol_parameters, evol::Evol_parameters, dm::Direct_method)
    dt = 0.
    event = false
    fails = zeros(Int, 4)
    maxs = zeros(4)
    smpl_c = 0
    while !event
        dt += next_time(dm)
        smpl_c = sample_c(dm)
        if smpl_c == 1
            update_gain(patch.resource)
            if patch.resource.gain > maxs[smpl_c]
                maxs[smpl_c] = patch.resource.gain
            end
            event = increase_level(patch.resource, dm)
        elseif smpl_c == 2
            update_loss(patch.resource, patch.total_uptake)
            if patch.resource.loss > maxs[smpl_c]
                maxs[smpl_c] = patch.resource.loss
            end
            event = decrease_level(patch.resource, dm)
        elseif smpl_c == 3
            i = rand(1:patch.total_N)
            update_uptake(patch.individuals[i], patch.resource.level, patch.N)
            update_assimilation(patch.individuals[i])
            update_gain(patch.individuals[i])
            if patch.individuals[i].gain > maxs[smpl_c]
                maxs[smpl_c] = patch.individuals[i].gain
            end
            if reproduce(patch.individuals[i], dm)
                patch.total_N += 1
                patch.individuals[patch.total_N] = Individual(ecol, evol, patch.individuals[i], patch.environment)
                patch.N[patch.individuals[patch.total_N].trophic_level] += 1
                event = true
            end
        else
            i = rand(1:sum(patch.total_N))
            update_loss(patch.individuals[i], patch.total_uptake, patch.N)
            if patch.individuals[i].loss > maxs[smpl_c]
                maxs[smpl_c] = patch.individuals[i].loss
            end
            if die(patch.individuals[i], dm)
                patch.N[patch.individuals[i].trophic_level] -= 1
                patch.individuals[i] = patch.individuals[patch.total_N]
                patch.total_N -= 1
                event = true
            end
        end
        if !event
            fails[smpl_c] += 1
        end
    end
    update_uptake(ecol, patch, dm)
    return dt, smpl_c, fails, maxs
end


function evolving_foodweb_dm(; time_steps = 1)
    ecol = Ecol_parameters(; in_rate = 200., scale_uptake = 5.)
    evol = Evol_parameters()
    dm = Direct_method(ecol; c_b_resource = 500., c_d_resource = 500., c_b = 0.7, c_d = 0.7)
    patch = Patch(ecol, evol; r = 200., m = 100000, init_N = 10000)
    update_uptake(ecol, patch, dm)
    # update_uptake(patch.individuals[1], patch.resource.level, patch.N)
    # patch.individuals[1].uptake*patch.N[1]
    time = 0
    events = zeros(Int, 4)
    fails = zeros(Int, 4)
    maxs = zeros(4)
    # update_flows(patch, dm)
    while time < time_steps
        # update_flows(patch, dm)
        dt, event, new_fails, new_maxs = next_event(patch, ecol, evol, dm)
        time += dt
        events[event] += 1
        fails .+= new_fails
        for i in 1:4
            if new_maxs[i] > maxs[i]
                maxs[i] = new_maxs[i]
            end
        end
    end
    return time, events, fails, maxs, patch, ecol
end

@time time, events, fails, maxs, patch, ecol = evolving_foodweb_dm(; time_steps = 250);
time
events
fails
@. events/(fails+events)
maxs

@profiler time, events, fails, maxs, patch, ecol = evolving_foodweb_dm(; time_steps = 100);


ind = Individual(ecol, evol, patch.individuals[1], patch.environment)



patch.individuals[patch.total_N].gain
patch.individuals[i].species_ID
patch.individuals[i].gain/dm.c_b_ind
patch.individuals[i].loss/dm.c_d_ind

next_time(dm, [20,39,13])
sample_c(dm, [7,9,3])
