
using Random
using StatsBase
using LinearAlgebra


struct Comm
    # patches
    grid::NamedTuple{(:X, :Y),Tuple{Int,Int}}
    torus::NamedTuple{(:X, :Y),Tuple{Symbol,Symbol}}
    env_range::NamedTuple{(:X, :Y),Tuple{Tuple{Float64, Float64},Tuple{Float64, Float64}}}
    env_step::Float64
    microsites::Int
    patches::Int
    tot_sites::Int
    patch_XY::Array{Int, 2}
    # species
    trophic_levels::Int
    species::Int
    rep_type::Symbol
    # bodymass classes
    bc_species::Array{Int, 1}
    bodymass_bc::Array{Float64, 1}
    # dispersal
    m::Float64
    rho::Float64
    m_bc::Symbol
    m_power::Float64
    neighbours::Array{Float64, 2}
    m_neighbours::Array{Float64, 3}
    # mortality
    d::Float64
    d_power::Float64
    mortality_bc::Array{Float64, 1}
    # resource
    in_rate::Float64
    out_rate::Float64
    resource_conversion::Float64
    # feeding
    uptake_pars::Array{Float64, 1}
    i_power::Float64
    uptake_bc::Array{Float64, 2}
    assimilation_eff::Float64
    scale_uptake::Float64

    function Comm(;
        grid = (X = 5, Y = 2),
        torus = (X = :NO, Y = :NO),
        env_range = (X = (-1., 1.), Y = (0., 0.)),
        env_step = 0.,
        microsites = 100,
        m = 0.1,
        rho = 2.,
        m_bc = :EQUAL,
        trophic_levels = 3,
        species = 30,
        rep_type = :ASEXUAL,
        d = 0.1,
        in_rate = 200.,
        dt = 0.1,
        scale_uptake = 1.
    )

        # patches
        patches = grid.X * grid.Y
        tot_sites = patches * microsites
        patch_XY = zeros(Int, 2, patches)
        for x in (1:grid.X) .- 1, y in (1:grid.Y) .- 1
            patch_XY[1, y*grid.X+x+1] = x + 1
            patch_XY[2, y*grid.X+x+1] = y + 1
        end
        # species
        # bodymass classes
        bc_species = [(s - 1) % trophic_levels + 1 for s = 1:species]
        bodymass_bc = [10^(c - 1) for c = 1:trophic_levels]
        # dispersal
        m_power = m_bc == :DECR ? -0.25 : (m_bc == :INCR ? 0.25 : 0.)
        neighbours = zeros(patches, patches)
        for x in (1:grid.X), y in (1:grid.Y), x2 in (1:grid.X), y2 in (1:grid.Y)
            d_X = abs(x - x2)
            d_X = (torus.X == :YES ? min(d_X, grid.X - d_X) : d_X)^2
            d_Y = abs(y - y2)
            d_Y = (torus.Y == :YES ? min(d_Y, grid.Y - d_Y) : d_Y)^2
            neighbours[(y-1)*grid.X+x, (y2-1)*grid.X+x2] = sqrt(d_X + d_Y)
        end
        temp_m_neighbours = @. rho * exp(-rho * neighbours)
        temp_m_neighbours[diagind(temp_m_neighbours)] .= 0.
        temp_m_neighbours ./= sum(temp_m_neighbours; dims = 1)
        m_neighbours = repeat(temp_m_neighbours, outer = [1, 1, trophic_levels])
        disp_bc = m .* bodymass_bc .^ m_power
        m_neighbours .= [
            i == j ? 1 - disp_bc[c] : disp_bc[c] * m_neighbours[i, j, c]
            for i = 1:patches, j = 1:patches, c = 1:trophic_levels
        ]
        # mortality
        d_power = -0.25
        mortality_bc = @. d * bodymass_bc^d_power
        # resource
        in_rate *= scale_uptake
        out_rate = 0.1
        resource_conversion = 10.0
        #  feeding
        uptake_pars = [1e-4, 0.4, 1.0]
        uptake_pars[1] /= scale_uptake
        i_power = 0.75
        uptake_bc = zeros(trophic_levels, 3)
        uptake_bc[:, 1] .= @. uptake_pars[1] * bodymass_bc^i_power
        uptake_bc[:, 2] .= @. uptake_pars[2] * bodymass_bc^(-i_power)
        uptake_bc[:, 3] .= uptake_pars[3]
        assimilation_eff = 0.7
        # dt
        in_rate *= dt
        out_rate *= dt
        d *= dt
        mortality_bc .*= dt
        env_step *= dt
        uptake_pars[1] *= dt
        uptake_pars[2] /= dt
        uptake_bc[:, 1] .*= dt
        uptake_bc[:, 2] ./= dt

        return new(
            grid,
            torus,
            env_range,
            env_step,
            microsites,
            patches,
            tot_sites,
            patch_XY,
            trophic_levels,
            species,
            rep_type,
            bc_species,
            bodymass_bc,
            m,
            rho,
            m_bc,
            m_power,
            neighbours,
            m_neighbours,
            d,
            d_power,
            mortality_bc,
            in_rate,
            out_rate,
            resource_conversion,
            uptake_pars,
            i_power,
            uptake_bc,
            assimilation_eff,
            scale_uptake
        )
    end
end

struct Evol
    omega_e::Float64
    trait_loci::Int
    tot_genes::Int
    mu::Float64
    sigma_z::Float64
    div_f::Float64
    all_mother::Array{Int, 1}
    all_father::Array{Int, 1}
    all_genes::Array{Int, 1}

    function Evol(; omega_e = 4.0, trait_loci = 20, mu = 1e-4, sigma_z = 0.1)
        div_f = 2 * omega_e^2
        tot_genes = 2 * trait_loci
        all_mother = 1:trait_loci
        all_father = all_mother .+ trait_loci
        all_genes = [all_mother; all_father]
        return new(
            omega_e,
            trait_loci,
            tot_genes,
            mu,
            sigma_z,
            div_f,
            all_mother,
            all_father,
            all_genes,
        )
    end
end

struct Sites
    comm::Comm
    evol::Evol
    environment::Array{Float64,2}
    species_ID::Array{UInt8,2}
    alive::BitArray{2}
    genotype::Array{Int8,3}
    phenotype::Array{Float64,2}
    fitness::Array{Float64,2}
    resource::Array{Float64,1}
    mass_abundance::Array{Int,2}
    uptake_prey::Array{Float64,2}
    consumption::Array{Float64,2}
    nbr_prey::Array{Int,2}
    prey_pos::Array{Int,3}
    nbr_empty::Array{Int,1}
    empty_pos::Array{Int,2}
    nbr_newborns::Array{Float64,1}
    prob_mothers::Array{Float64,2}
    prob_fathers::Array{Float64,3}

    function Sites(comm, evol)
        # environment
        step_X =
            comm.grid.X == 1 ? 0.0 : (comm.env_range.X[2] - comm.env_range.X[1]) / (comm.grid.X - 1)
        step_Y =
            comm.grid.Y == 1 ? 0.0 : (comm.env_range.Y[2] - comm.env_range.Y[1]) / (comm.grid.Y - 1)
        environment = [
            comm.env_range.X[1] +
            step_X * (comm.patch_XY[1, p] - 1) +
            comm.env_range.Y[1] +
            step_Y * (comm.patch_XY[2, p] - 1) for i in [1], p = 1:comm.patches
        ]

        # individuals
        species_ID = ones(UInt8, comm.microsites, comm.patches)
        alive = falses(comm.microsites, comm.patches)
        genotype =
            Int8.(round.(
                rand(evol.tot_genes, comm.microsites, comm.patches) .* 0.5 .*
                rand([-1.0, 1.0], evol.tot_genes, comm.microsites, comm.patches) .+
                reshape(environment, 1, 1, comm.patches),
            ))
        phenotype =
            reshape(mean(genotype; dims = 1), comm.microsites, comm.patches) .+
            randn(comm.microsites, comm.patches) .* evol.sigma_z
        fitness = @. exp(-(phenotype - environment)^2 / evol.div_f)

        # feeding
        resource = fill((comm.in_rate / comm.out_rate) / 2.0, comm.patches)
        mass_abundance = zeros(Int, comm.trophic_levels, comm.patches)
        uptake_prey = zeros(comm.trophic_levels, comm.patches)
        consumption = zeros(comm.trophic_levels, comm.patches)
        nbr_prey = zeros(Int, comm.trophic_levels, comm.patches)
        prey_pos = zeros(Int, comm.microsites, comm.trophic_levels, comm.patches)

        # mortality & reproduction
        nbr_empty = zeros(Int, comm.patches)
        empty_pos = zeros(Int, comm.microsites, comm.patches)
        nbr_newborns = zeros(comm.patches)
        prob_mothers = zeros(comm.tot_sites, comm.patches)
        prob_fathers = zeros(comm.microsites, comm.species, comm.patches)

        prob_species = [
            ((s - 1) ÷ (comm.species / comm.patches) == (p - 1) ? 1 : 0) /
            comm.bodymass_bc[comm.bc_species[s]] for s = 1:comm.species, p = 1:comm.patches
        ]
        # prob_species = [1,0,0]
        prob_species ./= sum(prob_species; dims = 1)
        init_N = Int(comm.microsites / 2)
        for p = 1:comm.patches
            init_inds = shuffle(1:comm.microsites)[1:init_N]
            species_ID[init_inds, p] .= wsample(1:comm.species, prob_species[:, p], init_N)
            alive[init_inds, p] .= true
            mass_abundance[:, p] .= addcounts!(
                mass_abundance[:, p],
                comm.bc_species[species_ID[init_inds, p]],
                1:comm.trophic_levels,
            )
            for i in init_inds
                bc = comm.bc_species[species_ID[i, p]]
                prey_pos[nbr_prey[bc, p]+1, bc, p] = i
                nbr_prey[bc, p] += 1
            end
        end

        return new(
            comm,
            evol,
            environment,
            species_ID,
            alive,
            genotype,
            phenotype,
            fitness,
            resource,
            mass_abundance,
            uptake_prey,
            consumption,
            nbr_prey,
            prey_pos,
            nbr_empty,
            empty_pos,
            nbr_newborns,
            prob_mothers,
            prob_fathers
        )
    end
end

struct Run
    dt::Float64
    runs::Int
    time_steps::Int
    pre_change::Int
    post_change::Int
    print_steps::Int
    save_steps::Int
    file::String

    function Run(;
        dt = 0.1,
        runs = 1,
        time_steps = 100,
        pre_change = 10,
        post_change = 10,
        print_steps = 10,
        save_steps = 100,
        file = "output_evolving_foodweb.csv",
    )
        return new(dt, runs, time_steps, pre_change, post_change, print_steps, save_steps, file)
    end
end


function new_body(sites::Sites, i, p, s)
    comm = sites.comm
    evol = sites.evol
    sites.alive[i, p] = true
    sites.species_ID[i, p] = s
    bc = comm.bc_species[s]
    bm = comm.bodymass_bc[bc]
    sites.mass_abundance[bc, p] += 1
    sites.phenotype[i, p] = mean(sites.genotype[:, i, p]) + randn() * evol.sigma_z
    sites.fitness[i, p] = exp(-(sites.phenotype[i, p] - sites.environment[p])^2 / evol.div_f)
    sites.prey_pos[sites.nbr_prey[bc, p]+1, bc, p] = i
    sites.nbr_prey[bc, p] += 1
end

function dead_body(sites::Sites, i, p)
    comm = sites.comm
    sites.alive[i, p] = false
    sites.mass_abundance[comm.bc_species[sites.species_ID[i, p]], p] -= 1
end

function change_environment(sites::Sites)
    comm = sites.comm
    sites.environment .+= comm.env_step
    adjust_fitness(sites)
end

function adjust_fitness(sites::Sites)
    evol = sites.evol
    sites.fitness[sites.alive] .=
        @. exp(-(sites.phenotype-sites.environment)[sites.alive]^2 / evol.div_f)
end

function uptake(sites::Sites)
    comm = sites.comm
    fill!(sites.consumption, 0.0)
    fill!(sites.uptake_prey, 0.0)
    abundance = Float64.(deepcopy(sites.mass_abundance))

    # uptake resource
    ing = @. comm.uptake_bc[1, 1] * sites.resource^comm.uptake_bc[1, 3]
    ingested = @. ing / (1 + ing * comm.uptake_bc[1, 2])
    sites.resource .+=
        @. (comm.in_rate - sites.resource * comm.out_rate) - ingested * abundance[1, :]
    sites.resource .= max.(0, sites.resource)
    sites.consumption[1, :] .+= @. comm.resource_conversion * ingested

    # uptake prey
    ing = @. comm.uptake_bc[2:comm.trophic_levels, 1] *
       (
        comm.bodymass_bc[1:comm.trophic_levels-1] .* abundance[1:comm.trophic_levels-1, :]
    )^comm.uptake_bc[2:comm.trophic_levels, 3]
    ingested = @. ing / (1 + ing * comm.uptake_bc[2:comm.trophic_levels, 2])
    abundance[1:comm.trophic_levels-1, :] .-= @. ingested * abundance[2:comm.trophic_levels, :] /
       comm.bodymass_bc[1:comm.trophic_levels-1]
    sites.consumption[2:comm.trophic_levels, :] .+=
        @. comm.assimilation_eff * ingested / comm.bodymass_bc[2:comm.trophic_levels]

    sites.uptake_prey .= @. min(sites.mass_abundance - abundance, sites.mass_abundance)
end

function update_prey(sites::Sites)
    comm = sites.comm
    for p = 1:comm.patches, l = 1:comm.trophic_levels
        # println("patch = ", p, "  tr lev =", l)
        pos_prey = sample(
            sites.prey_pos[1:sites.nbr_prey[l, p], l, p],
            round(Int, sites.uptake_prey[l, p]);
            replace = false,
        )
        map(i -> dead_body(sites, i, p), pos_prey)
    end
    fill!(sites.uptake_prey, 0)
end

function contribution_adults(sites::Sites)
    comm = sites.comm
    contr = 0.0
    fill!(sites.nbr_empty, 0)
    fill!(sites.nbr_prey, 0)
    # fill!(sites.prob_fathers, 0.0)
    # p = 1
    # i = 234
    for p = 1:comm.patches, i = 1:comm.microsites
        s = sites.species_ID[i, p]
        bc = comm.bc_species[s]
        if sites.alive[i, p] && rand() > (1 - comm.mortality_bc[bc]) * sites.fitness[i, p]
            dead_body(sites, i, p)
        end
        if sites.alive[i, p]
            sites.prey_pos[sites.nbr_prey[bc, p]+1, bc, p] = i
            sites.nbr_prey[bc, p] += 1
            contr = 1.0
        else
            sites.empty_pos[sites.nbr_empty[p]+1, p] = i
            sites.nbr_empty[p] += 1
            contr = 0.0
        end
        if comm.rep_type == :SEXUAL
            sites.prob_fathers[i, s, p] = contr
            # if i > 1
            #     sites.prob_fathers[i,:,p] .+= sites.prob_fathers[i-1,:,p]
            # end
        end
        contr *= sites.consumption[bc, p]
        i_long = i + comm.microsites * (p - 1)
        # p2 = 1
        for p2 in 1:comm.patches
            sites.prob_mothers[i_long, p2] = contr * comm.m_neighbours[p2, p, bc]
            if i_long > 1
                sites.prob_mothers[i_long, p2] += sites.prob_mothers[i_long-1, p2]
            end
            # if i_long == 1
            #     sites.prob_mothers[i_long, p2] = contr * comm.m_neighbours[p2, p, bc]
            # else
            #     sites.prob_mothers[i_long, p2] = sites.prob_mothers[i_long-1, p2] + contr * comm.m_neighbours[p2, p, bc]
            # end
        end
    end
    # sites.nbr_newborns .= reduce(+, sites.prob_mothers; dims = 1)[:]
    sites.nbr_newborns .= sites.prob_mothers[end, :]
    # sites.prob_mothers ./= sites.prob_mothers[end, :]'
    # for p in 1:comm.patches, i in 1:comm.tot_sites
    #     sites.prob_mothers[i,p] /= sites.prob_mothers[end, p]
    # end
    if comm.rep_type == :SEXUAL
        sites.prob_fathers .= cumsum(sites.prob_mothers; dims = 1)
        sites.prob_fathers ./= reshape(sites.prob_fathers[end, :, :], 1, comm.species, comm.patches)
    end
end

function reproduction(sites::Sites)
    comm = sites.comm
    for p = 1:comm.patches
        newborns = sites.nbr_newborns[p] < 1 ? (rand() <= sites.nbr_newborns[p] ? 1 : 0) :
            Int(round(sites.nbr_newborns[p]))
        nbr_settled = min(newborns, sites.nbr_empty[p])
        # pos_offspring = shuffle(sites.empty_pos[1:sites.nbr_empty[p], p])[1:nbr_settled]
        pos_offspring = sample(sites.empty_pos[1:sites.nbr_empty[p], p], nbr_settled; replace = false)
        # mothers_long = wsample(1:comm.tot_sites, (@view sites.prob_mothers[:, p]), nbr_settled)
        # rands = rand(nbr_settled)
        mothers_long = [binary_search((@view sites.prob_mothers[:, p]), r) for r in rand(nbr_settled).*sites.nbr_newborns[p]]
        ps_parents = @. div(mothers_long - 1, comm.microsites) + 1
        mothers = @. mothers_long - comm.microsites * (ps_parents - 1)
        for i in 1:nbr_settled
            offspring = pos_offspring[i]
            p_parents = ps_parents[i]
            mother = mothers[i]
            s = sites.species_ID[mother]
            if comm.rep_type == :SEXUAL
                father = wsample(sites.prob_fathers[:, s, p_parents])
                inherit(sites, offspring, mother, father, p, p_parents)
            else
                inherit(sites, offspring, mother, p, p_parents)
            end
            mutate(sites, i, p)
            new_body(sites, offspring, p, s)
        end
    end
end

function inherit(sites::Sites, i, mother, p, p_mother)
    sites.genotype[:, i, p] .= sites.genotype[:, mother, p_mother]
end

function inherit(sites::Sites, i, mother, father, p, p_mother)
    evol = sites.evol
    genes = rand(evol.trait_loci) .< 0.5
    sites.genotype[evol.all_mother[genes], i, p] .=
        sites.genotype[evol.all_mother[genes], mother, p_mother]
    sites.genotype[evol.all_mother[.!genes], i, p] .=
        sites.genotype[evol.all_father[.!genes], mother, p_mother]
    genes = rand(evol.trait_loci) .< 0.5
    sites.genotype[evol.all_father[genes], i, p] .=
        sites.genotype[evol.all_mother[genes], father, p_mother]
    sites.genotype[evol.all_father[.!genes], i, p] .=
        sites.genotype[evol.all_father[.!genes], father, p_mother]
end

function mutate(sites::Sites, i, p)
    evol = sites.evol
    muts = rand(evol.tot_genes) .<= evol.mu
    # muts = findall(rand(evol.tot_genes) .<= evol.mu)
    muts_sum = sum(muts)
    if muts_sum > 0
        sites.genotype[muts, i, p] .+= rand([-1,1] , muts_sum)
    end
end


# function binary_search(lst::Vector{T}, val::T) where T
#     low = 1
#     high = length(lst)
#     # mid = 0
#     while low ≤ high
#         mid = (low + high) ÷ 2
#         if lst[mid] >= val
#             high = mid - 1
#         else lst[mid] < val
#             low = mid + 1
#         end
#     end
#     return low
# end

function binary_search(lst, val)
    low = 1
    high = length(lst)
    # mid = 0
    while low ≤ high
        mid = (low + high) ÷ 2
        if lst[mid] >= val
            high = mid - 1
        else lst[mid] < val
            low = mid + 1
        end
    end
    return low
end

# function binary_search(lst::Array{Float64, 1}, val::Float64)
#     low = 1
#     high = length(lst)
#     # mid = 0
#     while low ≤ high
#         mid = (low + high) ÷ 2
#         if lst[mid] >= val
#             high = mid - 1
#         else lst[mid] < val
#             low = mid + 1
#         end
#     end
#     return low
# end


function evolving_foodweb(; time_steps = 1, print_steps = 1)

    run = Run(; time_steps = time_steps, print_steps = print_steps, dt = 0.1)
    comm = Comm(; grid = (X = 5, Y = 2), env_step = 0.001, microsites = 20000, m = 0.1, species = 30, rep_type = :ASEXUAL, dt = run.dt, scale_uptake = 1.)
    evol = Evol()
    sites = Sites(comm, evol)

    for t in Int.(1:(run.time_steps/run.dt))
        if (t * run.dt) % run.print_steps == 0
            println("time = ", t * run.dt)
        end
        change_environment(sites)
        uptake(sites)
        update_prey(sites)
        contribution_adults(sites)
        reproduction(sites)

        # comm.env_step
        # comm.bodymass_bc
        # comm.uptake_pars
        # comm.uptake_bc
        # sites.resource
        # sites.mass_abundance

    end

    return sites
end


@profiler sites = evolving_foodweb(; time_steps = 1)
@profiler sites = evolving_foodweb(; time_steps = 100, print_steps = 100);
@progress sites = evolving_foodweb()
@time sites = evolving_foodweb()
@time sites = evolving_foodweb(; time_steps = 100, print_steps = 100);

comm = sites.comm

sites.resource
uptake(sites)
sites.resource
sites.mass_abundance
update_prey(sites)
sites.mass_abundance
contribution_adults(sites)
sites.mass_abundance
reproduction(sites)
sites.mass_abundance


comm.env_step
comm.bodymass_bc
comm.uptake_pars
comm.uptake_bc
