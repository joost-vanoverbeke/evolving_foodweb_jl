

# rand(1:100, 10)

# a = [x * y for x in 1:10, y in 1:10]
# b = collect(1:10)
# a
# a*b
# b'*a
# b*b'
# b'*b

resource = 150.
in_rate = 150.
out_rate = 0.1

species = 10
match_sd = 1.
trophic_levels = 2

bm_offset = 1.
bm_power = 2.
# mortality
d = 0.2
d_power = -0.5 #-0.25
# feeding
uptake_pars = [0.0001, 0.4, 1.]
i_power = 0.75 #0.75
resource_conversion = 1
resource_assimilation = 10. 
assimilation_eff = 0.75 #0.7
scale_uptake = 2
scale_assim = 1.

in_rate = in_rate*scale_uptake
uptake_pars[1] /= scale_uptake

tl_species = [(s - 1) % trophic_levels + 1 for s = 1:species]
species_per_tl = species/trophic_levels
match_resource = [tl_species[s] == 1 ? 1 : 0 for s in 1:species]
# match_shuffle = round.(rand(convert(Int, species/trophic_levels)), digits=2)
match_shuffle = [0.11,0.3,0.42,0.57,0.87]
species_match = [match_shuffle[(s - 1) ÷ trophic_levels + 1] for s = 1:species]
match_matr = [tl_species[s1] == tl_species[s2] + 1 ? exp(-(species_match[s1] - species_match[s2])^2 / (2*match_sd^2)) : 0 for s1 in 1:species, s2 in 1:species]
match_rel = match_matr./[x == 0 ? 1 : x for x in sum(match_matr; dims=2)]

# s1 = 2
# s2 = 9

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
        # conversion_tl[tl] = resource_assimilation/bodymass_tl[tl]
        conversion_tl[tl] = resource_assimilation/(bodymass_tl[tl]^(1/scale_assim))
    else
        # conversion_tl[tl] = assimilation_eff*(1-scale_assim*bodymass_tl[tl]^(-i_power)) * bodymass_tl[tl-1]/bodymass_tl[tl]
        conversion_tl[tl] = assimilation_eff * bodymass_tl[tl-1]/(bodymass_tl[tl]^(1/scale_assim))
    end
end

N_s = zeros(Int64, species)
N_s[[1,3,5,7,9]] .= [1000,2500,1500,1300,2100]
N_s[[2,4,6,8,10]] .= [150,180,200,170,210]
# N_s[[1,4,7]] .= [2000,3000,1000]
# N_s[[2,5,8]] .= [100,200,300]
# N_s[[3,6,9]] .= [30,10,20]
N_tl = zeros(Int64, trophic_levels)
N_tl[1] =  sum(N_s[[1,3,5,7,9]])
N_tl[2] =  sum(N_s[[2,4,6,8,10]])
# N_tl[1] =  sum(N_s[[1,4,7]])
# N_tl[2] =  sum(N_s[[2,5,8]])
# N_tl[3] =  sum(N_s[[3,6,9]])
total_N = sum(N_tl)

gain_s = zeros(species)
loss_s = zeros(species)
gain_tl = zeros(trophic_levels)
loss_tl = zeros(trophic_levels)

uptake = zeros(species)
total_uptake_tl = zeros(trophic_levels)

match_N = match_matr .* N_s'
N_prey = sum(match_N; dims=2)
N_ratio = N_s./(N_tl[tl_species[1:species]]/(species/trophic_levels))
for s in 1:species
    tl = tl_species[s]
    if tl == 1
        upt = uptake_tl[tl, 1] * resource^uptake_pars[3]
    else
        upt = uptake_tl[tl, 1] * N_prey[s]^uptake_pars[3]
    end
    uptake[s] = upt/(1. + upt*uptake_tl[tl, 2])
    gain_s[s] = uptake[s] * conversion_tl[tl]
end
total_uptake = uptake .* N_s
match_loss = total_uptake .* match_rel .* N_ratio'
prey_loss = vec(sum(match_loss; dims = 1))
loss_s .= prey_loss ./ N_s
resource_uptake = sum(total_uptake .* match_resource)
resource_gain = in_rate
resource_loss = out_rate*resource + resource_uptake


s = 1
# total_uptake = zeros(Float64, species)
dir = :UP
N_s[s] += 1
N_tl[tl_species[s]] += 1

if s == 0
    resource_uptake = 0.
    tl = 1
    s_comp = (1:species)[tl_species .== tl]
    upt = uptake_tl[tl, 1] * resource^uptake_pars[3]
    upt = upt/(1. + upt*uptake_tl[tl, 2])
    for s2 in s_comp 
        uptake[s2] = upt
        gain_s[s2] = uptake[s2] * conversion_tl[tl]
        total_uptake[s2] = uptake[s2] .* N_s[s2]
        resource_uptake += total_uptake[s2]
    end
    resource_loss = out_rate*resource + resource_uptake
elseif s > 0
    tl_comp = tl_species[s]
    tl_prey = tl_comp - 1
    tl_pred = tl_comp + 1
    s_prey = (1:species)[tl_species .== tl_prey]
    s_comp = (1:species)[tl_species .== tl_comp]
    s_pred = (1:species)[tl_species .== tl_pred]
    N_ratio[s_comp] .= N_s[s_comp]./(N_tl[tl_comp]/species_per_tl)
    tl_ratio = species_per_tl/N_tl[tl_comp]
    if tl_comp == 1
        if dir == :UP
            resource_loss += uptake[s]
            N_prey += match_matr[:, s]
        else
            resource_loss -= uptake[s]
            N_prey -= match_matr[:, s]
        end
        for s2 in s_pred 
            upt = uptake_tl[tl_pred, 1] * N_prey[s2]^uptake_pars[3]
            uptake[s2] = upt/(1. + upt*uptake_tl[tl_pred, 2])
            gain_s[s2] = uptake[s2] * conversion_tl[tl_pred]
            total_uptake[s2] = uptake[s2] * N_s[s2] * tl_ratio
        end
        match_loss = total_uptake[s_pred] .* match_rel[s_pred, s_comp]
        prey_loss = vec(sum(match_loss; dims = 1))
        loss_s[s_comp] .= prey_loss .* N_ratio[s_comp] ./ N_s[s_comp]
        loss_s[s_comp] .= prey_loss .* tl_ratio


        
    elseif tl_comp > 1 && tl_comp < trophic_levels
        if dir == :UP
            N_prey += match_matr[:, s]
        else
            N_prey -= match_matr[:, s]
        end
        for s2 in s_pred 
            upt = uptake_tl[tl_pred, 1] * N_prey[s2]^uptake_pars[3]
            uptake[s2] = upt/(1. + upt*uptake_tl[tl_pred, 2])
            gain_s[s2] = uptake[s2] * conversion_tl[tl_pred]
            total_uptake[s2] = uptake[s2] .* N_s[s2]
        end
        match_loss_old = total_uptake[s] .* match_rel[s, s_prey] .* N_ratio[s_prey]
        @. loss_s[s_prey] -= match_loss_old / N_s[s_prey]
        total_uptake[s] = uptake[s] .* N_s[s]
        match_loss = total_uptake[s] .* match_rel[s, s_prey] .* N_ratio[s_prey]
        @. loss_s[s_prey] += match_loss / N_s[s_prey]
        match_loss = total_uptake[s_pred] .* match_rel[s_pred, s_comp]
        prey_loss = vec(sum(match_loss; dims = 1)) .* N_ratio[s_comp]
        @. loss_s[s_comp] = prey_loss / N_s[s_comp]
    else
        loss_s[s] = 0.
        match_loss_old = total_uptake[s] .* match_rel[s, s_prey] .* N_ratio[s_prey]
        @. loss_s[s_prey] -= match_loss_old / N_s[s_prey]
        total_uptake[s] = uptake[s] .* N_s[s]
        match_loss = total_uptake[s] .* match_rel[s, s_prey] .* N_ratio[s_prey]
        @. loss_s[s_prey] += match_loss / N_s[s_prey]
    end
end
