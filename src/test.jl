

# rand(1:100, 10)

# a = [x * y for x in 1:10, y in 1:10]
# b = collect(1:10)
# a
# a*b
# b'*a
# b*b'
# b'*b

resource = 220.
in_rate = 220.
out_rate = 0.1

species = 9
match_sd = 1000.
trophic_levels = 3

bm_offset = 1.
bm_power = 1.
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
match_resource = [tl_species[s] == 1 ? 1 : 0 for s in 1:species]
species_match = [(s - 1) ÷ trophic_levels + 1 for s = 1:species]
match_matr = [tl_species[s1] == tl_species[s2] + 1 ? exp(-(species_match[s1] - species_match[s2])^2 / (2*match_sd^2)) : 0 for s1 in 1:species, s2 in 1:species]
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
        # conversion_tl[tl] = resource_assimilation/bodymass_tl[tl]
        conversion_tl[tl] = resource_assimilation/(bodymass_tl[tl]^(1/scale_assim))
    else
        # conversion_tl[tl] = assimilation_eff*(1-scale_assim*bodymass_tl[tl]^(-i_power)) * bodymass_tl[tl-1]/bodymass_tl[tl]
        conversion_tl[tl] = assimilation_eff * bodymass_tl[tl-1]/(bodymass_tl[tl]^(1/scale_assim))
    end
end

N_s = zeros(Int64, species)
N_s[[1,4,7]] .= [2000,3000,1000]
N_s[[2,5,8]] .= [100,200,300]
N_s[[3,6,9]] .= [30,10,20]
N_tl = zeros(Int64, trophic_levels)
N_tl[1] =  sum(N_s[[1,4,7]])
N_tl[2] =  sum(N_s[[2,5,8]])
N_tl[3] =  sum(N_s[[3,6,9]])
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



for tl in 1:trophic_levels
    if tl == 1
        upt = uptake_tl[tl, 1]*resource^uptake_pars[3]
    else
        upt = uptake_tl[tl, 1]*N_tl[tl-1]^uptake_pars[3]
    end
    total_uptake[tl] = upt/(1. + upt*uptake_tl[tl, 2])*N_tl[tl]
end
resource_gain_tl = in_rate
resource_loss_tl = out_rate*resource + total_uptake[1]
for tl in 1:trophic_levels
    gain_tl[tl] = total_uptake[tl]/N_tl[tl] * conversion_tl[tl]
    if tl < trophic_levels
        loss_tl[tl] = total_uptake[tl+1]/N_tl[tl]
    end
end



