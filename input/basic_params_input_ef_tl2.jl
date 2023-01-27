
# include("evolving_foodweb_direct_method4.jl")

init_values = Init_values(;
    ## ecological input
    # patches
    grid = (X = 1, Y = 1),
    torus = (X = :NO, Y = :NO),
    env_range = (X = (-0.001, 0.001), Y = (-0.05, 0.05)),
   env_step_CC = 0.0,
   time_CC = 0,
    # env_step_CC = 0.0004,
    # time_CC = 1000,
    env_step_local = 0.002,
    dt_env = 1.,
    # dispersal
    m = 0.0,
    rho = 2.0,
    m_tl = :EQUAL,
    # resource
    resource = 220.,
    in_rate = 220.,
    out_rate = 0.1,
    # species
    N = 3000,
    spec_dist = :PATCH,
    rep_type = :SEXUAL,
    # trophic levels
    trophic_levels = 2,
    bm_offset = 1.,
    bm_power = 2.,
    # mortality
    d = 0.2,
    d_power = -0.5, #-0.25
    # feeding
    uptake_pars = [0.0001, 0.4, 1.],
    i_power = 0.75, #0.75
    resource_conversion = 1.,
    resource_assimilation = 10., 
    assimilation_eff = 0.75, #0.7
    scale_uptake = 2,
    scale_assim = 1.,

    ## evolutionary input
    omega_e = 0.25,
    trait_loci = 10,
    mu = 1e-4,
    # mu = 0.,
    sigma_z = 0.05,

    ## run input
    runs = 1,
    pre_post_change = 1000,
    print_steps = 1000,
    log_steps = 100,
    output_file = "results/output_test2.csv"
);

# @exportAll()
# end
# using .Evolving_foodweb

@time evolving_foodweb_dm(init_values);

# @profiler evolving_foodweb_dm(init);
