
init_1p_3tl = Init_values(;
    ## ecological input
    # patches
    grid = (X = 1, Y = 1),
    torus = (X = :NO, Y = :NO),
    env_range = (X = (-0.01, 0.01), Y = (-0.01, 0.01)),
    env_step_CC = 0. / 1000.,
    env_step_local = 0.0,
    dt_env = 1.,
    # dispersal
    m = 0.0,
    rho = 2.0,
    m_tl = :EQUAL,
    # resource
    resource = 200.,
    in_rate = 200.,
    out_rate = 0.1,
    # species
    N = 10000,
    rep_type = :SEXUAL,
    # trophic levels
    trophic_levels = 3,
    bm_offset = 1.,
    bm_power = 1.,
    # mortality
    d = 0.1,
    d_power = -0.25,
    # feeding
    uptake_pars = [0.0001, 0.4, 1.],
    i_power = 0.75,
    resource_conversion = 1.,
    resource_assimilation = 10.,
    assimilation_eff = 0.7,
    scale_uptake = 1,
    scale_assim = 0.,

    ## evolutionary input
    omega_e = 0.5,
    trait_loci = 10,
    mu = 1e-4,
    sigma_z = 0.1,

    ## run input
    runs = 1,
    time_steps = 10000,
    pre_change = 5000,
    post_change = 5000,
    print_steps = 1000,
    log_steps = 500,
    output_file = "output_1p_3tl.csv"
);

# @exportAll()
# end
# using .Evolving_foodweb

@time evolving_foodweb_dm(init_1p_3tl);

# @profiler evolving_foodweb_dm(init);
