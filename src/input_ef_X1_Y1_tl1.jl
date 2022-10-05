
init_X11_Y1_1tl = Init_values(;
    ## ecological input
    # patches
    grid = (X = 1, Y = 1),
    torus = (X = :NO, Y = :NO),
    env_range = (X = (-0.5, 0.5), Y = (-0.05, 0.05)),
    # env_step_CC = 0.0,
    env_step_CC = 0.0001,
    time_CC = 1000,
    env_step_local = 0.005,
    dt_env = 1.,
    # dispersal
    m = 0.01,
    rho = 2.0,
    m_tl = :EQUAL,
    # resource
    resource = 100.,
    in_rate = 100.,
    out_rate = 0.1,
    # species
    N = 3000,
    rep_type = :ASEXUAL,
    # trophic levels
    trophic_levels = 1,
    bm_offset = 1.,
    bm_power = 1.,
    # mortality
    d = 0.2,
    d_power = -0.5,
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
    # mu = 0.,
    sigma_z = 0.1,

    ## run input
    runs = 1,
    pre_post_change = 2000,
    print_steps = 500,
    log_steps = 100,
    output_file = "results/output_X1_Y1_tl1_d2e-1_dp5e-1_stCC1e-4_tCC1e+3.csv"
);

# @exportAll()
# end
# using .Evolving_foodweb

@time evolving_foodweb_dm(init_X11_Y1_1tl);

# @profiler evolving_foodweb_dm(init);
