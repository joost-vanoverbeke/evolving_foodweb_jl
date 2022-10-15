
include("evolving_foodweb_direct_method4.jl")

init_X11_Y1_1tl = Init_values(;
    ## ecological input
    # patches
    grid = (X = 11, Y = 4),
    torus = (X = :NO, Y = :YES),
    env_range = (X = (-1, 1), Y = (-0.1, 0.1)),
#    env_step_CC = 0.0,
#    time_CC = 0,
    env_step_CC = 0.0004,
    time_CC = 1000,
    env_step_local = 0.001,
    dt_env = 1.,
    # dispersal
    m = 0.01,
    rho = 2.0,
    m_tl = :EQUAL,
    # resource
    resource = 150.,
    in_rate = 150.,
    out_rate = 0.1,
    # species
    N = 3000,
    rep_type = :SEXUAL,
    # trophic levels
    trophic_levels = 2,
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
    omega_e = 1.,
    trait_loci = 10,
    mu = 1e-4,
    # mu = 0.,
    sigma_z = 0.1,

    ## run input
    runs = 10,
    pre_post_change = 10000,
    print_steps = 1000,
    log_steps = 100,
    output_file = "results/output_X11_Y4_tl2_loci10_oe1_d2e-1_dp5e-1_Yrange2e-1_stepLocal1e-3_sex_CC4e-4_1e3.csv"
);

# @exportAll()
# end
# using .Evolving_foodweb

@time evolving_foodweb_dm(init_X11_Y1_1tl);

# @profiler evolving_foodweb_dm(init);