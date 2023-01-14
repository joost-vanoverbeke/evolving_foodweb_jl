
include("evolving_foodweb_direct_method4.jl")

init_values = Init_values(;
    ## ecological input
    # patches
    grid = (X = 1, Y = 3),
    torus = (X = :NO, Y = :YES),
    env_range = (X = (-0.001, 0.001), Y = (-0.05, 0.05)),
#    env_step_CC = 0.0,
#    time_CC = 0,
    env_step_CC = 0.0002,
    time_CC = 1000,
    env_step_local = 0.0005,
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
    rep_type = :SEXUAL,
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
    runs = 10,
    pre_post_change = 5000,
    print_steps = 1000,
    log_steps = 100,
    output_file = "results/output_X1_Y3_tl1_loci10_oe5e-1_d2e-1_dp5e-1_Yrange1e-1_stepLocal5e-4_sex_CC2e-4_1e3.csv"
);

# @exportAll()
# end
# using .Evolving_foodweb

@time evolving_foodweb_dm(init_values);

# @profiler evolving_foodweb_dm(init);
