
include("../src/evolving_foodweb_direct_method4.jl")

init_values = Init_values(;
    ## ecological input
    # patches
    grid = (X = 5, Y = 5),
    torus = (X = :NO, Y = :YES),
    env_range = (X = (-0.4, 0.4), Y = (-0.05, 0.05)),
    CC_vec = [0.0, 0.1, 0.2, 0.4],
    time_CC_vec = [500, 1000, 2000, 5000],
    # env_step_CC = 0.0002,
    # time_CC = 100,
    env_step_local = 0.001,
    dt_env = 1.,
    # dispersal
    m_vec = [0.01, 0.03, 0.1],
    rho = 2.0,
    m_tl = :EQUAL,
    # resource
    resource = 40.,
    in_rate = 40.,
    out_rate = 0.1,
    # species
    N = 600,
    spec_dist = :PATCH,
    rep_type = :SEXUAL,
    # trophic levels
    trophic_levels = 1,
    bm_offset = 1.,
    bm_power = 2.,
    # mortality
    d = 0.2,
    d_power = -0.5,
    # feeding
    uptake_pars = [0.0001, 0.4, 1.],
    i_power = 0.75,
    resource_conversion = 1.,
    resource_assimilation = 10.,
    assimilation_eff = 0.75,
    scale_uptake = 2,
    scale_assim = 1.,

    ## evolutionary input
    omega_e = 0.25,
    trait_loci = 10,
    mu = 1e-4,
    # mu = 0.,
    sigma_z = 0.05,

    ## run input
    runs = 2,
    pre_post_change = 5000,
    print_steps = 1000,
    log_steps = 100,
    output_file = "results/output_X5_Y5_tl1_loci10_mu1e-4_oe25e-2_d2e-1_dp5e-1_Yrange1e-1_stepLocal1e-3_specPatch_sex_CCrange_TCCrange.csv"
);

# @exportAll()
# end
# using .Evolving_foodweb

@time evolving_foodweb_dm(init_values);

# @profiler evolving_foodweb_dm(init);
