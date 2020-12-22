include(joinpath("inc", "algorithm_help_msg.jl"));

using ArgParse;
using Plots;

s = ArgParseSettings();
@add_arg_table s begin
  "--E0", "-z"
    help = "magnitude of electric field"
    arg_type = Float64
    default = 0.0
  "--dk", "-D"
    help = "kappa2 - kappa1"
    arg_type = Float64
    default = 1.0
  "--kT", "-k"
    help = "dimensionless temperature"
    arg_type = Float64
    default = 1.0
  "--gz", "-g"
    help = "dimensionless stretch in the z-direction (direction of E-field)"
    arg_type = Float64
    default = 0.0
  "--gx", "-f"
    help = "dimensionless stretch in the x-direction"
    arg_type = Float64
    default = 0.0
  "--num-monomers", "-N"
    help = "number of monomers"
    arg_type = Int
    default = 100
  "--algorithm", "-A"
    help = "algorithm (ISRES | PRAXIS | SBPLX | SLSQP | LBFGS | AUGLAG | Newton | Anneal). Run with --print-algo-help for more detailed information"
    default = "try Newton then PRAXIS,SBPLX,Newton"
  "--print-algo-help"
    help = "Get detailed information regarding algorithms themselves, and how to run solvers in stages"
    action = :store_true
  "--print-every-percent-compl"
    help = "Print diagnostic info every percentage complete of maximum iterations (should be between 0 and 1)"
    arg_type = Float64
    default = 0.1
  "--verbose", "-v"
    help = "print diagnostic information"
    action = :store_true
  "--reltol-int", "-r"
    help = "relative tolerance for numerical integration"
    arg_type = Float64
    default = 1e-6
  "--maxevals-int", "-e"
    help = "maximum evaluations for numerical integration"
    arg_type = Int
    default = 2500
  "--reltol-opt", "-R"
    help = "relative tolerance for optimization"
    arg_type = Float64
    default = 0.0
  "--abstol-opt", "-B"
    help = "absolute tolerance for optimization"
    arg_type = Float64
    default = 1e-8
  "--maxevals-opt", "-E"
    help = "maximum evaluations for numerical integration"
    arg_type = String
    default = "try 1000 then 10000,5000,50"
  "--init-guess", "-X"
    help = "initial guess"
    arg_type = String
  "--anneal-factor", "-m"
    help = "annealing factor, gamma i.e. T_new = T * (1 - gamma)"
    arg_type = Float64
    default = 0.005
  "--anneal-init-T", "-t"
    help = "initial temperature"
    arg_type = Float64
    default = 25.0
  "--anneal-final-T", "-T"
    help = "final temperature"
    arg_type = Float64
    default = 0.0
  "--anneal-relax-steps", "-S"
    help = "number of steps before annealing and step size adaptation"
    arg_type = Int
    default = 100
  "--anneal-init-delta", "-s"
    help = "initial step size for annealing"
    arg_type = Float64
    default = 0.5
  "--anneal-adapt-factor", "-P"
    help = "step size adapation factor (should be greater than 1.0)"
    arg_type = Float64
    default = 1.2
  "--anneal-acc-thrshld", "-a"
    help = "step size adapation acceptance threshold (should be in [0, 1))"
    arg_type = Float64
    default = 0.7
  "--plot-dist", "-p"
    help = "plot monomer distribution"
    action = :store_true
end

pa = parse_args(s);
if pa["print-algo-help"]
  println(_ALGO_HELP_MSG_);
  exit();
end

const dk = pa["dk"];
const E0 = pa["E0"];
const kT = pa["kT"];
const ω0 = dk * E0*E0 / (kT); 
const _u0_ = (dk > 0) ? -dk*(E0^2) : 0.0; 
const γz = pa["gz"];
const γx = pa["gx"];
const N = pa["num-monomers"];
const reltol_opt = pa["reltol-opt"];
const abstol_opt = pa["abstol-opt"];

const verbose = pa["verbose"];
const print_every_percent_compl = pa["print-every-percent-compl"];
const RELTOL_INT = pa["reltol-int"];
const MAXEVALS_INT = pa["maxevals-int"];
const ANNEAL_FACTOR = pa["anneal-factor"];
const ANNEAL_INIT_T = pa["anneal-init-T"];
const ANNEAL_FINAL_T = pa["anneal-final-T"];
const ANNEAL_RELAX_STEPS = pa["anneal-relax-steps"];
const ANNEAL_INIT_DELTA = pa["anneal-init-delta"];
const ANNEAL_ADAPT_FACTOR = pa["anneal-adapt-factor"];
const ANNEAL_ACC_THRSHLD = pa["anneal-acc-thrshld"];

include(joinpath("inc", "solvers.jl"));
include(joinpath("inc", "free_energy.jl"));
include(joinpath("inc", "dipole.jl"));
include(joinpath("inc", "output.jl"));

_RHS_ = Float64[N; γz*N; γx*N];
ω(ϕ::Real, θ::Real) = (ω0 * (cos(θ))^2);

main_algos = parse_algostr(pa["algorithm"]);
main_max_iters = parse_max_iters(pa["maxevals-opt"]);
 DEFAULT_INIT_GUESS = (
  if pa["init-guess"] != nothing
    eval(parse(pa["init-guess"]));
  #=elseif ω == 0
    initial_guess(γx, γz);
  elseif abs(ω) < 1 
    ((1-abs(ω))*initial_guess(γx, γz) + 
     (abs(ω))*initial_guess(ω, (γx, γz)));
  else
    initial_guess(ω, (γx, γz));=#
  else
    initial_guess(γx, γz);
  end
);
(xs, rs, converged, iterations) = run_solver(main_algos, main_max_iters;
                                             x0s = DEFAULT_INIT_GUESS,
                                             reltol = reltol_opt,
                                             abstol = abstol_opt);
if verbose
  println((converged) ? "Solution converged." : "Solution did NOT converge.");
  println("Iterations: ", iterations);
end
residuals_output(rs);
standard_solution_output(xs);
if pa["plot-dist"]
  ρ_polar(θ) = (pquadrature(ϕ -> ρ(ϕ, θ, xs[1], xs[2], xs[3]), 0.0, 2*π))[1];
  p = plot(range(0, 2*π; length=1000), ρ_polar; proj=:polar);
  display(p);
  println("Press <Enter> to continue...");
  readline();
end
