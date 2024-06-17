using ArgParse;
using Distributions;
using LinearAlgebra;
using Logging;
using DelimitedFiles;
#using ProfileView;
using DecFP;
using Quadmath;

include(joinpath("inc", "eap_chain.jl"));
include(joinpath("inc", "average.jl"));
include(joinpath("inc", "acceptance.jl"));

s = ArgParseSettings();
@add_arg_table! s begin
  "--E0", "-e"
    help = "magnitude of electric field"
    arg_type = Float64
    default = 0.0
  "--chain-type", "-T"
    help = "chain type (dielectric|polar)"
    arg_type = String
    default = "dielectric"
  "--K1", "-J"
    help = "dipole susceptibility along the monomer axis (dielectric chain)"
    arg_type = Float64
    default = 1.0
  "--K2", "-K"
    help = "dipole susceptibility orthogonal to the monomer axis (dielectric chain)"
    arg_type = Float64
    default = 0.0
  "--mu", "-m"
    arg_type = Float64
    default = 1e-2
    help = "dipole magnitude (electret chain)"
  "--energy-type", "-u"
    help = "energy type (noninteracting|interacting|Ising)"
    arg_type = String
    default = "noninteracting"
  "--kT", "-k"
    help = "dimensionless temperature"
    arg_type = Float64
    default = 1.0
  "--Fz", "-F"
    help = "force in the z-direction (direction of E-field; force ensemble)"
    arg_type = Float64
    default = 0.0
  "--Fx", "-G"
    help = "force in the x-direction (force ensemble)"
    arg_type = Float64
    default = 0.0
  "--mlen", "-b"
    help = "monomer length"
    arg_type = Float64
    default = 1.0
  "--num-monomers", "-n"
    help = "number of monomers"
    arg_type = Int
    default = 100
  "--num-steps", "-N"
    help = "number of steps"
    arg_type = Int
    default = convert(Int, 1e6)
  "--phi-step", "-p"
    help = "maximum ϕ step length"
    arg_type = Float64;
    default = 3*π / 8;
#  "--do-clustering"
#    help = "trial moves with clustering"
#    action = :store_true
  "--cluster-prob"
    help = "probability of flipping a cluster"
    arg_type = Float64;
    default = 0.5;
  "--step-adjust-lb", "-L"
    help = "adjust step sizes if acc. ratio below this threshold"
    arg_type = Float64
    default = 0.15
  "--step-adjust-ub", "-U"
    help = "adjust step sizes if acc. ratio above this threshold"
    arg_type = Float64
    default = 0.40
  "--step-adjust-scale", "-A"
    help = "scale factor for adjusting step sizes (> 1.0)"
    arg_type = Float64
    default = 1.1
  "--steps-per-adjust", "-S"
    help = "steps between step size adjustments"
    arg_type = Int
    default = 2500
  "--umbrella-sampling", "-B"
    help = "use umbrella sampling (w/ electrostatic weight function)"
    action = :store_true
  "--update-freq"
    help = "update frequency (seconds)"
    arg_type = Float64;
    default = 15.0;
  "--verbose", "-v"
    help = "verbosity level: 0-nothing, 1-errors, 2-warnings, 3-info"
    arg_type = Int
    default = 3
  "--prefix", "-P"
    help = "prefix for output files"
    arg_type = String
    default = "eap-mcmc"
  "--postfix", "-Q"
    help = "postfix for output files"
    arg_type = String
    default = ""
  "--stepout", "-s"
    help = "steps between storing microstates"
    arg_type = Int
    default = 500
  "--numeric-type"
    help = "numerical data type for averaging (float64|float128|dec128|big)"
    arg_type = String
    default = "float64"
  "--profile", "-Z"
    help = "profile the program"
    action = :store_true
  "--burn-in"
    help = "steps for burn-in; i.e. steps before averaging"
    arg_type = Int
    default = 50000
end

pargs = parse_args(s);

if pargs["verbose"] == 3
  global_logger(ConsoleLogger(stderr, Logging.Info));
elseif pargs["verbose"] == 2
  global_logger(ConsoleLogger(stderr, Logging.Warn));
elseif pargs["verbose"] == 1
  global_logger(ConsoleLogger(stderr, Logging.Error));
else
  global_logger(Logging.NullLogger());
end

function mcmc(nsteps::Int, pargs)
  chain = EAPChain(pargs);
  return mcmc(nsteps, pargs, chain);
end

function mcmc(nsteps::Int, pargs, chain::EAPChain)  
  ϕstep = pargs["phi-step"];
  dϕ_dist = Uniform(-ϕstep, ϕstep);
  chain = EAPChain(pargs);
  chain.U = U(chain);
  wf = AntiDipoleWeightFunction(chain);
  acceptor = Metropolis(
               chain, 
               (pargs["umbrella-sampling"]) ? wf : WeightlessFunction()
              );
  numeric_type = if pargs["numeric-type"] == "float64"
    Float64;
  elseif pargs["numeric-type"] == "float128"
    Float128;
  elseif pargs["numeric-type"] == "dec128"
    Dec128;
  elseif pargs["numeric-type"] == "big"
    BigFloat;
  else
    error("numeric-type '$(pargs["numeric-type"])' not understood");
    exit(-1);
  end

  # setup averagers
  Avg, avgcons, avgconss = if pargs["umbrella-sampling"]
    (UmbrellaAverager, 
      ((acc, chain) -> begin;
        UmbrellaAverager(
                         StandardAverager(
                                          numeric_type(0.0*acc(chain)),
                                          numeric_type(0.0),
                                          acc
                                         ),
                         wf
                        );
     end),
     ((acc, chain) -> begin;
        UmbrellaAverager(
                         StandardAverager(
                                          Vector{numeric_type}(0.0*acc(chain)),
                                          numeric_type(0.0),
                                          acc
                                         ),
                         wf
                        );
     end)
     );
  else
    (StandardAverager,
      ((acc, chain) -> begin;
             StandardAverager(
                              numeric_type,
                              acc, chain
                             )
     end),
     ((acc, chain) -> begin;
             StandardAverager(
                              numeric_type,
                              Vector{numeric_type},
                              acc, chain
                             );
     end)
     );
  end
  scalar_averagers = (Avg{numeric_type,numeric_type}[
                       (avgcons(chain -> dot(end_to_end(chain), 
                                             end_to_end(chain)), chain)),
                       (avgcons(chain -> dot(chain_μ(chain), chain_μ(chain)), 
                                chain)),
                       (avgcons(chain -> chain.U, chain)),
                       (avgcons(chain -> chain.U*chain.U, chain)),
                      ]);
  vector_averagers = (Avg{Vector{numeric_type},numeric_type}[
                       avgconss(end_to_end, chain),
                       avgconss(chain -> map(x -> x*x, end_to_end(chain)), chain),
                       avgconss(chain_μ, chain),
                       avgconss(chain -> map(x -> x*x, chain_μ(chain)), chain)
                      ]);

  outfile = open("$(pargs["prefix"])_trajectory.csv", "w");
  writedlm(outfile, ["step" "r1" "r3" "p1" "p3" "U"], ',');
  rollfile = open("$(pargs["prefix"])_rolling.csv", "w");
  writedlm(rollfile, ["step" "r1" "r3" "r1sq" "r3sq" "rsq" "p1" "p3" "p1sq" "p3sq" "psq" "U" "Usq"], ',');

  start = time();
  last_update = start;
  nacc = 0;
  nacc_total = 0;
  natt = 0;
      
  for step=1:nsteps
    idx = rand(1:n(chain));
    dϕ = rand(dϕ_dist);
    trial_chain = EAPChain(chain);
    successful = move!(trial_chain, idx, dϕ);
    α = cluster_flip!(trial_chain, idx; ϵflip = pargs["cluster-prob"]);
    if successful && acceptor(trial_chain, rand(); α = α)
      chain = trial_chain;
      nacc += 1;
      nacc_total += 1;
    end
    natt += 1;

    if time() - last_update > pargs["update-freq"]
      @info "elapsed: $(time() - start)";
      @info "step:    $step / $nsteps";
      last_update = time();
    end

    if (
        pargs["step-adjust-scale"] != 1.0 &&
        step % pargs["steps-per-adjust"] == 0
       ) # adjust step size?
      ar = nacc / natt;
      if (ar > pargs["step-adjust-ub"] &&
          ϕstep != π)
        @info "acceptance ratio is high; increasing step size";
        nacc = 0;
        natt = 0;
        ϕstep = min(π, ϕstep*pargs["step-adjust-scale"]);
      elseif ar < pargs["step-adjust-lb"]
        @info "acceptance ratio is low; decreasing step size";
        nacc = 0;
        natt = 0;
        ϕstep /= pargs["step-adjust-scale"];
      end
      dϕ_dist = Uniform(-ϕstep, ϕstep);
    end

    foreach(a -> record!(a, chain), scalar_averagers);
    foreach(a -> record!(a, chain), vector_averagers);
    if step % pargs["stepout"] == 0
      writedlm(outfile, 
               hcat(step, transpose(chain.r), 
                    transpose(chain_μ(chain)), chain.U), 
               ',');
      writedlm(rollfile, 
               hcat(
                    step,
                    transpose(get_avg(vector_averagers[1])), # r
                    transpose(get_avg(vector_averagers[2])), # rj2
                    get_avg(scalar_averagers[1]), # r2
                    transpose(get_avg(vector_averagers[3])), # p
                    transpose(get_avg(vector_averagers[4])), # pj2
                    get_avg(scalar_averagers[2]), # p2
                    get_avg(scalar_averagers[3]), # U
                    get_avg(scalar_averagers[4])  # U2
                   ),
               ',');

    end
  
  end # steps

  AR = nacc_total / pargs["num-steps"];
  @info "total time elapsed: $(time() - start)";
  @info "acceptance rate: $AR";

  close(outfile);
  close(rollfile);

  return (chain, scalar_averagers, vector_averagers, AR);
end

(chain, sas, vas, ar) = if pargs["profile"]
  @info "Profiling the mcmc code...";
  mcmc(5, pargs); # run first to compile code
  #@profview mcmc(pargs["num-steps"], pargs);
  error("Not currently implemented...");
  #=
  println("Press ENTER to continue...");
  readline(stdin);
  exit(1);
  =#
else
  burned_in_chain = mcmc(pargs["burn-in"], pargs)[1];
  mcmc(pargs["num-steps"], pargs, burned_in_chain);
end

println("<r>    =   $(get_avg(vas[1]))");
println("<r/nb> =   $(get_avg(vas[1]) / (pargs["mlen"]*pargs["num-monomers"]))");
println("<rj2>  =   $(get_avg(vas[2]))");
println("<r2>   =   $(get_avg(sas[1]))");
println("<p>    =   $(get_avg(vas[3]))");
println("<pj2>  =   $(get_avg(vas[4]))");
println("<p2>   =   $(get_avg(sas[2]))");
println("<U>    =   $(get_avg(sas[3]))");
println("<U2>   =   $(get_avg(sas[4]))");
println("AR     =   $ar");
