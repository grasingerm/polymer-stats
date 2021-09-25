using ArgParse;
using Distributions;
using LinearAlgebra;
using Logging;
using DelimitedFiles;
using ProfileView;
using DecFP;
using Quadmath;

include(joinpath("inc", "eap_chain.jl"));
include(joinpath("inc", "average.jl"));
include(joinpath("inc", "acceptance.jl"));

function flip_n!(ψs::Vector)
  ψs[1] += π;
  ψs[2] = π - ψs[2];
end

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
    help = "energy type (noninteracting|interacting)"
    arg_type = String
    default = "noninteracting"
  "--kT", "-k"
    help = "dimensionless temperature"
    arg_type = Float64
    default = 1.0
  "--ensemble-type", "-E"
    help = "ensemble type (force|end-to-end)"
    arg_type = String
    default = "force"
  "--Fz", "-F"
    help = "force in the z-direction (direction of E-field; force ensemble)"
    arg_type = Float64
    default = 0.0
  "--Fx", "-G"
    help = "force in the x-direction (force ensemble)"
    arg_type = Float64
    default = 0.0
  "--rz", "-z"
    help = "end-to-end vector in the z-direction (direction of E-field; etoe ensemble)"
    arg_type = Float64
    default = 0.0
  "--rx", "-x"
    help = "end-to-end vector in the x-direction (etoe ensemble)"
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
    default = convert(Int, 1e5)
  "--num-inits", "-M"
    help = "number of random initializations"
    arg_type = Int
    default = 1
  "--force-init", "-I"
    help = "force each random initialization (false to use metro.)"
    action = :store_true
  "--phi-step", "-p"
    help = "maximum ϕ step length"
    arg_type = Float64;
    default = 3*π / 8;
  "--do-flips"
    help = "trial moves with flipping monomers"
    action = :store_true
  "--theta-step", "-q"
    help = "maximum θ step length"
    arg_type = Float64;
    default = 3*π / 16;
  "--chain-frac-step", "-f"
    help = "fraction of monomers to step (end-to-end ensemble)"
    arg_type = Float64
    default = 0.15
  "--step-adjust-lb", "-L"
    help = "adjust step sizes if acc. ratio below this threshold"
    arg_type = Float64
    default = 0.15
  "--step-adjust-ub", "-U"
    help = "adjust step sizes if acc. ratio above this threshold"
    arg_type = Float64
    default = 0.55
  "--step-adjust-scale", "-A"
    help = "scale factor for adjusting step sizes (> 1.0)"
    arg_type = Float64
    default = 1.1
  "--steps-per-adjust", "-S"
    help = "steps between step size adjustments"
    arg_type = Int
    default = 2500
  "--acc", "-a"
    help = "acceptance function (metropolis|kawasaki)"
    arg_type = String
    default = "metropolis"
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

if pargs["ensemble-type"] == "end-to-end"
  @warn "'end-to-end' ensemble is an experimental option; it has not been validated.";
end

function mcmc(nsteps::Int, pargs)
  ϕstep, θstep = pargs["phi-step"], pargs["theta-step"];
  dϕ_dist = Uniform(-ϕstep, ϕstep);
  dθ_dist = Uniform(-θstep, θstep);
  chain = EAPChain(pargs);
  chain.U = U(chain);
  wf = AntiDipoleWeightFunction(chain);
  acceptor = if pargs["acc"] == "metropolis"
    Metropolis(
               chain, 
               (pargs["umbrella-sampling"]) ? wf : WeightlessFunction();
              );
  else
    error("'$(pargs["acc"])' acceptance criteria has not yet been implemented.");
  end
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

  force_ensemble_flag = pargs["ensemble-type"] == "force";
  r0 = [pargs["rx"]; 0.0; pargs["rz"]];

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
  writedlm(outfile, ["step" "r1" "r2" "r3" "p1" "p2" "p3" "U"], ',');
  rollfile = open("$(pargs["prefix"])_rolling.csv", "w");
  writedlm(rollfile, ["step" "r1" "r2" "r3" "r1sq" "r2sq" "r3sq" "rsq" "p1" "p2" "p3" "p1sq" "p2sq" "p3sq" "psq" "U" "Usq"], ',');

  start = time();
  last_update = start;
  nacc = 0;
  nacc_total = 0;
  natt = 0;
  for init=1:pargs["num-inits"]
      
    if !force_ensemble_flag
      for i=1:10 
        if move!(chain, 1, 0, 0, r0; frac_mv=1.0) # start with admissible chain
          break;
        end
      end
    end

    for step=1:nsteps
      idx = rand(1:n(chain));
      dϕ = rand(dϕ_dist);
      dθ = (((pargs["do-flips"] && rand(Bool)) ? π - 2*chain.θs[idx] : 0) 
            + rand(dθ_dist));
      trial_chain = EAPChain(chain);
      successful = if force_ensemble_flag
        move!(trial_chain, idx, dϕ, dθ);
      else
        move!(trial_chain, idx, dϕ, dθ, r0; frac_mv=pargs["chain-frac-step"]);
      end
      if successful && acceptor(trial_chain, rand())
        chain = trial_chain;
        nacc += 1;
        nacc_total += 1;
      end
      natt += 1;

      if time() - last_update > pargs["update-freq"]
        @info "elapsed: $(time() - start)";
        @info "init:    $init / $(pargs["num-inits"])";
        @info "step:    $step / $nsteps";
        last_update = time();
      end

      if (
          pargs["step-adjust-scale"] != 1.0 &&
          step % pargs["steps-per-adjust"] == 0
         ) # adjust step size?
        acc_ratio = nacc / natt;
        if (acc_ratio > pargs["step-adjust-ub"] &&
            ϕstep != π && θstep != π/2)
          @info "acceptance ratio is high; increasing step size";
          nacc = 0;
          natt = 0;
          ϕstep = min(π, ϕstep*pargs["step-adjust-scale"]);
          θstep = min(π/2, θstep*pargs["step-adjust-scale"]);
        elseif acc_ratio < pargs["step-adjust-lb"]
          @info "acceptance ratio is low; decreasing step size";
          nacc = 0;
          natt = 0;
          ϕstep /= pargs["step-adjust-scale"];
          θstep /= pargs["step-adjust-scale"];
        end
        dϕ_dist = Uniform(-ϕstep, ϕstep);
        dθ_dist = Uniform(-θstep, θstep);
      end

      #for callback in callbacks
      #  callback(chain, step, false, false, pargs);
      #end
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

    # reinitialize polymer chain
    new_chain = EAPChain(pargs);
    new_chain.U = U(new_chain);
    if (
        pargs["force-init"] || 
        metropolis_acc(chain.kT, new_chain.U - chain.U, 
                       prod(chain.sθs), prod(new_chain.sθs), rand())
       )
      chain = new_chain;
    end

  end
  
  ar = nacc_total / (pargs["num-inits"]*pargs["num-steps"]);
  @info "total time elapsed: $(time() - start)";
  @info "acceptance rate: $ar";

  #for callback in callbacks
  #  callback(chain, pargs["num-steps"], false, true, pargs);
  #end
  close(outfile);
  close(rollfile);

  return (scalar_averagers, vector_averagers, ar);
end

(sas, vas, ar) = if pargs["profile"]
  @info "Profiling the mcmc code...";
  mcmc(5, pargs); # run first to compile code
  @profview mcmc(pargs["num-steps"], pargs);
  println("Press ENTER to continue...");
  readline(stdin);
  exit(1);
else
  mcmc(pargs["num-steps"], pargs);
end
total_steps = pargs["num-steps"] * pargs["num-inits"];

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
