using ArgParse;
using Distributions;
using LinearAlgebra;
using Logging;
using DelimitedFiles;

include(joinpath("inc", "eap_chain.jl"));
include(joinpath("inc", "average.jl"));
include(joinpath("inc", "acceptance.jl"));

s = ArgParseSettings();
@add_arg_table! s begin
  "--E0", "-z"
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
  "--energy-type", "-e"
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
  "--theta-step", "-q"
    help = "maximum θ step length"
    arg_type = Float64;
    default = 3*π / 16;
  "--chain-frac-step", "-f"
    help = "fraction of monomers to step (end-to-end ensemble)"
    arg_type = Float64
    default = 0.10
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
  "--acc", "-x"
    help = "acceptance function (metropolis|kawasaki)"
    arg_type = String
    default = "metropolis"
  "--umbrella-sampling", "-B"
    help = "use umbrella sampling (w/ electrostatic weight function)"
    action = :store_true
  "--update-freq", "-u"
    help = "update frequency (seconds)"
    arg_type = Float64;
    default = 15.0;
  "--verbose", "-v"
    help = "verbosity level: 0-nothing, 1-errors, 2-warnings, 3-info"
    arg_type = Int
    default = 3
  "--callbacks", "-C"
    help = "file with callback functions"
    arg_type = String
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

callbacks = Any[];
try
  callbacks = if pargs["callbacks"] !== nothing
                evalfile(pargs["callbacks"])
              end
catch e
  @error "error occured loading $(pargs["callbacks"])";
  rethrow();
end

function mcmc(nsteps::Int, pargs, callbacks)
  ϕstep, θstep = pargs["phi-step"], pargs["theta-step"];
  dϕ_dist = Uniform(-ϕstep, ϕstep);
  dθ_dist = Uniform(-θstep, θstep);
  chain = EAPChain(pargs);
  chain.U = U(chain);
  acceptor = if pargs["acc"] == "metropolis"
    Metropolis(
               chain, 
               (pargs["umbrella-sampling"]) ? weight_anti_dipole : chain -> 1.0
              );
  else
    error("'$(pargs["acc"])' acceptance criteria has not yet been implemented.");
  end
  Avg, avgcons = if pargs["umbrella-sampling"]
    UmbrellaAverager, (acc, chain) -> UmbrellaAverager(acc, weight_anti_dipole, 
                                                       chain);
  else
    StandardAverager, StandardAverager
  end
  scalar_averagers = (Avg[
                       avgcons(chain -> dot(end_to_end(chain), 
                                            end_to_end(chain)), chain),
                       avgcons(chain -> dot(chain_μ(chain), chain_μ(chain)), 
                                        chain),
                       avgcons(chain -> chain.U, chain),
                       avgcons(chain -> chain.U*chain.U, chain),
                      ]);
  vector_averagers = (Avg{Vector{Float64}}[
                       avgcons(end_to_end, chain),
                       avgcons(chain -> map(x -> x*x, end_to_end(chain)), chain),
                       avgcons(chain_μ, chain),
                       avgcons(chain -> map(x -> x*x, chain_μ(chain)), chain)
                      ]);
  outfile = open("$(pargs["prefix"])_trajectory.csv", "w");

  start = time();
  last_update = start;
  num_accepted = 0;
  for init=1:pargs["num-inits"]

    for step=1:nsteps
      dϕ = rand(dϕ_dist);
      dθ = rand(dθ_dist);
      idx = rand(1:n(chain));
      trial_chain = EAPChain(chain);
      move!(trial_chain, idx, dϕ, dθ);
      if acceptor(trial_chain, rand())
        chain = trial_chain;
        num_accepted += 1;
      end

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
        acc_ratio = num_accepted / (step + (init - 1)*nsteps);
        if (acc_ratio > pargs["step-adjust-ub"] &&
            ϕstep != π && θstep != π/2)
          @info "acceptance ratio is high; increasing step size";
          ϕstep = min(π, ϕstep*pargs["step-adjust-scale"]);
          θstep = min(π/2, θstep*pargs["step-adjust-scale"]);
        elseif acc_ratio < pargs["step-adjust-lb"]
          @info "acceptance ratio is low; decreasing step size";
          ϕstep /= pargs["step-adjust-scale"];
          θstep /= pargs["step-adjust-scale"];
        end
        dϕ_dist = Uniform(-ϕstep, ϕstep);
        dθ_dist = Uniform(-θstep, θstep);
      end

      #for callback in callbacks
      #  callback(chain, step, false, false, pargs);
      #end
      if step % pargs["stepout"] == 0
        writedlm(outfile, 
                 hcat(step, transpose(chain.r), 
                      transpose(chain_μ(chain)), chain.U), 
                 ',');
      end
      foreach(a -> record!(a, chain), scalar_averagers);
      foreach(a -> record!(a, chain), vector_averagers);
    
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
  
  ar = num_accepted / (pargs["num-inits"]*pargs["num-steps"]);
  @info "total time elapsed: $(time() - start)";
  @info "acceptance rate: $ar";

  #for callback in callbacks
  #  callback(chain, pargs["num-steps"], false, true, pargs);
  #end
  close(outfile);

  return (scalar_averagers, vector_averagers, ar);

end

(sas, vas, ar) = mcmc(pargs["num-steps"], pargs, callbacks);
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
