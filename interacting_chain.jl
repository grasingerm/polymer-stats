using ArgParse;
using Distributions;
using LinearAlgebra;
using Logging;
using DelimitedFiles;

include("eap_chain.jl");

s = ArgParseSettings();
@add_arg_table! s begin
  "--mlen", "-b"
    help = "monomer length"
    arg_type = Float64
    default = 1.0
  "--E0", "-z"
    help = "magnitude of electric field"
    arg_type = Float64
    default = 0.0
  "--K1", "-K"
    help = "dipole susceptibility along the monomer axis"
    arg_type = Float64
    default = 1.0
  "--K2", "-L"
    help = "dipole susceptibility orthogonal to the monomer axis"
    arg_type = Float64
    default = 0.0
  "--kT", "-k"
    help = "dimensionless temperature"
    arg_type = Float64
    default = 1.0
  "--Fz", "-F"
    help = "force in the z-direction (direction of E-field)"
    arg_type = Float64
    default = 0.0
  "--Fx", "-G"
    help = "force in the x-direction"
    arg_type = Float64
    default = 0.0
  "--num-monomers", "-n"
    help = "number of monomers"
    arg_type = Int
    default = 100
  "--num-steps", "-N"
    help = "number of steps"
    arg_type = Int
    default = convert(Int, 1e5)
  "--phi-step", "-p"
    help = "maximum ϕ step length"
    arg_type = Float64;
    default = 3*π / 8;
  "--theta-step", "-q"
    help = "maximum θ step length"
    arg_type = Float64;
    default = 3*π / 16;
  "--update-freq", "-U"
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
  "--stepout", "-S"
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

metropolis_acc(kT::Real, dU::Real, ϵ::Real) = ( (dU <= 0) ? true : 
                                               ( ϵ <= exp(-dU / kT) ) );

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
  dϕ_dist = Uniform(-pargs["phi-step"], pargs["phi-step"]);
  dθ_dist = Uniform(-pargs["theta-step"], pargs["theta-step"]);
  chain = EAPChain(pargs);
  chain.U = U(chain);
  end_to_end_sum = chain.r[:];
  r2_sum = dot(chain.r, chain.r);
  chain_μ_sum = chain_μ(chain);
  Usum = chain.U;
  outfile = open("$(pargs["prefix"])_trajectory.csv", "w");
  #for callback in callbacks
  #  callback(chain, 0, true, false, pargs);
  #end

  start = time();
  last_update = start;
  num_accepted = 0;
  for step=1:nsteps
    dϕ = rand(dϕ_dist);
    dθ = rand(dϕ_dist);
    idx = rand(1:n(chain));
    ds = move!(chain, idx, dϕ, dθ);
    Ucurr = U(chain);
    if metropolis_acc(chain.kT, Ucurr - chain.U, rand())
      chain.U = Ucurr;
      num_accepted += 1;
    else # reverse move
      @inbounds move!(chain, idx, -ds[1], -ds[2]);
    end

    if time() - last_update > pargs["update-freq"]
      @info "elapsed: $(time() - start)";
      @info "step:    $step";
      last_update = time();
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
    @inbounds end_to_end_sum += chain.r[:];
    r2_sum += dot(chain.r, chain.r);
    chain_μ_sum += chain_μ(chain);
    Usum += chain.U;
  
  end

  @info "total time elapsed: $(time() - start)";
  @info "acceptance rate: $(num_accepted / pargs["num-steps"])";

  #for callback in callbacks
  #  callback(chain, pargs["num-steps"], false, true, pargs);
  #end
  close(outfile);

  return (end_to_end_sum, r2_sum, chain_μ_sum, Usum);

end

data = mcmc(pargs["num-steps"], pargs, callbacks);

println("<r>    =   $(data[1] / pargs["num-steps"])");
println("<r/nb> =   $(data[1] / 
                      (pargs["num-steps"]*pargs["mlen"]*pargs["num-monomers"])
                     )");
println("<r2>   =   $(data[2] / pargs["num-steps"])");
println("<p>    =   $(data[3] / pargs["num-steps"])");
println("<U>    =   $(data[4] / pargs["num-steps"])");
