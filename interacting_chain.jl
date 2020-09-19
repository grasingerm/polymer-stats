using ArgParse;
using Distributions;
using LinearAlgebra;
using Logging;

CnMatrix = AbstractMatrix{Int};
FdMatrix = AbstractMatrix{Float64};
CnVector = AbstractVector{Int};
FdVector = AbstractVector{Float64};

ϕ_dist = Uniform(0.0, 2*π);
θ_dist = Uniform(0.0, π);

struct EAPChain
  b::Real;
  E0::Real;
  K1::Real;
  K2::Real;
  kT::Real;
  Fz::Real;
  Fx::Real;
  ϕs::FdVector;
  cϕs::FdVector;
  sϕs::FdVector;
  θs::FdVector;
  cθs::FdVector;
  sθs::FdVector;
  μs::FdMatrix;
  us::FdVector;
  xs::FdMatrix;
end

@inline n(chain::EAPChain) = length(chain.ϕs);

@inline n̂(cϕ, sϕ, cθ, sθ) = [cϕ * sθ; sϕ * sθ; cθ];

function μ(E0::Real, K1::Real, K2::Real, cϕ::Real, sϕ::Real, cθ::Real, sθ::Real)
  n̂i = n̂(cϕ, sϕ, cθ, sθ);
  return ((K1 - K2) * E0 * cθ * n̂i) + (K2 * [0.0; 0.0; E0]);
end

@inline n̂j(chain::EAPChain, idx::Int) = n̂(chain.cϕs[idx], chain.sϕs[idx],
                                          chain.cθs[idx], chain.sθs[idx]);

function update_xs!(chain::EAPChain, idx::Int=1)
  if idx == 1
    chain.xs[:, idx] = chain.b / 2 * n̂j(chain, idx);
    idx += 1;
  end
  for i=idx:n(chain)
    chain.xs[:, i] = (chain.xs[:, i-1] + 
                      chain.b / 2 * (n̂j(chain, i) + n̂j(chain, i-1))); 
  end
end

@inline u(E0::Real, μ::FdVector) = -1/2*E0*μ[3];

function EAPChain(pargs::Dict)
  ϕs = rand(ϕ_dist, pargs["num-monomers"]);
  θs = rand(θ_dist, pargs["num-monomers"]);

  ret =  EAPChain(
                  pargs["mlen"],
                  pargs["E0"],
                  pargs["K1"],
                  pargs["K2"],
                  pargs["kT"],
                  pargs["Fz"],
                  pargs["Fx"],
                  ϕs,
                  map(cos, ϕs),
                  map(sin, ϕs),
                  θs,
                  map(cos, θs),
                  map(sin, θs),
                  zeros(3, pargs["num-monomers"]),
                  zeros(pargs["num-monomers"]),
                  zeros(3, pargs["num-monomers"]),
                 );
  for i=1:n(ret)
    ret.μs[:, i] = μ(ret.E0, ret.K1, ret.K2, 
                     ret.cϕs[i], ret.sϕs[i], ret.cθs[i], ret.sθs[i]);
  end
  ret.us[:] = map(i -> u(ret.E0, view(ret.μs, :, i)), 1:n(ret));
  update_xs!(ret);
  return ret;
end

# TODO: consider rewriting this to make better use of past calculations
# this is probably by far the slowest calculation *shrug emoji*
function U_interaction(chain::EAPChain)
  U = 0.0;
  for i=1:n(chain), j=i+1:n(chain) # once debugged, use inbounds
    r = chain.xs[:, i] - chain.xs[:, j];
    r2 = dot(r, r);
    rmag = sqrt(r2);
    r̂ = r / rmag;
    r3 = r2*rmag;
    μi = view(chain.μs, :, i);
    μj = view(chain.μs, :, i);
    U += (dot(μi, μj) - 3*dot(μi, r̂)*dot(μj, r̂)) / (4*π);
  end
  return U;
end

function move!(chain::EAPChain, idx::Int, dϕ::Real, dθ::Real)
  ϕprev = chain.ϕs[idx];
  chain.ϕs[idx] += dϕ;
  if chain.ϕs[idx] < 0
    chain.ϕs[idx] += 2*π;
  elseif chain.ϕs[idx] >= 2*π
    chain.ϕs[idx] -= 2*π;
  end
  chain.cϕs[idx] = cos(chain.ϕs[idx]);
  chain.sϕs[idx] = sin(chain.ϕs[idx]);

  θprev = chain.θs[idx]; 
  chain.θs[idx] += dθ;
  if chain.θs[idx] < 0
    chain.θs[idx] = 0;
  elseif chain.θs[idx] >= π
    chain.θs[idx] = π;
  end
  chain.cθs[idx] = cos(chain.θs[idx]);
  chain.sθs[idx] = sin(chain.θs[idx]);

  chain.μs[:, idx] = μ(chain.E0, chain.K1, chain.K2,
                       chain.cϕs[idx], chain.sϕs[idx], 
                       chain.cθs[idx], chain.sθs[idx]);
  chain.us[idx] = u(chain.E0, view(chain.μs, :, idx));
  update_xs!(chain, idx);
  return (chain.ϕs[idx] - ϕprev, chain.θs[idx] - θprev);
end

@inline end_to_end(chain::EAPChain) = (chain.xs[:, end] + 
                                       chain.b/2*n̂j(chain, n(chain)));

@inline function U(chain::EAPChain)
  return (sum(chain.us) + U_interaction(chain)
          - dot(end_to_end(chain), [chain.Fx; 0.0; chain.Fz]));
end

@inline chain_μ(chain::EAPChain) = sum(map(i -> chain.μs[:, i], 1:n(chain)));

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
    default = convert(Int, 1e6)
  "--phi-step", "-p"
    help = "maximum ϕ step length"
    arg_type = Float64;
    default = π / 16;
  "--theta-step", "-q"
    help = "maximum θ step length"
    arg_type = Float64;
    default = π / 32;
  "--update-freq", "-U"
    help = "update frequency (seconds)"
    arg_type = Float64;
    default = 15.0;
  "--verbose", "-v"
    help = "verbosity level: 0-nothing, 1-errors, 2-warnings, 3-info"
    arg_type = Int
    default = 2
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

function mcmc(nsteps, pargs)
  dϕ_dist = Uniform(-pargs["phi-step"], pargs["phi-step"]);
  dθ_dist = Uniform(-pargs["theta-step"], pargs["theta-step"]);
  chain = EAPChain(pargs);
  Uprev = U(chain);
  r = end_to_end(chain);
  end_to_end_sum = r[:];
  r2_sum = dot(r, r);
  chain_μ_sum = chain_μ(chain);
  Usum = Uprev;

  start = time();
  last_update = start;
  num_accepted = 0;
  for step=1:nsteps
    dϕ = rand(dϕ_dist);
    dθ = rand(dϕ_dist);
    idx = rand(1:n(chain));
    ds = move!(chain, idx, dϕ, dθ);
    Ucurr = U(chain);
    if metropolis_acc(chain.kT, Ucurr - Uprev, rand())
      Uprev = Ucurr;
      num_accepted += 1;
    else # reverse move
      Ucurr = Uprev;
      move!(chain, idx, -ds[1], -ds[2]);
    end

    if time() - last_update > pargs["update-freq"]
      @info "elapsed: $(time() - start)";
      @info "step:    $step";
      last_update = time();
    end

    r = end_to_end(chain);
    end_to_end_sum += r[:];
    r2_sum += dot(r, r);
    chain_μ_sum += chain_μ(chain);
    Usum += Ucurr;
  end

  @info "total time elapsed: $(time() - start)";
  @info "acceptance rate: $(num_accepted / pargs["num-steps"])";
  return (end_to_end_sum, r2_sum, chain_μ_sum, Usum);

end

@show data = mcmc(pargs["num-steps"], pargs);

println("<r>    =   $(data[1] / pargs["num-steps"])");
println("<r/nb> =   $(data[1] / 
                      (pargs["num-steps"]*pargs["mlen"]*pargs["num-monomers"])
                     )");
println("<r2>   =   $(data[2] / pargs["num-steps"])");
println("<p>    =   $(data[3] / pargs["num-steps"])");
println("<U>    =   $(data[4] / pargs["num-steps"])");
