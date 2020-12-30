using Distributions;

include(joinpath(@__DIR__, "dipole_response.jl"));

ϕ_dist = Uniform(0.0, 2*π);
θ_dist = Uniform(0.0, π);

# the way that the energy functions are organized is a mess; I'm sorry
abstract type Energy end

mutable struct EAPChain
  b::Float64;
  E0::Float64;
  μ::DipoleResponse;
  UFunction::Energy;
  kT::Float64;
  Fz::Float64;
  Fx::Float64;
  ϕs::FdVector;
  cϕs::FdVector;
  sϕs::FdVector;
  θs::FdVector;
  cθs::FdVector;
  sθs::FdVector;
  Ω::Float64;
  μs::FdMatrix;
  us::FdVector;
  xs::FdMatrix;
  r::FdVector;
  U::Float64;
end

@inline n(chain::EAPChain) = length(chain.ϕs);

@inline n̂(cϕ, sϕ, cθ, sθ) = [cϕ * sθ; sϕ * sθ; cθ];

@inline n̂j(chain::EAPChain, idx::Int) = n̂(chain.cϕs[idx], chain.sϕs[idx],
                                          chain.cθs[idx], chain.sθs[idx]);

function update_xs!(chain::EAPChain, idx::Int=1)
  if idx == 1
    @inbounds chain.xs[:, idx] = chain.b / 2 * n̂j(chain, idx);
    idx += 1;
  end
  for i=idx:n(chain)
    @inbounds chain.xs[:, i] = (chain.xs[:, i-1] + 
                                chain.b / 2 * (n̂j(chain, i) + n̂j(chain, i-1))); 
  end
end

@inline u(E0::Real, μ::FdVector) = -1/2*E0*μ[3];

function EAPChain(pargs::Dict)
  ϕs = rand(ϕ_dist, pargs["num-monomers"]);
  θs = rand(θ_dist, pargs["num-monomers"]);

  dr = if pargs["chain-type"] == "dielectric"
    DielectricResponse(pargs["K1"], pargs["K2"]);
  elseif pargs["chain-type"] == "polar"
    PolarResponse(pargs["mu"] * [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]);
  else
    error("chain-type is not understood.");
  end

  ret =  EAPChain(
                  pargs["mlen"],
                  pargs["E0"],
                  dr,
                  if pargs["energy-type"] == "noninteracting"
                    NonInteractingEnergy();
                  elseif pargs["energy-type"] == "interacting"
                    InteractingEnergy();
                  else
                    error("energy-type is not understood.");
                  end,
                  pargs["kT"],
                  pargs["Fz"],
                  pargs["Fx"],
                  ϕs,
                  map(cos, ϕs),
                  map(sin, ϕs),
                  θs,
                  map(cos, θs),
                  map(sin, θs),
                  log(prod(map(sin, θs))), # store logarithm of solid angle instead
                  zeros(3, pargs["num-monomers"]),
                  zeros(pargs["num-monomers"]),
                  zeros(3, pargs["num-monomers"]),
                  zeros(3),
                  0.0
                 );
  @simd for i=1:n(ret)
    @inbounds ret.μs[:, i] = dr(ret.E0, ret.cϕs[i], ret.sϕs[i], 
                                ret.cθs[i], ret.sθs[i]);
  end
  ret.us[:] = map(i -> u(ret.E0, view(ret.μs, :, i)), 1:n(ret));
  update_xs!(ret);
  ret.r[:] = end_to_end(ret);
  ret.U = U(ret);
  return ret;
end

function EAPChain(chain::EAPChain)
  return EAPChain(
                  chain.b,
                  chain.E0,
                  chain.μ,
                  chain.UFunction,
                  chain.kT,
                  chain.Fz,
                  chain.Fx,
                  chain.ϕs[:],
                  chain.cϕs[:],
                  chain.sϕs[:],
                  chain.θs[:],
                  chain.cθs[:],
                  chain.sθs[:],
                  chain.Ω,
                  chain.μs[:, :],
                  chain.us[:],
                  chain.xs[:, :],
                  chain.r[:],
                  chain.U
                 );
end

# TODO: consider rewriting this to make better use of past calculations
# this is probably by far the slowest calculation *shrug emoji*
function U_interaction(chain::EAPChain)
  U = 0.0;
  @inbounds for i=1:n(chain)
    @simd for j=i+1:n(chain) # once debugged, use inbounds
      r = chain.xs[:, i] - chain.xs[:, j];
      r2 = dot(r, r);
      rmag = sqrt(r2);
      r̂ = r / rmag;
      r3 = r2*rmag;
      μi = view(chain.μs, :, i);
      μj = view(chain.μs, :, i);
      U += (dot(μi, μj) - 3*dot(μi, r̂)*dot(μj, r̂)) / (4*π);
    end
  end
  return U;
end

function move!(chain::EAPChain, idx::Int, dϕ::Real, dθ::Real)
  @inbounds begin
    chain.ϕs[idx] += dϕ;
    chain.cϕs[idx] = cos(chain.ϕs[idx]);
    chain.sϕs[idx] = sin(chain.ϕs[idx]);

    chain.θs[idx] = min(π, max(0.0, chain.θs[idx]+dθ));
    sθ = sin(chain.θs[idx]);
    chain.Ω += log(sθ / chain.sθs[idx]); # update solid angle
    chain.cθs[idx] = cos(chain.θs[idx]);
    chain.sθs[idx] = sθ;

    # the order of these updates is important
    chain.μs[:, idx] = chain.μ(chain.E0, chain.cϕs[idx], chain.sϕs[idx], 
                               chain.cθs[idx], chain.sθs[idx]);
    chain.us[idx] = u(chain.E0, view(chain.μs, :, idx));
    update_xs!(chain, idx);
    chain.r[:] = end_to_end(chain);
    chain.U = U(chain);
  end
  return (dϕ, dθ);
end

@inline end_to_end(chain::EAPChain) = @inbounds (chain.xs[:, end] + 
                                                 chain.b/2.0*n̂j(chain, n(chain)));

@inline function chain_μ(chain::EAPChain)
  @inbounds begin;
    sum(map(i -> chain.μs[:, i], 1:n(chain)));
  end
end

include(joinpath(@__DIR__, "energy.jl"));
@inline U(chain::EAPChain) = chain.UFunction(chain); # forward energy call to chain
