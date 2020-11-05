using Distributions;

CnMatrix = AbstractMatrix{Int};
FdMatrix = AbstractMatrix{Float64};
CnVector = AbstractVector{Int};
FdVector = AbstractVector{Float64};

ϕ_dist = Uniform(0.0, 2*π);
θ_dist = Uniform(0.0, π);

mutable struct EAPChain
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
  r::FdVector;
  U::Real;
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
                  zeros(3),
                  0.0
                 );
  @simd for i=1:n(ret)
    @inbounds ret.μs[:, i] = μ(ret.E0, ret.K1, ret.K2, 
                               ret.cϕs[i], ret.sϕs[i], ret.cθs[i], ret.sθs[i]);
  end
  ret.us[:] = map(i -> u(ret.E0, view(ret.μs, :, i)), 1:n(ret));
  update_xs!(ret);
  ret.r[:] = end_to_end(ret);
  ret.U = U(ret);
  return ret;
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
    chain.r[:] = end_to_end(chain);
  end
  return (chain.ϕs[idx] - ϕprev, chain.θs[idx] - θprev);
end

@inline end_to_end(chain::EAPChain) = @inbounds (chain.xs[:, end] + 
                                                 chain.b/2*n̂j(chain, n(chain)));

@inline function chain_μ(chain::EAPChain)
  @inbounds begin;
    sum(map(i -> chain.μs[:, i], 1:n(chain)));
  end
end
