using Distributions;
using NLsolve;
using ForwardDiff;

include(joinpath(@__DIR__, "dipole_response.jl"));

ϕ_dist = Uniform(0.0, 2*π);
θ_dist = Uniform(0.0, π);

# the way that the energy functions are organized is a mess; I'm sorry
# also... not sorry
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
  n̂s::FdMatrix;
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

function update_xs!(chain::EAPChain, idx::Int)
  @warn "This implementation of update_xs! is generally slower than the vectorized implementation";
  if idx == 1
    @inbounds chain.xs[:, idx] = chain.b / 2 * n̂j(chain, idx);
    idx += 1;
  end
  for i=idx:n(chain)
    @inbounds chain.xs[:, i] = (chain.xs[:, i-1] + 
                                chain.b / 2 * (n̂j(chain, i) + n̂j(chain, i-1))); 
  end
end

@inline function update_xs!(chain::EAPChain)
  chain.xs[:, :] = chain.b*(cumsum(chain.n̂s, dims=2) - 0.5*chain.n̂s);
end

#=
function update_xs!(chain::EAPChain)
  @show chain.n̂s
  @show cumsum(chain.n̂s, dims=2);
  @show (cumsum(chain.n̂s, dims=2) - 0.5*chain.n̂s);
  idx = rand(1:n(chain));
  @assert(dot(chain.n̂s[:, idx], chain.n̂s[:, idx]) ≈ 1.0);
  chain.xs[:, :] = chain.b*(cumsum(chain.n̂s, dims=2) - 0.5*chain.n̂s);
end
=#

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
                  elseif pargs["energy-type"] == "Ising"
                    IsingEnergy();
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
                  zeros(3, pargs["num-monomers"]),
                  log(prod(map(sin, θs))), # store logarithm of solid angle instead
                  zeros(3, pargs["num-monomers"]),
                  zeros(pargs["num-monomers"]),
                  zeros(3, pargs["num-monomers"]),
                  zeros(3),
                  0.0
                 );
  ret.n̂s[:, :] = hcat(map(j -> n̂j(ret, j), 1:n(ret))...);
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
                  chain.n̂s[:, :],
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
      μj = view(chain.μs, :, j);
      U += (dot(μi, μj) - 3*dot(μi, r̂)*dot(μj, r̂)) / (4*π*r3);
    end
  end
  return U;
end

# TODO: consider rewriting this to make better use of past calculations
# this is probably by far the slowest calculation *shrug emoji*
function U_Ising(chain::EAPChain)
  U = 0.0;
  @inbounds for i=1:n(chain)-1
    r = chain.xs[:, i] - chain.xs[:, i+1];
    r2 = dot(r, r);
    rmag = sqrt(r2);
    r̂ = r / rmag;
    r3 = r2*rmag;
    μi = view(chain.μs, :, i);
    μj = view(chain.μs, :, i+1);
    U += (dot(μi, μj) - 3*dot(μi, r̂)*dot(μj, r̂)) / (4*π*r3);
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
    chain.n̂s[:, idx] = n̂j(chain, idx);
    chain.μs[:, idx] = chain.μ(chain.E0, chain.cϕs[idx], chain.sϕs[idx], 
                               chain.cθs[idx], chain.sθs[idx]);
    chain.us[idx] = u(chain.E0, view(chain.μs, :, idx));
    update_xs!(chain);
    chain.r[:] = end_to_end(chain);
    chain.U = U(chain);
  end
  return true;
end

function flip_n!(chain::EAPChain, idx::Int)
  move!(chain, idx, π, π - 2*chain.θs[idx])
end

function refl_n!(chain::EAPChain, idx::Int)
  move!(chain, idx, 0, π - 2*chain.θs[idx])
end

pflip_linear(ip::Real) = (1 + ip) / 2;

function cluster_flip!(chain::EAPChain, idx::Int; 
                       pflip::Function = pflip_linear,
                       ϵflip::Real = 1.0,
                       ηerr::Real = 1e-10,
                       flip_f!::Function = refl_n!
                      )
  
  # grow the cluster to the right
  upper_p = NaN;
  upper_idx = idx;
  while true
    if upper_idx >= n(chain)
      upper_p = 0.0;
      break;
    end
    upper_p = pflip_linear(dot(n̂j(chain, upper_idx), n̂j(chain, upper_idx+1)));
    if rand() <= upper_p
      upper_idx += 1
    else
      break;
    end
  end
  @assert(!isnan(upper_p) && upper_idx <= n(chain));

  # grow the cluster to the left
  lower_p = NaN;
  lower_idx = idx;
  while true
    if lower_idx <= 1
      lower_p = 0.0;
      break;
    end
    lower_p = pflip_linear(dot(n̂j(chain, lower_idx), n̂j(chain, lower_idx-1)));
    if rand() <= lower_p
      lower_idx -= 1
    else
      break;
    end
  end
  @assert(!isnan(lower_p) && lower_idx >= 1);

  if rand() <= ϵflip
    prev_n̂s = copy(chain.n̂s[:, lower_idx:upper_idx]);
    for i in lower_idx:upper_idx
      flip_f!(chain, i);
    end

    if norm(prev_n̂s + chain.n̂s[:, lower_idx:upper_idx], Inf) > ηerr
      println("prev");
      display(prev_n̂s)
      println();
      println("curr");
      display(chain.n̂s[:, lower_idx:upper_idx])
      println();
      @show(norm(prev_n̂s + chain.n̂s[:, lower_idx:upper_idx], Inf))
    end
    @assert(norm(prev_n̂s + chain.n̂s[:, lower_idx:upper_idx], Inf) < ηerr);

    new_upper_p = if upper_idx < n(chain)
      pflip_linear(dot(n̂j(chain, upper_idx), n̂j(chain, upper_idx+1)));
    else
      0
    end
    new_lower_p = if lower_idx > 1
      pflip_linear(dot(n̂j(chain, lower_idx), n̂j(chain, lower_idx-1)));
    else
      0
    end
    return ( 
             ((1 - new_upper_p)*(1 - new_lower_p)) /
             ((1 - upper_p)*(1 - lower_p))
            );
  else
    return 1.0;
  end

end

function move!(chain::EAPChain, idx::Int, dϕ::Real, dθ::Real, r0::AbstractVector;
               frac_mv::Real = 0.15, acc_tol::Real = 1e-1)

  @assert(frac_mv > 0.0 && frac_mv <= 1.0, 
          "Fraction of chain to move must be between 0.0 and 1.0");

  L0 = chain.b * n(chain);
  
  trial_idxs = if (idx / n(chain)) >= (1 - frac_mv)
    max(1, round(Int, (1 - frac_mv)*n(chain))):n(chain);
  elseif (idx / n(chain)) <= frac_mv
    1:(min(n(chain), round(Int, frac_mv*n(chain))));
  else
    (max(1, round(Int, idx - frac_mv*n(chain)/2.0))):round(Int, idx + frac_mv*n(chain)/2.0)
  end

  r1 = chain.b * ([
                   sum(chain.cϕs[trial_idxs] .* chain.sθs[trial_idxs]);
                   sum(chain.sϕs[trial_idxs] .* chain.sθs[trial_idxs]);
                   sum(chain.cθs[trial_idxs])
                  ]);

  # internal lambda function for solving nonlinear system
  f! = (F::Vector, x::Vector) -> begin;
    n = convert(Int, length(x) / 2);
    cϕs = map(cos, x[1:n]);
    sϕs = map(sin, x[1:n]);
    sθs = map(sin, x[(n+1):(2*n)]);
    F[1] = chain.b * sum(cϕs .* sθs) - r1[1];
    F[2] = chain.b * sum(sϕs .* sθs) - r1[2];
    F[3] = chain.b * sum(map(cos, x[(n+1):(2*n)])) - r1[3];
  end

  #=
  j! = (J::Matrix, x::Vector) -> begin;
    n = convert(Int, length(x) / 2);
    cϕs = map(cos, x[1:n]);
    sϕs = map(sin, x[1:n]);
    cθs = map(cos, x[(n+1):(2*n)]);
    sθs = map(sin, x[(n+1):(2*n)]);
    J[:, :] = chain.b * hcat(map(i -> [cϕs[i]*cθs[i] -sϕs[i]*sθs[i];
                                       sϕs[i]*cθs[i]  cϕs[i]*sθs[i];
                                       -sθs[i]        0.0;], 1:n)...);
  end
  =#

  # trial move initial guess
  x0 = vcat(chain.ϕs[trial_idxs], chain.θs[trial_idxs]);
  local_idx = findnext(i -> i==idx, collect(trial_idxs), 1);
  x0[local_idx] += dϕ;
  x0[2*local_idx] += dθ;

  # solve nonlinear system to ensure end to end vector constrain is satisfied
  result = nlsolve(f!, x0);
  if result.residual_norm / L0 > acc_tol
    @warn "residual norm is high", result.residual_norm / L0;
    return false;
  end

  # update chain based on result
  m = length(trial_idxs);
  chain.ϕs[trial_idxs] = result.zero[1:m];
  chain.θs[trial_idxs] = result.zero[(m+1):(2*m)];
  chain.cϕs[trial_idxs] = map(cos, chain.ϕs[trial_idxs]);
  chain.sϕs[trial_idxs] = map(sin, chain.ϕs[trial_idxs]);
  chain.cθs[trial_idxs] = map(cos, chain.θs[trial_idxs]);
  chain.sθs[trial_idxs] = map(sin, chain.θs[trial_idxs]);

  # the order of these updates is important
  chain.n̂s[:, trial_idxs] = hcat(map(idx -> n̂j(chain, idx), trial_idxs)...);
  chain.μs[:, trial_idxs] = hcat(map(idx -> chain.μ(chain.E0, chain.cϕs[idx], 
                                                    chain.sϕs[idx], chain.cθs[idx], 
                                                    chain.sθs[idx]), 
                                     trial_idxs)...);
  chain.us[trial_idxs] = map(idx -> u(chain.E0, view(chain.μs, :, idx)), 
                             trial_idxs);
  update_xs!(chain);
  chain.r[:] = end_to_end(chain);
  chain.U = U(chain);
  
  return true;
end

@inline end_to_end(chain::EAPChain) = @inbounds (chain.xs[:, end] + 
                                                 chain.b/2.0*n̂j(chain, n(chain)));

@inline chain_μ(chain::EAPChain) = reshape(sum(chain.μs, dims=2), 3);

include(joinpath(@__DIR__, "energy.jl"));
@inline U(chain::EAPChain) = chain.UFunction(chain); # forward energy call to chain
