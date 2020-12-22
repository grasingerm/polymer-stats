include(joinpath(@__DIR__, "types.jl"));

abstract type DipoleResponse end

(dr::DipoleResponse)(E0::Real) = begin; error("Not yet implemented"); E0; end

function μ_dielectric(E0::Real, K1::Real, K2::Real, cϕ::Real, sϕ::Real, cθ::Real, 
                      sθ::Real)
  n̂i = n̂(cϕ, sϕ, cθ, sθ);
  return ((K1 - K2) * E0 * cθ * n̂i) + (K2 * [0.0; 0.0; E0]);
end

struct DielectricResponse <: DipoleResponse
  K1::Real;
  K2::Real;
end

function (dr::DielectricResponse)(E0::Real, cϕ::Real, sϕ::Real, cθ::Real, 
                                  sθ::Real)
  return μ_dielectric(E0, dr.K1, dr.K2, cϕ, sϕ, cθ, sθ);
end

struct PolarResponse <: DipoleResponse
  M::FdMatrix;
end

function (pr::PolarResponse)(::Real, cϕ::Real, sϕ::Real, cθ::Real, sθ::Real)
  return pr.M * n̂(cϕ, sϕ, cθ, sθ);
end
