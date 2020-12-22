function μ(E0::Real, K1::Real, K2::Real, cϕ::Real, sϕ::Real, cθ::Real, sθ::Real)
  n̂i = n̂(cϕ, sϕ, cθ, sθ);
  return ((K1 - K2) * E0 * cθ * n̂i) + (K2 * [0.0; 0.0; E0]);
end
