function μ(E0::Real, K1::Real, K2::Real, cϕ::Real, sϕ::Real)
  n̂i = n̂(cϕ, sϕ);
  return ((K1 - K2) * E0 * cϕ * n̂i) + (K2 * [E0; 0.0]);
end
