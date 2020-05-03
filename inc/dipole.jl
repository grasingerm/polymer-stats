using LinearAlgebra;

function dipole(κ1::Real, κ2::Real, E0::Vector{Float64}, ϕ::Real, θ::Real)
  cϕ = cos(ϕ);
  sϕ = sin(ϕ);
  cθ = cos(θ);
  sθ = sin(θ);
  n = [cϕ*sθ; sϕ*sθ; cθ];

  #return (κ1 - κ2) * n * dot(E0, n) + κ2 * E0;
  χ = κ1 * (n * transpose(n)) + κ2 * (I - n * transpose(n));
  return χ*E0;
end

function polarization(xs::Vector)
  k2 = (dk > 0) ? dk : 0.0;
  k1 = (dk < 0) ? -dk : 0.0;
  I1 = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*dipole(k1, k2, [0.0; 0.0; E0], θs[1], θs[2])[1],
                 _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  I2 = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*dipole(k1, k2, [0.0; 0.0; E0], θs[1], θs[2])[2],
                 _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  I3 = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*dipole(k1, k2, [0.0; 0.0; E0], θs[1], θs[2])[3],
                 _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  return [I1[1], I2[1], I3[1]];
end

function polarization(xs::Vector, E0::Real, dk::Real)
  k2 = (dk > 0) ? dk : 0.0;
  k1 = (dk < 0) ? -dk : 0.0;
  ω0 = E0*E0*dk;
  I1 = pcubature(θs -> ρ(ω0, θs, xs)*sin(θs[2])*dipole(k1, k2, [0.0; 0.0; E0], θs[1], θs[2])[1],
                 _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  I2 = pcubature(θs -> ρ(ω0, θs, xs)*sin(θs[2])*dipole(k1, k2, [0.0; 0.0; E0], θs[1], θs[2])[2],
                 _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  I3 = pcubature(θs -> ρ(ω0, θs, xs)*sin(θs[2])*dipole(k1, k2, [0.0; 0.0; E0], θs[1], θs[2])[3],
                 _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  return [I1[1], I2[1], I3[1]];
end
