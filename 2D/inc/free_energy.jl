function free_energy(xs::Vector)
  I1 = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*u(θs),
                 _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  I2 = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*log(ρ(θs, xs)),
                 _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  return I1[1] + kT*I2[1] - kT*N*log(N);
end

function free_energy_smallomega(x::Vector)
  c, λ, α = x[:]; 
  I1 = pcubature(θs -> c*(1 - ω(θs) + α*cos(θs[1])*sin(θs[2]))*exp(λ*cos(θs[2]))*sin(θs[2])*u(θs),
                 _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  I2 = pcubature(θs -> c*(1 - ω(θs) + α*cos(θs[1])*sin(θs[2]))*exp(λ*cos(θs[2]))*sin(θs[2])*log(ρ(θs, xs)),
                 _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  return I1[1] + kT*I2[1] - kT*N*log(N);
end

function free_energy_smalllambda(x::Vector)
  c, λ, α = x[:]; 
  I1 = pcubature(θs -> c*(1 + λ*cos(θs[2]) + α*cos(θs[1])*sin(θs[2]))*exp(-ω(θs))*sin(θs[2])*u(θs),
                 _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
                 I2 = pcubature(θs -> c*(1 + λ*cos(θs[2]) + α*cos(θs[1])*sin(θs[2]))*exp(-ω(θs))*sin(θs[2])*log(ρ(θs, xs)),
                 _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  return I1[1] + kT*I2[1] - kT*N*log(N);
end
