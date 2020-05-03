function langevin_inv(y::Real)
  return (3*y - y/5 * (6*y^2 + y^4 - 2*y^6))/(1 - y^2);
end

function free_energy_langevin(k::Real, T::Real, N::Int, γ::Real)
  linvγ = langevin_inv(γ);
  return N*k*T*( log( linvγ / (4*pi*sinh(linvγ)) ) + γ*linvγ );
end

function force_langevin(k::Real, T::Real, b::Real, γ::Real)
  return k*T/b * langevin_inv(γ);
end

function forces(As::Vector, γs::Vector, N::Int)
  n = length(As)
  grads = zeros(n);
  grads[1] = (As[2]-As[1])/(γs[2]-γs[1]);
  grads[end] = (As[n]-As[n-1])/(γs[n]-γs[n-1]);
  for i=2:(n-1)
    grads[i] = (As[i+1] - As[i-1])/(γs[i+1] - γs[i-1]);
  end
  return grads / N;
end
