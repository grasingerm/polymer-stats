@inline function metropolis_acc(kT::Real, dU::Real, sθa::Real, sθb::Real, ϵ::Real)
  return ϵ <= (exp(-dU / kT) * sθb / sθa);
end
  
function kawasaki_acc(kT::Real, dU::Real, ϵ::Real)
  boltz = exp(-dU / (2*kT));
  anti_boltz = exp(dU / (2*kT));
  return ( ϵ <= (boltz / (boltz + anti_boltz)) );
end

abstract type Acceptor end;

mutable struct Metropolis <: Acceptor
  π_prev::Float64;
  weight_function::Function;
end

π_chain(chain::EAPChain) = exp(-chain.U / chain.kT) * chain.Ω;

function Metropolis(chain::EAPChain, wf::Function)
  return Metropolis(wf(chain) * π_chain(chain), wf);
end

function (metro::Metropolis)(chain::EAPChain, ϵ::Real)
  metro.weight_function(chain);
  π_curr = metro.weight_function(chain) * π_chain(chain);
  if ϵ < π_curr / metro.π_prev
    metro.π_prev = π_curr;
    return true;
  else
    return false;
  end
  return false;
end
