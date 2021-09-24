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
  logπ_prev::Float64;
  weight_function::WeightFunction;
end

logπ_chain(chain::EAPChain) = -chain.U / chain.kT + chain.Ω;

function logπ_chain(chain::EAPChain, log_wf::WeightFunction)
  return -chain.U / chain.kT + chain.Ω + log_wf(chain);
end

function Metropolis(chain::EAPChain, wf::WeightFunction)
  return Metropolis(logπ_chain(chain, wf), wf);
end

# alpha is a ratio of trial move probabilities (related to clustering)
function (metro::Metropolis)(chain::EAPChain, ϵ::Real; α::Real = 1.0)
  metro.weight_function(chain);
  logπ_curr = logπ_chain(chain, metro.weight_function) + log(α);
  if (logπ_curr >= metro.logπ_prev) || (ϵ < exp(logπ_curr - metro.logπ_prev))
    metro.logπ_prev = logπ_curr;
    return true;
  else
    return false;
  end
  return false;
end
