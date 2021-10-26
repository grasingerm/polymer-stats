include(joinpath(@__DIR__, "types.jl"));

struct NonInteractingEnergy <: Energy end

(::Energy)(::EAPChain) = error("Not yet implemented.");

@inline function (::NonInteractingEnergy)(chain::EAPChain)
  return (sum(chain.us) - dot(end_to_end(chain), [chain.Fx; chain.Fz]));
end

struct InteractingEnergy <: Energy end

@inline function (::InteractingEnergy)(chain::EAPChain)
  return (sum(chain.us) + U_interaction(chain)
          - dot(end_to_end(chain), [chain.Fx; chain.Fz]));
end

struct IsingEnergy <: Energy end

@inline function (::IsingEnergy)(chain::EAPChain)
  return (sum(chain.us) + U_Ising(chain)
          - dot(end_to_end(chain), [chain.Fx; chain.Fz]));
end
