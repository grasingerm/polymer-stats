@inline function U(chain::EAPChain)
  return (sum(chain.us) + U_interaction(chain)
          - dot(end_to_end(chain), [chain.Fx; 0.0; chain.Fz]));
end
