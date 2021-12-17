@inline function U(chain::EAPChain)
  return (sum(chain.us) - dot(chain.r, [chain.Fx; 0.0; chain.Fz]));
end
