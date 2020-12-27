const Accessor = Function;

abstract type Averager end;

# TODO: consider creating an averager which uses arbitrary precision arithmetic
mutable struct StandardAverager{T} <: Averager
  value::T;
  normalizer::Float64;
  accessor::Accessor;
end

function StandardAverager(accessor::Accessor, chain::EAPChain)
  return StandardAverager(
                          0.0 * accessor(chain),
                          0.0,
                          accessor
                         );
end

@inline get_avg(stdavg::StandardAverager) = stdavg.value / stdavg.normalizer;

function record!(stdavg::StandardAverager, chain::EAPChain)
  stdavg.value += stdavg.accessor(chain);
  stdavg.normalizer += 1;
end

function record!(stdavg::StandardAverager{FdVector}, chain::EAPChain)
  stdavg.value[:] += stdavg.accessor(chain);
  stdavg.normalizer += 1;
end

const WeightFunction = Function;

struct UmbrellaAverager{T} <: Averager
  stdavg::StandardAverager{T};
  weightFunction::WeightFunction;
end

function UmbrellaAverager(ac::Accessor, wf::WeightFunction, chain::EAPChain)
  return UmbrellaAverager(StandardAverager(ac, chain), wf);
end

@inline get_avg(ua::UmbrellaAverager) = get_avg(ua.stdavg);

function record!(ua::UmbrellaAverager, chain::EAPChain)
  w = ua.weightFunction(chain);
  ua.stdavg.value += ua.stdavg.accessor(chain) / w;
  ua.stdavg.normalizer += 1.0 / w;
end

function record!(ua::UmbrellaAverager{FdVector}, chain::EAPChain)
  w = ua.weightFunction(chain);
  ua.stdavg.value[:] += ua.stdavg.accessor(chain) / w;
  ua.stdavg.normalizer += 1.0 / w;
end

function weight_anti_dipole(chain::EAPChain)
  return exp(sum(chain.us) / chain.kT);
end
