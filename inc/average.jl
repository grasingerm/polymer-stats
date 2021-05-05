

const Accessor = Function;

abstract type Averager end;

# TODO: consider creating an averager which uses arbitrary precision arithmetic
mutable struct StandardAverager{T,N} <: Averager
  value::T;
  normalizer::N;
  accessor::Accessor;
end

function StandardAverager(accessor::Accessor, chain::EAPChain)
  return StandardAverager(
                          0.0 * accessor(chain),
                          0.0,
                          accessor
                         );
end

function StandardAverager(dt::DataType, accessor::Accessor, chain::EAPChain)
  return StandardAverager(
                          dt(0.0 * accessor(chain)),
                          dt(0.0),
                          accessor
                         );
end

function StandardAverager(dt1::DataType, dt2::DataType, accessor::Accessor, chain::EAPChain)
  return StandardAverager(
                          dt2(0.0 * accessor(chain)),
                          dt1(0.0),
                          accessor
                         );
end

@inline get_avg(stdavg::StandardAverager) = stdavg.value / stdavg.normalizer;

function record!(stdavg::StandardAverager, chain::EAPChain)
  stdavg.value += stdavg.accessor(chain);
  stdavg.normalizer += 1;
end

function record!(stdavg::StandardAverager{<:Vector}, chain::EAPChain)
  stdavg.value[:] += stdavg.accessor(chain);
  stdavg.normalizer += 1;
end

abstract type WeightFunction end;

struct UmbrellaAverager{T,N} <: Averager
  stdavg::StandardAverager{T,N};
  weightFunction::WeightFunction;
end

function UmbrellaAverager(ac::Accessor, wf::WeightFunction, chain::EAPChain)
  return UmbrellaAverager(StandardAverager(ac, chain), wf);
end

@inline get_avg(ua::UmbrellaAverager) = get_avg(ua.stdavg);

function record!(ua::UmbrellaAverager, chain::EAPChain)
  expw = exp(ua.weightFunction(chain));
  ua.stdavg.value += ua.stdavg.accessor(chain) / expw;
  ua.stdavg.normalizer += 1.0 / expw;
end

function record!(ua::UmbrellaAverager{Vector{<:Real}}, chain::EAPChain)
  expw = exp(ua.weightFunction(chain));
  ua.stdavg.value[:] += ua.stdavg.accessor(chain) / expw;
  ua.stdavg.normalizer += 1.0 / expw;
end

function record!(ua::UmbrellaAverager{Dec128}, chain::EAPChain)
  expw = exp(Dec128(ua.weightFunction(chain)));
  ua.stdavg.value += ua.stdavg.accessor(chain) / expw;
  ua.stdavg.normalizer += 1.0 / expw;
end

function record!(ua::UmbrellaAverager{Vector{Dec128}}, chain::EAPChain)
  expw = exp(Dec128(ua.weightFunction(chain)));
  ua.stdavg.value[:] += ua.stdavg.accessor(chain) / expw;
  ua.stdavg.normalizer += 1.0 / expw;
end

function record!(ua::UmbrellaAverager{BigFloat}, chain::EAPChain)
  expw = exp(BigFloat(ua.weightFunction(chain)));
  ua.stdavg.value += ua.stdavg.accessor(chain) / expw;
  ua.stdavg.normalizer += 1.0 / expw;
end

function record!(ua::UmbrellaAverager{Vector{BigFloat}}, chain::EAPChain)
  expw = exp(BigFloat(ua.weightFunction(chain)));
  ua.stdavg.value[:] += ua.stdavg.accessor(chain) / expw;
  ua.stdavg.normalizer += 1.0 / expw;
end

struct WeightlessFunction <: WeightFunction
end

(wf::WeightlessFunction)(::EAPChain) = 1.0;

struct AntiDipoleWeightFunction <: WeightFunction
  log_gauge::Float64;
end

# forward call
@inline AntiDipoleWeightFunction(chain::EAPChain) = AntiDipoleWeightFunction(chain, chain.μ);
function AntiDipoleWeightFunction(chain::EAPChain, dr::DielectricResponse)
  return AntiDipoleWeightFunction(-(dr.K1 + 2 * dr.K2) * chain.E0*chain.E0 *
                                   n(chain) / (3 * chain.kT) + chain.Ω);
end

function AntiDipoleWeightFunction(chain::EAPChain, dr::PolarResponse)
  return AntiDipoleWeightFunction(-eigvals(dr.M)[end] * chain.E0 *
                                   n(chain) / (3 * chain.kT) + chain.Ω);
end

function (wf::AntiDipoleWeightFunction)(chain::EAPChain)
  return (sum(chain.us) / chain.kT * 
          (0.2 + 0.8*(exp(-(chain.Fx^2+chain.Fz^2)/chain.kT))) 
          - wf.log_gauge);
end
