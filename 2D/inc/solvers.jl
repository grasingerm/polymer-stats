include("langevin.jl");

using NLopt;
using Cubature;
using LinearAlgebra;

 _ALGO_DICT_ = Dict("PRAXIS" => :LN_PRAXIS, "SBPLX" => :LN_SBPLX, 
                         "ISRES" => :GN_ISRES, "Newton" => :NEWTON, 
                         "SLSQP" => :LD_SLSQP, "LBFGS" => :LD_LBFGS,
                         "AUGLAG" => :AUGLAG, "Anneal" => :ANNEAL);

function algostr_to_symbol(algostr::AbstractString)
  if algostr in keys(_ALGO_DICT_)
    return _ALGO_DICT_[algostr];
  else
    error("Do not understand the algorithm: \"" * algostr * "\"");
  end
  return :AJKLLKJALDSKJFAS
end

function parse_algostr(algostr::AbstractString)
  if startswith(algostr, "try")
    steps = map(strip, split(algostr[4:end], "then"));
    return map(parse_algostr, steps);
  end
  return map(algostr_to_symbol, map(strip, split(algostr, ",")));
end

function parse_max_iters(max_iters_str::AbstractString)
  if startswith(max_iters_str, "try")
    steps = map(strip, split(max_iters_str[4:end], "then"));
    return map(parse_max_iters, steps);
  end
  return map(x -> parse(Int, x), map(strip, split(max_iters_str, ",")));
end

function initial_guess(γx::Real, γz::Real)
  γ = sqrt(γx^2 + γz^2);
  return [(γ > 0) ? N*langevin_inv(γ)*csch(langevin_inv(γ))/(4*pi) : 0;
          (γz > 0) ? γz / γ * langevin_inv(γ) : 0; 
          (γx > 0) ? γx / γ * langevin_inv(γ) : 0];
end

function initial_guess(γ::Real)
  return [γ > 0 ? N*langevin_inv(γ)*csch(langevin_inv(γ))/(4*pi) : 0;
          γ > 0 ? langevin_inv(γ) : 0; 
          0];
end

function initial_guess(ω::Real, γs::Tuple{Real,Real})
  sqrtω = sqrt(abs(ω));
  erfsω = ω > 0 ? erf(sqrtω) : erfi(sqrtω);
  expω = exp(ω);
  sqrtpi = sqrt(pi);
  γx, γz = γs;
  return [
          N * sqrtω / (2*pi*sqrtpi*erfsω);
          (2*expω*sqrtpi*γz*ω*erfω) / (-2*sqrtω + expω*sqrtpi*erfsω);
          (4*expω*sqrtpi*γx*ω*erfω) / (2*sqrtω + expω*sqrtpi*(-1+2*ω)*erfsω)
         ];
end

function initial_guess(ωs::Tuple{Real,Real}, γ::Real)
  ωx, ωz = ωs;
  ω = ωx + ωz;
  sqrtωx = sqrt(abs(ωx));
  sqrtωz = sqrt(abs(ωz));
  sqrtω = sqrt(abs(ω));
  erfsω = ω > 0 ? erf(sqrtω) : erfi(sqrtω);
  expω = exp(ω);
  sqrtpi = sqrt(pi);
  return [
          N * sqrtω / (2*pi*sqrtpi*erfsω);
          (2*expω*sqrtpi*sqrtωz*sqrtω*erfsω) / (-2*sqrtω + expω*sqrtpi*erfsω);
          (4*expω*sqrtpi*sqrtωx*sqrtω*erfsω) / (2*sqrtω + expω*sqrtpi*(-1+2*ω)*erfsω)
         ];
end

 _RHS_ECOORD_ = Float64[N; γz*N; γx*N];

function residuals(Is::Vector)
  global _RHS_;
  map(i -> Is[i][1] - _RHS_[i], 1:3)
end

ω(θs::Vector) = ω(θs[1], θs[2]);

ρ(ϕ::Real, θ::Real, c::Real, λ::Real, α::Real) = c * exp(-ω(ϕ, θ) + λ*cos(θ) 
                                                         + α*cos(ϕ)*sin(θ));
ρ(θs::Vector, xs::Vector) = ρ(θs[1], θs[2], xs[1], xs[2], xs[3]);
u(ϕ::Real, θ::Real) = kT * ω(ϕ, θ) + _u0_;
u(θs::Vector) = kT * ω(θs[1], θs[2]) + _u0_;

integrand_functions(xs::Vector) = (
  [
    (θs::Vector) -> ρ(θs, xs)*sin(θs[2]),
    (θs::Vector) -> ρ(θs, xs)*sin(θs[2])*cos(θs[2]),
    (θs::Vector) -> ρ(θs, xs)*sin(θs[2])*sin(θs[2])*cos(θs[1])
  ]
);

 _LOWER_BOUNDS_ = [0.0, 0.0];
 _UPPER_BOUNDS_ = [2*pi, pi];

function integrate_sys(xs::Vector)
  return map(integrand -> pcubature(integrand, _LOWER_BOUNDS_, _UPPER_BOUNDS_; 
                                    reltol=RELTOL_INT, maxevals=MAXEVALS_INT),
             integrand_functions(xs));
end

function gradc(xs::Vector)
  global RELTOL_INT;
  global MAXEVALS_INT;
  
  grad = zeros(3);
  grad[1], = pcubature(θs -> ρ(θs, xs)/xs[1]*sin(θs[2]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  grad[2], = pcubature(θs -> ρ(θs, xs)/xs[1]*sin(θs[2])*cos(θs[2]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  grad[3], = pcubature(θs -> ρ(θs, xs)/xs[1]*sin(θs[2])*sin(θs[2])*cos(θs[1]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);

  return grad;
end

function gradλ(xs::Vector)
  global RELTOL_INT;
  global MAXEVALS_INT;
  
  grad = zeros(3);
  grad[1], = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*cos(θs[2]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  grad[2], = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*cos(θs[2])*cos(θs[2]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  grad[3], = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*sin(θs[2])*cos(θs[1])*cos(θs[2]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);

  return grad;
end

function gradα(xs::Vector)
  global RELTOL_INT;
  global MAXEVALS_INT;
  
  grad = zeros(3);
  grad[1], = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*sin(θs[2])*cos(θs[1]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  grad[2], = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*cos(θs[2])*sin(θs[2])*cos(θs[1]), _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);
  grad[3], = pcubature(θs -> ρ(θs, xs)*sin(θs[2])*(sin(θs[2])*cos(θs[1]))^2, _LOWER_BOUNDS_, _UPPER_BOUNDS_; reltol=RELTOL_INT, maxevals=MAXEVALS_INT);

  return grad;
end

function hessian(xs::Vector)
  h = zeros(3, 3);

  h[:, 1] = gradc(xs);
  h[:, 2] = gradλ(xs);
  h[:, 3] = gradα(xs);

  return h;
end

function newton_raphson_iteration(xs::Vector, rs::Vector)
  h = hessian(xs);
  return xs - inv(h) * rs;
end

function newton_raphson(xs::Vector; reltol::Real=1e-10, abstol::Real=1e-8, max_iters::Int=100000)
  initial_Is = integrate_sys(xs);
  current_residuals = residuals(initial_Is);
  residual_norm = norm(current_residuals);
  init_residual_norm = residual_norm;
    
  if residual_norm < abstol
    return (xs, current_residuals, true, 0);
  end

  global print_every_percent_compl;
  iters_per_msg = ceil(Int, print_every_percent_compl * max_iters);

  for iter = 1:max_iters

    lin_alg_error_flag = false;
    try
      xs = newton_raphson_iteration(xs, current_residuals);
    catch
      lin_alg_error_flag = true;
    end

    # if Hessian is singular, or Newton iteration fails, perturb guess
    if lin_alg_error_flag
      @warn "hessian was singular; perturbing guess"
      xs += rand(3);
    end

    Is = integrate_sys(xs);
    current_residuals = residuals(Is);
    residual_norm = norm(current_residuals);

    relnorm = residual_norm / init_residual_norm;

    if residual_norm < abstol || relnorm < reltol
      return (xs, current_residuals, true, iter);
    end

    if iter % iters_per_msg == 0
      println("iter: $iter / $max_iters; Rnorm = $residual_norm");
    end

  end

  return (xs, current_residuals, false, max_iters);
end

function simulated_annealing(xs::Vector; anneal_factor::Real=1e-3,
                             init_T::Real=10.0, final_T::Real=1e-7,
                             relax_steps::Int=100,
                             init_delta::Real=0.5,
                             delta_adapt_factor::Real=1.2,
                             adapt_acc_thrshld::Real=0.7,
                             reltol::Real=1e-10, abstol::Real=1e-8, 
                             max_iters::Int=1000000)

  initial_Is = integrate_sys(xs);
  current_residuals = residuals(initial_Is);
  residual_norm = norm(current_residuals);
  init_residual_norm = residual_norm;
  Is = initial_Is;
    
  if residual_norm < abstol
    return (xs, current_residuals, true, 0);
  end

  acc = 0;
  delta = init_delta;
  T = init_T;

  for iter = 1:max_iters

    # check to see if relaxation should occur
    if iter % relax_steps == 0
      T *= (1.0 - anneal_factor);
      if T < final_T
        break;
      end

      # adapt step size to current acceptance ratio
      if acc / relax_steps > adapt_acc_thrshld
        delta *= delta_adapt_factor;
      elseif acc / relax_steps < (1 - adapt_acc_thrshld)
        delta /= delta_adapt_factor;
      end

      # reset acceptance count
      acc = 0;
    end

    # trial move
    choice = rand(1:3);
    dx = (rand(Float64)*2 - 1) * delta;
    new_xs = copy(xs);
    new_xs[choice] += dx;

    new_Is = integrate_sys(new_xs);
    new_residuals = residuals(new_Is);

    # metropolis acceptance criteria
    if rand(Float64) < exp(-(dot(new_residuals, new_residuals) - 
                             dot(current_residuals, current_residuals)) / T)

      # accept move
      acc += 1;
      copy!(xs, new_xs);
      current_residuals = new_residuals;
      Is = copy(new_Is);

      # check for convergence
      residual_norm = norm(current_residuals);
      relnorm = residual_norm / init_residual_norm;
      if residual_norm < abstol || relnorm < reltol
        return (xs, current_residuals, true, iter);
      end
      
    end
   
  end

  return (xs, current_residuals, false, iter);

end

_GRADS_ = [gradc, gradλ, gradα];

NLopt_raint_functions = map(i -> begin
  (xs::Vector, grad::Vector) -> begin
    if length(grad) > 0
      grad = gradc(xs);
    end

    return (pcubature(integrand_functions(xs)[i], _LOWER_BOUNDS_, _UPPER_BOUNDS_; 
                      reltol=RELTOL_INT, maxevals=MAXEVALS_INT)[1] - _RHS_[i]);
  end
end, 1:3);

function NLopt_optimization_function!(xs::Vector, grad::Vector)

  Is = integrate_sys(xs);
  current_residuals = residuals(Is);

  if length(grad) > 0
    h = hessian(xs);
    grad = vec(transpose(current_residuals) * h);
  end

  return 0.5 * dot(current_residuals, current_residuals);
end

 _CONV_SMBLS_ = [:SUCCESS; :STOPVAL_REACHED; :FTOL_REACHED; :XTOL_REACHED];

function run_solver(algo::Symbol; 
                    x0s::Vector = DEFAULT_INIT_GUESS, reltol::Real=1e-10,
                    abstol::Real=1e-8, max_iters::Int=10000)
  solver = init_solver(algo, x0s, reltol, abstol, max_iters);
  return solve(solver);
end

function run_solver(algos::Vector{Symbol}, maxs_iters::Vector{Int}; 
                    x0s::Vector = DEFAULT_INIT_GUESS, reltol::Real=1e-10,
                    abstol::Real=1e-8)

  if length(algos) != length(maxs_iters)
    error("Length of algorithm and maximum iteration vectors must be equal");
  end

  total_iterations = 0;
  xs, rs, converged = x0s, [NaN; NaN; NaN], false;

  for (algo, max_iters) in zip(algos, maxs_iters)
    (xs, rs, converged, iterations) = run_solver(algo; x0s=xs, 
                                                        reltol=reltol, 
                                                        abstol=abstol,
                                                        max_iters=max_iters);
    total_iterations += iterations;
  end

  return (xs, rs, converged, total_iterations);
end

function run_solver(algo_steps::Vector{Vector{Symbol}}, 
                    max_iters_by_step::Vector{Vector{Int}};
                    x0s::Vector = DEFAULT_INIT_GUESS, reltol::Real=1e-10,
                    abstol::Real=1e-6)
  
  if length(algo_steps) != length(max_iters_by_step)
    error("Length of algorithm and maximum iteration vectors must be equal");
  end

  total_iterations = 0;
  xs, rs, converged = x0s, [NaN; NaN; NaN], false;

  for (algo_step, max_iters_step) in zip(algo_steps, max_iters_by_step)
    (xs, rs, converged, iterations) = run_solver(algo_step, max_iters_step; 
                                                        x0s=x0s, 
                                                        reltol=reltol, 
                                                        abstol=abstol);
    total_iterations += iterations;
    if converged
      return (xs, rs, converged, total_iterations);
    end
  end
  
  return (xs, rs, converged, total_iterations);
end

 _NLOPT_SOLVERS_ = [:LN_PRAXIS; :LN_SBPLX; :GN_ISRES; :LD_SLSQP;
                         :LD_LBFGS; :AUGLAG];
 _NLOPT_REQ_BOUNDS_ = [:GN_ISRES];
 _NLOPT_REQ_RAINTS_ = [:GN_ISRES; :LD_SLSQP; :AUGLAG];
 _NLOPT_DEFAULT_LOWER_BOUNDS_ = [-N; -100.0; -100.0];
 _NLOPT_DEFAULT_UPPER_BOUNDS_ = [100.0*N; 100.0; 100.0];

function init_solver(algo::Symbol, x0s::Vector, reltol::Real, abstol::Real,
                     max_iters::Int)
  
  global _NUMEVALS_;
  global _ITERS_PER_PRINT_;
  global print_every_percent_compl;
  
  _NUMEVALS_ = 0;
  _ITERS_PER_PRINT_ = ceil(Int, print_every_percent_compl * max_iters);

  if algo in _NLOPT_SOLVERS_
    opt = Opt(algo, 3);
    min_objective!(opt, NLopt_optimization_function!);
    if algo in _NLOPT_REQ_BOUNDS_
      lower_bounds!(_NLOPT_DEFAULT_LOWER_BOUNDS_);
      upper_bounds!(_NLOPT_DEFAULT_UPPER_BOUNDS_);
    end
    if algo in _NLOPT_REQ_RAINTS_
      for i=1:3
        equality_raint!(opt, NLopt_raint_functions[i]); 
      end
    end
    if algo == :AUGLAG
      local_opt = init_solver(:LD_LBFGS, x0s, reltol, abstol, max_iters);
      local_optimizer!(opt, local_opt[1]);
    end
    ftol_rel!(opt, reltol);
    ftol_abs!(opt, abstol);
    maxeval!(opt, max_iters);
    return (opt, x0s);
  elseif algo == :NEWTON
    return ((xs) -> newton_raphson(xs; reltol=reltol, abstol=abstol, max_iters=max_iters),
            x0s);
  elseif algo == :ANNEAL
    return ((xs) -> simulated_annealing(xs; anneal_factor=ANNEAL_FACTOR,
                             init_T=ANNEAL_INIT_T, final_T=ANNEAL_FINAL_T,
                             relax_steps=ANNEAL_RELAX_STEPS,
                             init_delta=ANNEAL_INIT_DELTA,
                             delta_adapt_factor=ANNEAL_ADAPT_FACTOR,
                             adapt_acc_thrshld=ANNEAL_ACC_THRSHLD,
                             reltol=reltol, abstol=abstol,
                             max_iters=max_iters),
            x0s);
  else
    error("Do not understand algorithm symbol: ", algo);
  end
end

function solve(opt_pair::Tuple{NLopt.Opt, Vector})
  (minf, minx, ret) = optimize(opt_pair[1], opt_pair[2]);
  if ret == :XTOL_REACHED; @warn "Convergence given by x-tolerance"; end
  return (minx, sqrt(2 * minf), ret in _CONV_SMBLS_, _NUMEVALS_);
end

solve(newton_pair::Tuple{Function, Vector}) = newton_pair[1](newton_pair[2]);
