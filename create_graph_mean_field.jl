include(joinpath("inc", "algorithm_help_msg.jl"));

using ArgParse;

s = ArgParseSettings();
@add_arg_table s begin
  "--E0", "-z"
    help = "magnitude of electric field"
    arg_type = Float64
    default = 0.0
  "--dk", "-D"
    help = "kappa1 - kappa2"
    arg_type = Float64
    default = 1.0
  "--kT", "-k"
    help = "dimensionless temperature"
    arg_type = Float64
    default = 1.0
  "--theta", "-t"
    help = "angle of stretch direction (with respect to x-axis; e-field direction is pi/2)"
    arg_type = AbstractString
    default = "pi"
  "--num-monomers", "-N"
    help = "number of monomers"
    arg_type = Int
    default = 100
  "--algorithm", "-A"
    help = "algorithm (ISRES | PRAXIS | SBPLX | SLSQP | LBFGS | AUGLAG | Newton | Anneal). Run with --print-algo-help for more detailed information"
    default = "try Newton then PRAXIS,SBPLX,Newton"
  "--print-algo-help"
    help = "Get detailed information regarding algorithms themselves, and how to run to solvers in stages"
    action = :store_true
  "--reltol-int", "-r"
    help = "relative tolerance for numerical integration"
    arg_type = Float64
    default = 1e-6
  "--maxevals-int", "-e"
    help = "maximum evaluations for numerical integration"
    arg_type = Int
    default = 2500
  "--reltol-opt", "-R"
    help = "relative tolerance for optimization"
    arg_type = Float64
    default = 0.0
  "--abstol-opt", "-B"
    help = "absolute tolerance for optimization"
    arg_type = Float64
    default = 1e-8
  "--maxevals-opt", "-E"
    help = "maximum evaluations for numerical integration"
    arg_type = AbstractString
    default = "try 1000 then 10000,5000,50"
  "--num-data-points"
    help = "number of data points to collect in free energy vs. stretch graph"
    arg_type = Int
    default = 10
  "--outfile", "-O"
    help = "filename"
    default = "energy-vs-stretch.csv"
  "--verbose", "-v"
    help = "print diagnostic information"
    action = :store_true
end

pa = parse_args(s);
if pa["print-algo-help"]
  println(_ALGO_HELP_MSG_);
  exit();
end
const dk = pa["dk"];
const E0 = pa["E0"];
const theta = eval(parse(pa["theta"]));
const st = sin(theta);
const ct = cos(theta);
const kT = pa["kT"];
const N = pa["num-monomers"];
const reltol = pa["reltol-int"];
const maxevals = pa["maxevals-int"];
const numpoints = pa["num-data-points"];
const outfile = pa["outfile"];
const algo = pa["algorithm"];
const verbose = pa["verbose"];
const reltol_opt = pa["reltol-opt"];
const abstol_opt = pa["abstol-opt"];
const maxevals_opt = pa["maxevals-opt"];

x0 = 1.0 / numpoints / 10;
dx = 1.0 / (numpoints - 1);
xs = vcat([x0], collect(dx:dx:1.0));
xs[end]-=x0;

open(outfile, "w") do w
  for x in xs
    println("x = $x");
    if verbose
      command = `julia mean_field.jl -A $algo --E0 $E0 --dk $dk --kT $kT -N $N --reltol-int $reltol --maxevals-int $maxevals --reltol-opt $reltol_opt --abstol-opt $abstol_opt --maxevals-opt $maxevals_opt --gz $(x*st) --gx $(x*ct) --verbose`;
    else
      command = `julia mean_field.jl -A $algo --E0 $E0 --dk $dk --kT $kT -N $N --reltol-int $reltol --maxevals-int $maxevals --reltol-opt $reltol_opt --abstol-opt $abstol_opt --maxevals-opt $maxevals_opt --gz $(x*st) --gx $(x*ct)`;
    end

    process_output = readlines(command);
    free_energy = NaN;
    lambda = NaN;
    alpha = NaN;
    C = NaN;
    mu1 = NaN;
    mu2 = NaN;
    mu3 = NaN;
    l2 = NaN;

    if verbose; println("output:"); end
    for line in process_output
      if verbose
        println(">>  ", strip(line));
      end
      splits = map(strip, split(line, '='));
      if length(splits) > 1 
        if splits[1] == "A"
          free_energy = parse(Float64, splits[2]);
        elseif splits[1] == "lambda"
          lambda = parse(Float64, splits[2]);
        elseif splits[1] == "alpha"
          alpha = parse(Float64, splits[2]);
        elseif splits[1] == "C"
          C = parse(Float64, splits[2]);
        elseif splits[1] == "mu1"
          mu1 = parse(Float64, splits[2]);
        elseif splits[1] == "mu2"
          mu2 = parse(Float64, splits[2]);
        elseif splits[1] == "mu3"
          mu3 = parse(Float64, splits[2]);
        elseif splits[1] == "L2"
          l2 = parse(Float64, splits[2]);
        end
      end
    end
    if verbose; println(); end
    println("================ RESULTS ================");
    
    println("        x:  $x");
    println("        A:  $free_energy");
    println("   lambda:  $lambda");
    println("    alpha:  $alpha");
    println("        C:  $C");
    println("      mu1:  $mu1");
    println("      mu2:  $mu2");
    println("      mu3:  $mu3");
    println("       L2:  $l2");
    println("=========================================\n");
    write(w, "$x, $free_energy, $lambda, $alpha, $C, $mu1, $mu2, $mu3, $l2\n");
  end
end
