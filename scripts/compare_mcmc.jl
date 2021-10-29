using Glob;
using DelimitedFiles;
using GLM;
using DataFrames;
using Logging;
using Printf;

if length(ARGS) != 2
  println("usage: julia statistical_analysis_mcmc.jl <indir> <subpattern>");
  exit(1);
end

global_logger(ConsoleLogger(stderr, Logging.Info));

fmt(x) = @sprintf("%+0.4f", x);

indir = ARGS[1];
subpattern = ARGS[2];
pattern = Glob.GlobMatch(subpattern * ".out");

chunks = map(piece -> split(piece, "-"), split(subpattern, "_"));
param_map = Dict();
for chunk in chunks
  try
    param_map[chunk[1]] = Meta.parse(chunk[end]) / 1000.0;
  catch e
    println("Could not process $(chunk)");
  end
end
@show param_map;

data_serieses = Dict();
for subdir in ["standard"; "clustering"; "umbrella"]
  data_series = Dict();
  for infile in readdir(pattern, joinpath(indir, subdir))
    for line in readlines(infile)
      var, val = map(strip, split(line, "="));
      val = eval(Meta.parse(val));
      if haskey(data_series, var)
        push!(data_series[var], val);
      else
        data_series[var] = [val];
      end
    end
  end
  data_serieses[subdir] = data_series;
end
@show map(keys, values(data_serieses));

colidxs = Dict();
for (idx, colname) in enumerate(["step","r1","r2","r3","r1sq","r2sq","r3sq","rsq","p1","p2","p3","p1sq","p2sq","p3sq","psq","U","Usq"])
  colidxs[colname] = idx;
end
Ustar = (param_map["E0"])^2 * abs(param_map["K2"] - param_map["K1"]);
rstar = param_map["b"];
pstar = (param_map["E0"])*(2*param_map["K2"] + param_map["K1"]);
for (name, data_series) in data_serieses
  @show name;
  avgs = Dict();
  vars = Dict();
  std_norms = Dict();
  for (k, v) in data_series
    avgs[k] = sum(v) / length(v);
    vars[k] = sum(map(x -> (x - avgs[k]).^2, v)) / length(v);
    std_norms[k] = map(sqrt, vars[k]) ./ avgs[k];
  end
  for k in keys(avgs)
    @info "$k: μ = $(avgs[k]); σ^2 = $(vars[k]); σ / μ = $(std_norms[k])";
  end
  rolling_data = map(datafile -> readdlm(datafile, ','; skipstart=1),
                     readdir(Glob.GlobMatch(subpattern * "_rolling.csv"), 
                             joinpath(indir, name))
                    );
  L1s = Dict();
  αs = [];
  ϵs = [];
  for (varavg, cname, star) in zip([avgs["<r>"][1], avgs["<r>"][3], avgs["<p>"][1], avgs["<p>"][3], avgs["<U>"]],
                             ["r1", "r3", "p1", "p3", "U"],
                             [rstar, rstar, pstar, pstar, Ustar])
    L1s[cname] = map(i -> begin;
                       sum(map(rd -> abs(varavg - rd[i, colidxs[cname]]), rolling_data)) / (star*length(rolling_data));
              end,
              1:size(rolling_data[1], 1));
    df = DataFrame(X = map(log, rolling_data[1][:, 1]),
                   Y = map(log, L1s[cname]));
    fit = lm(@formula(Y ~ X), df);
    α = coef(fit)[2];
    ϵ = L1s[cname][end];
    println("$cname: $(fmt(α)), $(fmt(ϵ))");
    push!(αs, α);
    push!(ϵs, ϵ);
  end
  println("$(fmt(minimum(αs))); $(fmt(maximum(αs)))");
  println("$(fmt(minimum(ϵs))); $(fmt(maximum(ϵs)))");
end
