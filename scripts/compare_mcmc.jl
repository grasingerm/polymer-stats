using Glob;
using DelimitedFiles;
using GLM;
using DataFrames;
using Logging;

if length(ARGS) != 2
  println("usage: julia statistical_analysis_mcmc.jl <indir> <subpattern>");
  exit(1);
end

global_logger(ConsoleLogger(stderr, Logging.Info));

indir = ARGS[1];
subpattern = ARGS[2];
pattern = Glob.GlobMatch(subpattern * ".out");

chunks = map(piece -> split(piece, "-"), split(subpattern, "_"));
param_map = Dict();
for chunk in chunks
  if chunk[1] == "run"; continue; end
  param_map[chunk[1]] = Meta.parse(chunk[end]) / 1000.0;
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

Ustar = (param_map["E0"])^2 * abs(param_map["K2"] - param_map["K1"]);
rstar = param_map["b"];
μstar = (param_map["E0"])*(2*param_map["K2"] + param_map["K1"]);
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
end
