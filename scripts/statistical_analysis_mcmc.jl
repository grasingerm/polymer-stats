using Glob;
using DelimitedFiles;
using Quadmath;
using DecFP;

if length(ARGS) != 3
  println("usage: julia statistical_analysis_mcmc.jl <outfile> <indir> <pattern>");
  exit(1);
end

outfile = open(ARGS[1], "w");
indir = ARGS[2];
pattern = Glob.GlobMatch(ARGS[3]);

data_series = Dict();
for infile in readdir(pattern, indir)
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

close(outfile);
