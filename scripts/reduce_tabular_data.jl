using Glob;
using DelimitedFiles;

USAGE_MSG = "usage: julia reduce_tabular_data.jl <indir> <outdir> <dielectric|polar>"

if length(ARGS) < 3
  println(USAGE_MSG);
  exit(1);
end

indir = ARGS[1];
outdir = ARGS[2];

nparams = if ARGS[3] == "dielectric"
  length(["E0" "K1" "K2" "kT" "Fz" "Fx" "n" "b" "kappa"]);
elseif ARGS[3] == "polar"
  length(["E0" "mu" "kT" "Fz" "Fx" "n" "b" "kappa"]);
else
  println("I don't understand the second to last argument");
  exit(1);
end

for infile in readdir("*.csv", indir)
  println("    processing $infile ... ")
  raw_data, headers = readdlm(infile, ','; header=true)
  pooled_data = Dict()
  for rowidx in 1:size(raw_data, 1)
    k = raw_data[rowidx, 1:nparams]
    if haskey(pooled_data, k)
        pooled_data[k] += [raw_data[nparams+1:end], 1]
    else
        pooled_data[k] = [raw_data[nparams+1:end], 1]
    end
  end
  reduced_data = [headers]
  for (k, v) in pooled_data
    push!(reduced_data, hcat(k, v[1] / v[2]))
  end
  reduced_data = vcat(reduced_data...)
  writedlm(joinpath(outdir, infile), reduced_data, ',')
end
