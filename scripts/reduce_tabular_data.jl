using Glob;
using DelimitedFiles;

USAGE_MSG = "usage: julia reduce_tabular_data.jl <indir> <outdir> <dielectric|polar> <kappaflag>"

if length(ARGS) < 3
  println(USAGE_MSG);
  exit(1);
end

indir = ARGS[1];
outdir = ARGS[2];

nparams = if ARGS[3] == "dielectric"
  length(["E0" "K1" "K2" "kT" "Fz" "Fx" "n" "b"]);
elseif ARGS[3] == "polar"
  length(["E0" "mu" "kT" "Fz" "Fx" "n" "b"]);
else
  println("I don't understand the second to last argument");
  exit(1);
end

if length(ARGS) >= 4 && ARGS[4] == "true"
    nparams += 1
end

@show nparams

for infile in readdir(glob"*.csv", indir)
  println("    processing $infile ... ")
  raw_data, headers = readdlm(infile, ','; header=true)
  pooled_data = Dict()
  @show headers
  for rowidx in 1:size(raw_data, 1)
    @show k = raw_data[rowidx, 1:nparams]
    @show [raw_data[rowidx, nparams+1:end], 1]
    try
        if haskey(pooled_data, k)
            pooled_data[k] += [raw_data[rowidx, nparams+1:end], 1]
        else
            pooled_data[k] = [raw_data[rowidx, nparams+1:end], 1]
        end
    catch e
        @show e
    end
  end
  reduced_data = Any[headers]
  for k in sort(collect(keys(pooled_data)))
      v = pooled_data[k]
      @show transpose(k), transpose(v[1]), v[2]
      push!(reduced_data, hcat(transpose(k), transpose(v[1] / v[2])))
  end
  reduced_data = vcat(reduced_data...)
  writedlm(joinpath(outdir, basename(infile)), reduced_data, ',')
end
