using Glob;

if length(ARGS) < 4
    println("usage: julia aggregate_by.jl <outdir> <indir> <param> <dielectric|polar> [<3D|2D>] [runflag]");
  exit(1);
end

dims = (length(ARGS) == 5) ? ARGS[5] : "3D";
runflag = (length(ARGS) == 6) ? parse(Bool, ARGS[6]) : false;

outdir = ARGS[1];
indir = ARGS[2];
param = ARGS[3];
chain_type = ARGS[4];

mkpath(outdir);

datafiles = readdir(glob"*.out", indir);
params_ran = [];
for datafile in datafiles
  datafile = basename(datafile);
  fileparams = split(string(split(basename(datafile), ".")[1:end-1]...), "_");
  if runflag; pop!(fileparams); end
  @show filtered_params = filter(x -> !startswith(x, param), fileparams);
  if filtered_params in params_ran
    println("-- param set already run -- \n");;
    continue;
  else
    push!(params_ran, filtered_params);
  end
  @show datafile, "$param-"
  @show strspan = findnext("$param-", datafile, 1);
  if strspan == nothing; continue; end
  value_start_index = strspan[end]+1;
  @show strspan2 = findnext("_", datafile, value_start_index);
  value_end_index = (strspan2 != nothing) ? strspan2[1] : length(datafile)-3;
  pattern = string(datafile[1:value_start_index-1], "*", datafile[value_end_index:end])
  @show pattern = (runflag) ? pattern[1:end-5] * "*" * pattern[end-3:end] : pattern
  @show value = datafile[value_start_index:value_end_index-1];
  @show outfile = joinpath(outdir, join(filtered_params, "_")*".csv");
  @show command = `julia scripts/aggregate_mcmc.jl $outfile $indir "$pattern" $chain_type $dims`;
  @show process_output = readlines(command);
end

display(params_ran);
println();
@show length(params_ran);
