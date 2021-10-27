using Glob;
using DelimitedFiles;
using Quadmath;
using DecFP;

if length(ARGS) != 4
  println("usage: julia aggregate.jl <outfile> <indir> <pattern> <dielectric|polar>");
  exit(1);
end

outfile = open(ARGS[1], "w");
indir = ARGS[2];
pattern = Glob.GlobMatch(ARGS[3]);

if ARGS[4] == "dielectric"
  writedlm(outfile, ["E0" "K1" "K2" "kT" "Fz" "Fx" "n" "b" "r1" "r2" "r3" "lambda1" "lambda2" "lambda3" "r1sq" "r2sq" "r3sq" "rsquared" "p1" "p2" "p3" "p1sq" "p2sq" "p3sq" "psquared" "U" "Usquared" "AR"], ',');
elseif ARGS[4] == "polar"
  writedlm(outfile, ["E0" "mu" "kT" "Fz" "Fx" "n" "b" "r1" "r2" "r3" "lambda1" "lambda2" "lambda3" "r1sq" "r2sq" "r3sq" "rsquared" "p1" "p2" "p3" "p1sq" "p2sq" "p3sq" "psquared" "U" "Usquared" "AR"], ',');
else
  println("I don't understand the last argument");
  exit(1);
end

for infile in readdir(pattern, indir)
  writedlm(outfile, 
           hcat(
                transpose(map(pair -> eval(Meta.parse(split(pair, "-"; limit=2)[2]))*1e-3,
                              split(split(basename(infile), ".")[1], "_"))),
                map(line -> transpose(eval(Meta.parse(split(line, "=")[2]))),
                    readlines(infile))...
               ),
           ',');
end

close(outfile);
