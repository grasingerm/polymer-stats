using Glob;
using DelimitedFiles;

if length(ARGS) != 3
  println("usage: julia aggregate.jl <outfile> <indir> <pattern>");
  exit(1);
end

outfile = open(ARGS[1], "w");
indir = ARGS[2];
pattern = Glob.GlobMatch(ARGS[3]);

writedlm(outfile, ["E0" "K1" "K2" "kT" "Fz" "Fx" "n" "b" "r1" "r2" "r3" "lambda1" "lambda2" "lambda3" "rsquared" "p1" "p2" "p3" "U" "AR"], ',');

for infile in readdir(pattern, indir)
  writedlm(outfile, 
           hcat(
                transpose(map(pair -> eval(Meta.parse(split(pair, "-")[2]))*1e-3,
                          split(split(infile, ".")[1], "_"))),
                map(line -> transpose(eval(Meta.parse(split(line, "=")[2]))),
                    readlines(infile))...
               ),
           ',');
end

close(outfile);