using Glob;
using DelimitedFiles;
using Quadmath;
using DecFP;

if length(ARGS) < 4
    println("usage: julia aggregate.jl <outfile> <indir> <pattern> <dielectric|polar> [<3D|2D>] [<kappaflag>] [<runflag>]");
  exit(1);
end

kappaflag = if length(ARGS) >= 6
    parse(Bool, ARGS[6])
else
    false
end

runflag = if length(ARGS) >= 7
    parse(Bool, ARGS[7])
else
    false
end

dims = if length(ARGS) >= 5
         if ARGS[5] == "2D"
           2
         elseif ARGS[5] == "3D"
           3
         else
           println("usage: julia aggregate.jl <outfile> <indir> <pattern> <dielectric|polar> [<3D|2D>]");
           exit(1);
         end
       else
         3
       end

outfile = open(ARGS[1], "w");
indir = ARGS[2];
pattern = Glob.GlobMatch(ARGS[3]);

input_headers = if ARGS[4] == "dielectric"
    ["E0" "K1" "K2" "kT" "Fz" "Fx" "n" "b"]
elseif ARGS[4] == "polar"
    ["E0" "mu" "kT" "Fz" "Fx" "n" "b"]
else
  println("I don't understand the second to last argument");
  exit(1);
end

if kappaflag
    input_headers = hcat(input_headers, "kappa")
end

output_headers = if dims == 3
    ["r1" "r2" "r3" "lambda1" "lambda2" "lambda3" "r1sq" "r2sq" "r3sq" "rsquared" "p1" "p2" "p3" "p1sq" "p2sq" "p3sq" "psquared" "U" "Usquared" "Ealign" "psi" "AR"]
else
    ["r1" "r2" "lambda1" "lambda2" "r1sq" "r2sq" "rsquared" "p1" "p2" "p1sq" "p2sq" "psquared" "U" "Usquared" "Ealign" "psi" "AR"]
end

writedlm(outfile, hcat(input_headers, output_headers), ',');

for infile in readdir(pattern, indir)
  println("    processing $infile ... ")
  input_fields = split(split(basename(infile), ".")[1], "_")
  if runflag
      pop!(input_fields)
  end
  writedlm(outfile, 
           hcat(
                transpose(map(pair -> eval(Meta.parse(split(pair, "-"; limit=2)[2]))*1e-3,
                              input_fields)),
                map(line -> transpose(eval(Meta.parse(split(line, "=")[2]))),
                    readlines(infile))...
               ),
           ',');
end

close(outfile);
