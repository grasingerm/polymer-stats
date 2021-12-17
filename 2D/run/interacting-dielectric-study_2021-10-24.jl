println("Hello world!");

@everywhere using Printf;

@show workdir = if length(ARGS) > 0
  ARGS[1];
else
  ".";
end

@everywhere fmt(x) = @sprintf("%07d", round(Int, 1e3*x));
@everywhere fmt_int(x) = @sprintf("%03d", x);

@everywhere function prefix(case)
  "E0-$(fmt(case[:E0]))_K1-$(fmt(case[:K1]))_K2-$(fmt(case[:K2]))_kT-$(fmt(case[:kT]))_Fz-$(fmt(case[:Fz]))_Fx-$(fmt(case[:Fx]))_n-$(fmt(case[:n]))_b-$(fmt(case[:b]))_run-$(fmt_int(case[:run]))";
end

cases = Any[];
#E0s = [1.0; 5.0];
#K1s = [1e-1; 1.0]; 
#kTs = [1.0; 0.1];
#Fzs = [0.0; 0.5; 1.0];
#Fxs = [0.0; 0.5; 1.0];
#ns = Int[100];
#bs = [1.0];
runs = 1:25;
for run in runs
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 0,
                    :Fx => 0, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 5, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 0,
                    :Fx => 0, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 5, :K1 => 0.1, :K2 => 0.0,
                    :kT => 1, :Fz => 0,
                    :Fx => 0, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 0.1, :K2 => 0.0,
                    :kT => 1, :Fz => 0,
                    :Fx => 0, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 1,
                    :Fx => 0, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 2,
                    :Fx => 0, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 0,
                    :Fx => 1, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 0.1, :Fz => 0,
                    :Fx => 1, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 0.1, :Fz => 0,
                    :Fx => 0, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 0.1, :Fz => 1,
                    :Fx => 0, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 0.1, :K2 => 0.0,
                    :kT => 0.1, :Fz => 1,
                    :Fx => 0, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 0,
                    :Fx => 0, :n => 200, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 0,
                    :Fx => 0, :n => 100, :b => 0.5, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 0,
                    :Fx => 0, :n => 100, :b => 2.0, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 1,
                    :Fx => 0, :n => 100, :b => 0.5, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 0,
                    :Fx => 1, :n => 100, :b => 0.5, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 1,
                    :Fx => 1, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 3, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 3,
                    :Fx => 3, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 3,
                    :Fx => 0, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 5,
                    :Fx => 0, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 0,
                    :Fx => 3, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 1, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 0,
                    :Fx => 5, :n => 100, :b => 1, :run => run));
  push!(cases, Dict(:E0 => 5, :K1 => 1, :K2 => 0.0,
                    :kT => 1, :Fz => 5,
                    :Fx => 5, :n => 100, :b => 1, :run => run));
end

@info "total number of cases to run: $(length(cases))";

mkpath(workdir);

pmap(case -> begin;

  outfile = joinpath(workdir, "$(prefix(case)).out");
  if !isfile(outfile)
    println("Running case: $case.");
    command = `julia -O 3 mcmc_clustering_eap_chain.jl --chain-type dielectric --energy-type Ising -b $(case[:b]) --E0 $(case[:E0]) --K1 $(case[:K1]) --K2 $(case[:K2]) --kT $(case[:kT]) --Fz $(case[:Fz]) --Fx $(case[:Fx]) -n $(case[:n]) --num-steps 1000000 -v 2 --prefix $(joinpath(workdir, prefix(case)))`;
    output = read(command, String);
    write(outfile, output); 
  else
    println("Case: $case has already been run.");
  end

end, cases);
