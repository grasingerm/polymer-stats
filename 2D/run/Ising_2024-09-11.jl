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
    "E0-$(fmt(case[:E0]))_K1-$(fmt(case[:K1]))_K2-$(fmt(case[:K2]))_kT-$(fmt(case[:kT]))_Fz-$(fmt(case[:Fz]))_Fx-$(fmt(case[:Fx]))_n-$(fmt(case[:n]))_b-$(fmt(case[:b]))_run-$(case[:run])";
end

cases = Any[];
E0s = [1e-1;];
K1s = [0.0]; 
#K2s = [0.5; 1.0; 2.0]; 
K2s = 1e-2*[1.0; 4.0]; 
kTs = [1.0];
Fxs = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 7.5; 10.0; 12.5; 15.0; 20.0; 25.0; 30.0; 
       35.0; 40.0; 45.0; 50.0; 60.0; 70.0; 80.0; 90.0; 100.0; 125.0; 150.0; 
       175.0; 200.0; 250.0; 300.0; 400.0; 500.0];
Fzs = [0.0];
#ns = Int[100; 200; 500];
ns = Int[25];
#bs = [0.5; 1.0; 2.0];
bs = [1.0; 4.0];
runs = 1:10
for b in bs, n in ns, Fx in Fxs, Fz in Fzs, kT in kTs, E0 in E0s, K1 in K1s, K2 in K2s, run in runs
  push!(cases, Dict(:E0 => E0, :K1 => K1, :K2 => K2,
                    :kT => kT, :Fz => Fz,
                    :Fx => Fx, :n => n, :b => b, :run => run));
end

@info "total number of cases to run: $(length(cases))";

mkpath(workdir);

pmap(case -> begin;

  outfile = joinpath(workdir, "$(prefix(case)).out");
  if !isfile(outfile)
    println("Running case: $case.");
    command = `/home/grasinmj/julia-1.10.4/bin/julia -O 3 mcmc_clustering_eap_chain.jl --chain-type dielectric --energy-type Ising -b $(case[:b]) --E0 $(case[:E0]) --K1 $(case[:K1]) --K2 $(case[:K2]) --kT $(case[:kT]) --Fz $(case[:Fz]) --Fx $(case[:Fx]) -n $(case[:n]) --num-steps 10000000 --burn-in 200000 -v 2 --prefix $(joinpath(workdir, prefix(case)))`;
    output = read(command, String);
    write(outfile, output); 
  else
    println("Case: $case has already been run.");
  end

end, cases);
