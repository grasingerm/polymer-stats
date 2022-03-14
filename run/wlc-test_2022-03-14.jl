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
  "E0-$(fmt(case[:E0]))_K1-$(fmt(case[:K1]))_K2-$(fmt(case[:K2]))_kT-$(fmt(case[:kT]))_Fz-$(fmt(case[:Fz]))_Fx-$(fmt(case[:Fx]))_n-$(fmt(case[:n]))_b-$(fmt(case[:b]))";
end

cases = Any[];
E0s = [0.0];
κs = [1.0];
Ks = zip([1.0], [0.0]);
Fs = vcat(0.0, 0.05:0.05:1.0, 1.5:0.5:5.0, 7.5, 10.0, 25.0, 50.0, 100.0);
Fhats = [[0.0; 1.0]];
ns = Int[100];
bs = [1.0];
kT = 1.0;
for b in bs, n in ns, κ in κs, Kvec in Ks, E0 in E0s, Fmag in Fs, Fhat in Fhats
  Fx, Fz = Fmag*Fhat;
  K1, K2 = Kvec;
  push!(cases, Dict(:E0 => E0, :K1 => K1, :K2 => K2,
                    :kT => kT, :Fz => Fz, :kappa => κ,
                    :Fx => Fx, :n => n, :b => b, :run => run));
end

@info "total number of cases to run: $(length(cases))";

mkpath(workdir);

pmap(case -> begin;

  outfile = joinpath(workdir, "$(prefix(case)).out");
  if !isfile(outfile)
    println("Running case: $case.");
    command = `julia -O 3 mcmc_clustering_eap_chain.jl --chain-type dielectric --energy-type noninteracting -b $(case[:b]) --bend-mod $(case[:kappa]) --E0 $(case[:E0]) --K1 $(case[:K1]) --K2 $(case[:K2]) --kT $(case[:kT]) --Fz $(case[:Fz]) --Fx $(case[:Fx]) -n $(case[:n]) --num-steps 2000000 -v 2 --prefix $(joinpath(workdir, prefix(case)))`;
    output = read(command, String);
    write(outfile, output); 
  else
    println("Case: $case has already been run.");
  end

end, cases);
