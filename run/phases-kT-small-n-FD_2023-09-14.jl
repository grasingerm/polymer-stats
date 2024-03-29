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
  "E0-$(fmt(case[:E0]))_K1-$(fmt(case[:K1]))_K2-$(fmt(case[:K2]))_kT-$(fmt(case[:kT]))_Fz-$(fmt(case[:Fz]))_Fx-$(fmt(case[:Fx]))_n-$(fmt(case[:n]))_b-$(fmt(case[:b]))_kappa-$(fmt(case[:kappa]))";
end

cases = Any[];
κs = [0.0];
Ks = zip([0.0], [1.0]);
#E0s = vcat(0.0, 0.005:0.005:1.0);
E0s = [1.0];
ns = Int[5; 10; 15; 20];
bs = [1.0; 0.1; 0.01];
kTs = vcat(0.005:0.005:1.1, 1.5:0.25:10.0, 20.0, 25.0);
#Fs = zip([0.0; 1.0; 0.0], [0.0; 0.0; 1.0]);
Fs = zip([0.0], [0.0]);
for b in bs, κ in κs, n in ns, Kvec in Ks, E0 in E0s, kT in kTs, F in Fs
  K1, K2 = Kvec;
  push!(cases, Dict(:E0 => E0, :K1 => K1, :K2 => K2,
                    :kT => kT, :Fz => F[2], :kappa => κ,
                    :Fx => F[1], :n => n, :b => b));
end

@info "total number of cases to run: $(length(cases))";

mkpath(workdir);

pmap(case -> begin;

  outfile = joinpath(workdir, "$(prefix(case)).out");
  if !isfile(outfile)
    println("Running case: $case.");
    command = `julia -O 3 -t 1 mcmc_clustering_eap_chain.jl --chain-type dielectric --energy-type interacting --x0 "[0.0; π/2]" -b $(case[:b]) --bend-mod $(case[:kappa]) --E0 $(case[:E0]) --K1 $(case[:K1]) --K2 $(case[:K2]) --kT $(case[:kT]) --Fz $(case[:Fz]) --Fx $(case[:Fx]) -n $(case[:n]) --num-steps 10000000 --burn-in 100000 -v 2 --prefix $(joinpath(workdir, prefix(case)))`;
    output = read(command, String);
    write(outfile, output); 
  else
    println("Case: $case has already been run.");
  end

end, cases);
