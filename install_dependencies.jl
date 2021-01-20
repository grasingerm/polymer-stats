import Pkg;

@show depends = ([
                  "ArgParse";
                  "Distributions";
                  "Logging";
                  "ProfileView";
                  "DecFP";
                  "Quadmath";
                  "NLsolve";
                  "ForwardDiff";
                  "Plots";
                  "Optim";
                  "NLopt";
                  "Cubature";
                 ]);

foreach(Pkg.add, depends);

Pkg.update();
