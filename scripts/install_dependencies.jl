import Pkg;

@show depends = ([
                  "ArgParse";
                  "Distributions";
                  "Logging";
                  "DecFP";
                  "Quadmath";
                  "NLsolve";
                  "Plots";
                  "Optim";
                  "NLopt";
                  "Cubature";
                  "Glob"
                 ]);

foreach(Pkg.add, depends);

Pkg.update();
