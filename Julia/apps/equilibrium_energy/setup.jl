using Pkg
Pkg.develop(path = raw"C:\GIT_REPOS\BELL\Julia\ring")
Pkg.instantiate()
using LieRing
println("OK ", pkgversion(LieRing))
