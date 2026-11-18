@testset "insert_cavity keeps the ring harmonic" begin
    ref = reference_from_kinetic(deuteron, 10.0)
    lat = Lattice("r", AbstractElement[Drift("D", 1.0), Marker("M")], ref, 3)
    cav = RFCavity("RF", 0.1, 7, 0.0)
    lat2, idx = insert_cavity(lat, cav; in_element = "D")
    @test lat2.harmonic == 3
    @test lat2.elements[idx] isa RFCavity
    @test lat2.elements[idx].harmonic == 7
    lat3, _ = insert_cavity(lat, cav; in_element = "D", harmonic = 9)
    @test lat3.harmonic == 9
    @test_throws ErrorException insert_cavity(lat, cav; in_element = "M")
    @test_throws ErrorException insert_cavity(lat, cav; in_element = "missing")
end

@testset "sextupole bounds must match the lattice" begin
    ref = reference_from_kinetic(deuteron, 10.0)
    lat = Lattice("r", AbstractElement[Drift("D", 1.0)], ref, 1)
    @test_throws ErrorException place_sextupole_families(lat, (bounds = Tuple{Float64,Float64,AbstractElement}[],))
    other = (bounds = [(0.0, 1.0, Drift("other", 1.0))],)
    @test_throws ErrorException place_sextupole_families(lat, other)
end

@testset "convert_optim close_ring" begin
    opt = tempname() * ".opt"
    write(opt, "OptiM\nEnergy[MeV]=270   Mass[MeV]=1875.612945\nbegin lattice. Number of periods=1\n" *
               join(fill("O1 BM", 16), " ") * "\nend lattice\nbegin list\nO1 L[cm]=10\n" *
               "BM L[cm]=532.94148 B[kG]=2.56106599 G[kG/cm]=0\nend list\n")
    for (cr, tol) in ((true, 1e-15), (false, Inf))
        bell = tempname() * ".bell"
        m = convert_optim(opt, bell; close_ring = cr)
        lat = load_structure(bell)
        tot = sum(e -> e isa SBend ? e.angle : 0.0, lat.elements)
        cr ? (@test abs(tot - 2pi) < 1e-13 && abs(lat.elements[2].angle - pi / 8) < 1e-15) :
             (@test 1e-9 < abs(tot - 2pi) < 1e-6)
        @test any(w -> occursin("close_ring", w), m.warnings) == cr
    end
end

@testset "structure round-trip magnetic_2" begin
    src = joinpath(pkgdir(Bell), "examples", "lattices", "magnetic_2.bell")
    lat = load_structure(src)
    @test abs(total_bend_angle(lat) - 2pi) < 1e-6
    mktempdir() do dir
        out = joinpath(dir, "magnetic_2.bell")
        write_structure(out, lat)
        lat2 = load_structure(out)
        @test path_length(lat) == path_length(lat2)
        @test abs(total_bend_angle(lat2) - 2pi) < 1e-6
        @test length(lat.elements) == length(lat2.elements)
        @test [typeof(e) for e in lat.elements] == [typeof(e) for e in lat2.elements]
        @test [e.name for e in lat.elements] == [e.name for e in lat2.elements]
        @test lat.ref.Ekin_MeV == lat2.ref.Ekin_MeV
    end
end
