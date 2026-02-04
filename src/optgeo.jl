const BOHR_TO_ANGSTROM = 0.5291772109
const ANGSTROM_TO_BOHR = 1.8897261246

function optgeo(mol::String, bset_ref, p_ref)

    # Move coordinates so center of mass is at origin
    mt, cm = center_of_mass(bset_ref)
    spin_mult, symbols, coords = string_to_xyz(mol, cm)
    mol_at_cm = xyz_data_to_string(spin_mult, symbols, coords)

    # Reference Single Point
    bset, _ = DoNOF.molecule(mol_at_cm, bset_ref.name, spherical = true)
    p = deepcopy(p_ref)
    E, C_ref, n_ref = DoNOF.energy(bset, p)

    # Optimize geometry (coords in bohr)
    coords = coords .* ANGSTROM_TO_BOHR
    res = optimize(
        c -> e_linearized(c, symbols, spin_mult, C_ref, n_ref, bset, p_ref),
        c -> g_linearized(c, symbols, spin_mult, C_ref, n_ref, bset, p_ref),
        coords,
        LBFGS(),
        Optim.Options(g_abstol = 10^-4),
        inplace = false,
    )
    
    # coords de Bohr a Angstrom
    coords = res.minimizer * BOHR_TO_ANGSTROM
    coords[abs.(coords) .< 10^-4] .= 0
    println(coords)

    mol_optimized = xyz_data_to_string(spin_mult, symbols, coords)

    return mol_optimized

end
