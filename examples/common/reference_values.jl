# reference_values.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Numbers copied from OptiM / COSY output files for magnetic_2 (= magnetic_2_md:
# full ring, deuteron mass 1875.612945 MeV). Used only in comparison
# tables in examples 10 and 11. No model code path reads these constants.

const OPTIM_FILE = normpath(joinpath(@__DIR__, "..", "lattices", "magnetic_2.opt"))
const STRUCTURE_FILE = joinpath(@__DIR__, "..", "lattices", "magnetic_2.bell")

# OptiM "Initial:" block of magnetic_2.opt (start of the lattice), listed
# here only for documentation; the parser reads the file at run time:
#   BetaX = 1936.75 cm, BetaY = 405.447 cm, AlfaX = -0.000251909,
#   AlfaY = -0.000188151, DispersX = 0, Dsp_PrimeX = -0.00050616

const COSY_REF = (
    Qx_frac = (0.4617201483757378, raw"COSY OV 3 2 1 + TP, FR 0 (cosy_ov321.txt, magnetic_2_md)"),
    Qy_frac = (0.4613373212771295, raw"COSY OV 3 2 1 + TP, FR 0 (cosy_ov321.txt, magnetic_2_md)"),
    alpha0  = (0.23047028248812707, raw"COSY OV 3 3 0, FR 0 (ov330_map_analysis.json alpha_c)"),
    chrom_x = (-1.461535714233974, raw"COSY OV 3 2 1 + TP, FR 0 (cosy_ov321.txt xi_x_p)"),
    chrom_y = (-2.225484930000162, raw"COSY OV 3 2 1 + TP, FR 0 (cosy_ov321.txt xi_y_p)"),
    betx0   = (19.3675, raw"OptiM Initial BetaX of magnetic_2.opt"),
    bety0   = (4.05447, raw"OptiM Initial BetaY of magnetic_2.opt"),
)
