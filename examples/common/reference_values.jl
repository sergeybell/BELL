# reference_values.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Numbers copied from OptiM / COSY output files. Used only in comparison
# tables in examples 10 and 11. No model code path reads these constants.

const OPTIM_FILE = joinpath(@__DIR__, "..", "lattices", "magnetic_2.opt")
const STRUCTURE_FILE = joinpath(@__DIR__, "..", "lattices", "magnetic_2.bell")

# OptiM "Initial:" block of magnetic_2.opt (start of the lattice), listed
# here only for documentation; the parser reads the file at run time:
#   BetaX = 1936.75 cm, BetaY = 405.447 cm, AlfaX = -0.000251909,
#   AlfaY = -0.000188151, DispersX = 0, Dsp_PrimeX = -0.00050616

const COSY_REF = (
    Qx_frac = (0.4462838512100592, raw"COSY\sct_study\dat\magnetic_2\wp_natural.dat (SGx=SGy=0)"),
    Qy_frac = (0.4766824430028650, raw"COSY\sct_study\dat\magnetic_2\wp_natural.dat (SGx=SGy=0)"),
    alpha0  = (0.2303626034898868, raw"COSY\src\dat\magnetic_2\integrals_Icirc.dat / sct_study audit.json alpha0_mean"),
    chrom_x = (-1.840800730134095, raw"COSY\src\dat\magnetic_2\chrom_check_I0_FR0.dat chrom_x_DA (all SG=0)"),
    chrom_y = (-1.847678833791661, raw"COSY\src\dat\magnetic_2\chrom_check_I0_FR0.dat chrom_y_DA (all SG=0)"),
    betx0   = (19.3675, raw"COSY\src\dat\magnetic_2\BETAX, first row (element 0)"),
    bety0   = (4.05447, raw"COSY\src\dat\magnetic_2\BETAY, first row (element 0)"),
)
