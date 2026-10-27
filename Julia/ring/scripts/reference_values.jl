# =============================================================================
# reference_values.jl - external reference numbers (OUTPUT files only)
# =============================================================================
#
# WHAT
#   Numbers copied verbatim from OptiM / COSY OUTPUT/data files, with the
#   file they come from.  They are used ONLY for the final comparison tables
#   in scripts/01 and scripts/02; no code path of the model uses them.
#
#   NOTE: structures/magnetic_2/magnetic_2.lis and magnetic_2_maps.lis
#   (requested as reference) contain only the echo of the COSY input deck and
#   "BIN FILE WRITTEN" lines - no numeric optics results.  The numbers below
#   are therefore taken from the COSY result files listed next to each value.
#
# WHY
#   Keeps the comparison honest and traceable.
#
# WHAT YOU SEE
#   Nothing printed by this file.
# =============================================================================
const OPTIM_FILE = raw"C:\GIT_REPOS\SCT\OptiM\magnetic\magnetic_2.opt"

# OptiM "Initial:" block of magnetic_2.opt (start of the lattice).  Parsed at
# run time from the file; listed here only for documentation:
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
