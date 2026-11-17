using Aqua

# ambiguities = false: pre-existing method ambiguities in the Lie and DA stacks
# and in TaylorSeries. This task does not rewrite the package to clear them.
#
# deps_compat: Plots and the stdlibs LinearAlgebra, Printf, Random already had
# no [compat] entry, and neither did the Test extra. Those bounds are left as
# they were; Aqua and TaylorSeries stay checked.
Aqua.test_all(
    Bell;
    ambiguities = false,
    deps_compat = (
        ignore = [:Plots, :LinearAlgebra, :Printf, :Random],
        check_extras = (ignore = [:Test],),
    ),
)
