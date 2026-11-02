# OptiM (OptiMX x86-64) — CLI

`optimx64.exe -s lattice.opt out.txt` (MAD8-like table, m) and `-f` (cm). Extra args: filter MatchCase HideFirstLetter Comment step CloseLattice (defaults `* Y Y # 0 Y`).
Other modes: -c coupled, -i space charge, -t trajectory, -o orbit, -p tracking, -w/-r snap files.
CLI has no chromaticity / alpha_c / integrals: these are only in GUI View | Integrals (help/htmldoc/view_menu.html); there is no script/macro mode in the htmldoc.
Dispersion in the CLI table starts from the file `Initial:` block (D=0), so it is not the periodic solution (END: D = 2.19 cm).
`gui_reference_manual.txt` = numbers typed by the user from the GUI (manual reference only).
