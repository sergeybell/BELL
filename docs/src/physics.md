# Physics Reference

The physics of Bell — coordinates, the s-Hamiltonian, every element including
the RF cavity, fringe fields, linear optics, the Lie product, map readings,
the normal form, RF and acceleration, sextupoles, Lie maps and differential
algebra, the computation paths, and the modelling assumptions — is described
in the **Beam Physics Manual**, with the same chapters in English and Russian:

```@raw html
<p>English: <code>docs/manual/en/</code>,
<code>julia --project=docs docs/make_pdf.jl en</code>
writes <code>docs/pdf/Bell_Beam_Physics_Manual_en.pdf</code>.
Russian: <code>docs/manual/ru/</code>,
<code>julia --project=docs docs/make_pdf.jl ru</code>
writes <code>docs/pdf/Bell_Beam_Physics_Manual_ru.pdf</code>.
Shared figures: <code>docs/manual/assets/</code>.</p>
```

This site does not repeat that text. How the formulas are implemented is in the
[Developer's Guide](dev/architecture.md).
