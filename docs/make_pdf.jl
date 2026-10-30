# =============================================================================
# docs/make_pdf.jl - PDF from docs/src Markdown (beam physics or programmer)
# Usage: julia docs/make_pdf.jl              -> Beam Physics PDF
#        julia docs/make_pdf.jl programmer   -> Programmer PDF
# =============================================================================
const DOC = @__DIR__
const MODE = length(ARGS) > 0 && ARGS[1] == "programmer" ? :programmer : :beam
const PAGES = MODE === :programmer ?
    ["programmer/index.md"] :
    ["architecture.md", "index.md", "canonical.md", "hamiltonians.md", "lie.md", "linear.md", "second_order.md",
     "chromaticity.md", "fringe.md", "rf.md", "acceleration.md", "longitudinal.md", "sextupoles.md", "equilibrium.md", "assumptions.md"]
const OUTNAME = MODE === :programmer ? "Bell_Programmer_Manual" : "Bell_Beam_Physics_Manual"
const TITLE = MODE === :programmer ?
    "Bell Programmer Manual\\\\\\large modules, examples, structure file" :
    "Bell Beam Physics Manual\\\\\\large canonical variables, Lie maps, RF, acceleration"

esc_text(s) = replace(s, "\\" => "\\textbackslash{}", "&" => "\\&", "%" => "\\%", "#" => "\\#",
                      "_" => "\\_", "\$" => "\\\$", "^" => "\\^{}", "~" => "\\~{}", "{" => "\\{", "}" => "\\}")
function inline(s)
    out = IOBuffer(); parts = split(s, "``")
    for (k, p) in enumerate(parts)
        if iseven(k); print(out, "\$", p, "\$")
        else
            q = replace(p, r"\[([^\]]*)\]\((@ref|[\w./]+\.md)\)" => s"\1")
            cs = split(q, "`")
            for (j, c) in enumerate(cs)
                if iseven(j); print(out, "\\texttt{", esc_text(c), "}")
                else
                    t = esc_text(c)
                    t = replace(t, r"\*\*([^*]+)\*\*" => s"\\textbf{\1}")
                    print(out, t)
                end
            end
        end
    end; String(take!(out))
end
function convert_md(md; srcfile=joinpath(DOC, "src"))
    L = split(replace(md, "\ufeff" => "", "\r" => ""), "\n")
    out = String[]; i = 1; list = ""
    closelist() = (list != "" && (push!(out, "\\end{$list}"); list = ""))
    while i <= length(L)
        l = L[i]
        if startswith(l, "```math")
            closelist(); j = i + 1; buf = String[]
            while !startswith(L[j], "```"); push!(buf, L[j]); j += 1; end
            push!(out, "\\[", join(buf, "\n"), "\\]"); i = j + 1; continue
        elseif startswith(l, "```")
            closelist(); j = i + 1; buf = String[]
            while !startswith(L[j], "```"); push!(buf, L[j]); j += 1; end
            push!(out, "\\begin{verbatim}", join(buf, "\n"), "\\end{verbatim}"); i = j + 1; continue
        elseif startswith(l, "## "); closelist(); push!(out, "\\subsection{" * inline(l[4:end]) * "}")
        elseif startswith(l, "# "); closelist(); push!(out, "\\section{" * inline(l[3:end]) * "}")
        elseif startswith(l, "|")
            closelist(); rows = String[]
            while i <= length(L) && startswith(L[i], "|"); push!(rows, L[i]); i += 1; end
            cells(r) = strip.(split(strip(r)[2:end-1], "|")); hdr = cells(rows[1]); n = length(hdr)
            w = round(0.95 / n, digits = 3)
            push!(out, "\\begin{center}\\small\\begin{tabular}{" * join(fill("p{$(w)\\linewidth}", n), "") * "}\\hline")
            push!(out, join(["\\textbf{" * inline(c) * "}" for c in hdr], " & ") * "\\\\\\hline")
            for r in rows[3:end]; push!(out, join(inline.(cells(r)), " & ") * "\\\\"); end
            push!(out, "\\hline\\end{tabular}\\end{center}"); continue
        elseif (m = match(r"^!\[(.*)\]\((.*)\)\s*$", l)) !== nothing
            closelist(); cap, path = m.captures
            if endswith(path, ".gif")
                push!(out, "\\begin{quote}\\emph{Animation (see HTML or GIF file): " * inline(cap) * "}\\end{quote}")
            else
                fp = replace(relpath(normpath(joinpath(dirname(srcfile), path)), joinpath(DOC, "pdf")), "\\" => "/")
                push!(out, "\\begin{figure}[htbp]\\centering\\includegraphics[width=\\linewidth]{$fp}\\caption{" * inline(cap) * "}\\end{figure}")
            end
        elseif (m = match(r"^(\d+)\.\s+(.*)$", l)) !== nothing
            list != "enumerate" && (closelist(); push!(out, "\\begin{enumerate}"); list = "enumerate")
            push!(out, "\\item " * inline(m.captures[2]))
        elseif (m = match(r"^\*\s+(.*)$", l)) !== nothing
            list != "itemize" && (closelist(); push!(out, "\\begin{itemize}"); list = "itemize")
            push!(out, "\\item " * inline(m.captures[1]))
        elseif strip(l) == ""; closelist(); push!(out, "")
        else; closelist(); push!(out, inline(l)); end
        i += 1
    end; closelist(); join(out, "\n")
end
body = join([convert_md(read(joinpath(DOC, "src", p), String);
                        srcfile = joinpath(DOC, "src", p)) for p in PAGES], "\n\\clearpage\n")
body = replace(body, r"(?<=\n)\*([^*$]+?)\*" => s"\\emph{\1}")
tex = """
\\documentclass[11pt,a4paper]{article}
\\usepackage{amsmath,amssymb,graphicx,geometry,hyperref,fontspec}
\\geometry{margin=2.2cm}
\\setlength{\\parindent}{0pt}\\setlength{\\parskip}{4pt}
\\title{$TITLE}
\\date{\\today}
\\begin{document}
\\maketitle
\\tableofcontents
\\clearpage
$body
\\end{document}
"""
pdfdir = joinpath(DOC, "pdf"); mkpath(pdfdir)
texfile = joinpath(pdfdir, OUTNAME * ".tex"); write(texfile, tex)
cd(pdfdir) do
    for _ in 1:2
        run(pipeline(`xelatex -interaction=nonstopmode -halt-on-error $(OUTNAME).tex`, stdout = "xelatex.log"))
    end
end
println("written: ", joinpath(pdfdir, OUTNAME * ".pdf"))
