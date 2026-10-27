# =============================================================================
# docs/make_pdf.jl - build the PDF Beam Physics Manual from the SAME Markdown
# =============================================================================
# WHAT : small Markdown -> LaTeX converter for the subset used in docs/src
#        (headings, paragraphs, math blocks, ``inline math``, `code`, **bold**,
#        lists, pipe tables, figures with captions), then xelatex (twice).
# WHY  : no pandoc on this machine; one source for HTML and PDF.
# HOW  : julia docs/make_pdf.jl  -> docs/pdf/LieRing_Beam_Physics_Manual.pdf
# NOTE : GIF animations cannot be embedded in a PDF; they are replaced by a note.
# =============================================================================
const DOC = @__DIR__
const PAGES = ["index.md", "coordinates.md", "transform.md", "elements.md", "linear.md",
               "second_order.md", "chromaticity.md", "rf.md", "acceleration.md", "longitudinal.md",
               "symplectic.md", "sextupoles.md", "assumptions.md", "codemap.md", "programmer.md"]

esc_text(s) = replace(s, "\\" => "\\textbackslash{}", "&" => "\\&", "%" => "\\%", "#" => "\\#",
                      "_" => "\\_", "\$" => "\\\$", "^" => "\\^{}", "~" => "\\~{}", "{" => "\\{", "}" => "\\}")

function inline(s)
    out = IOBuffer()
    parts = split(s, "``")
    for (k, p) in enumerate(parts)
        if iseven(k)
            print(out, "\$", p, "\$")
        else
            q = replace(p, r"\[([^\]]*)\]\((@ref|[\w.]+\.md)\)" => s"\1")
            cs = split(q, "`")
            for (j, c) in enumerate(cs)
                if iseven(j)
                    print(out, "\\texttt{", esc_text(c), "}")
                else
                    t = esc_text(c)
                    t = replace(t, r"\*\*([^*]+)\*\*" => s"\\textbf{\1}")
                    t = replace(t, r"(?<![\w*])\*([^*\s][^*]*)\*(?![\w*])" => s"\\emph{\1}")
                    print(out, t)
                end
            end
        end
    end
    String(take!(out))
end

function convert(md)
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
        elseif startswith(l, "## ")
            closelist(); push!(out, "\\subsection{" * inline(l[4:end]) * "}")
        elseif startswith(l, "# ")
            closelist(); push!(out, "\\section{" * inline(l[3:end]) * "}")
        elseif startswith(l, "|")
            closelist(); rows = String[]
            while i <= length(L) && startswith(L[i], "|"); push!(rows, L[i]); i += 1; end
            cells(r) = strip.(split(strip(r)[2:end-1], "|"))
            hdr = cells(rows[1]); n = length(hdr)
            w = round(0.95 / n, digits = 3)
            push!(out, "\\begin{center}\\small\\begin{tabular}{" * join(fill("p{$(w)\\linewidth}", n), "") * "}\\hline")
            push!(out, join(["\\textbf{" * inline(c) * "}" for c in hdr], " & ") * "\\\\\\hline")
            for r in rows[3:end]
                push!(out, join(inline.(cells(r)), " & ") * "\\\\")
            end
            push!(out, "\\hline\\end{tabular}\\end{center}"); continue
        elseif (m = match(r"^!\[(.*)\]\((.*)\)\s*$", l)) !== nothing
            closelist()
            cap, path = m.captures
            if endswith(path, ".gif")
                push!(out, "\\begin{quote}\\emph{Animation " * esc_text(path) * " (see the HTML manual or the GIF file): " * inline(cap) * "}\\end{quote}")
            else
                fp = replace(joinpath(DOC, "src", path), "\\" => "/")
                push!(out, "\\begin{figure}[htbp]\\centering\\includegraphics[width=\\linewidth]{" * fp)
                out[end] *= "}\\caption{" * inline(cap) * "}\\end{figure}"
            end
        elseif (m = match(r"^(\d+)\.\s+(.*)$", l)) !== nothing
            list != "enumerate" && (closelist(); push!(out, "\\begin{enumerate}"); list = "enumerate")
            push!(out, "\\item " * inline(m.captures[2]))
        elseif (m = match(r"^\*\s+(.*)$", l)) !== nothing
            list != "itemize" && (closelist(); push!(out, "\\begin{itemize}"); list = "itemize")
            push!(out, "\\item " * inline(m.captures[1]))
        elseif startswith(l, "   ") && list != "" && strip(l) != ""
            out[end] *= " " * inline(strip(l))
        elseif strip(l) == ""
            closelist(); push!(out, "")
        else
            closelist(); push!(out, inline(l))
        end
        i += 1
    end
    closelist()
    join(out, "\n")
end

body = join([convert(read(joinpath(DOC, "src", p), String)) for p in PAGES], "\n\\clearpage\n")
body = replace(body, r"(?<=\n)\*([^*$]+?)\*" => s"\\emph{\1}")   # *italic paragraph*
tex = """
\\documentclass[11pt,a4paper]{article}
\\usepackage{amsmath,amssymb,graphicx,geometry,hyperref,fontspec}
\\geometry{margin=2.2cm}
\\setlength{\\parindent}{0pt}\\setlength{\\parskip}{4pt}
\\title{LieRing Beam Physics Manual\\\\\\large symplectic ring modelling with Lie maps (v0.4)}
\\date{\\today}
\\begin{document}
\\maketitle
\\tableofcontents
\\clearpage
$body
\\end{document}
"""
pdfdir = joinpath(DOC, "pdf"); mkpath(pdfdir)
texfile = joinpath(pdfdir, "LieRing_Beam_Physics_Manual.tex")
write(texfile, tex)
cd(pdfdir) do
    for _ in 1:2
        run(pipeline(`xelatex -interaction=nonstopmode -halt-on-error LieRing_Beam_Physics_Manual.tex`, stdout = "xelatex.log"))
    end
end
println("written: ", joinpath(pdfdir, "LieRing_Beam_Physics_Manual.pdf"))
