# Shared console layout for Bell program examples.
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# banner("title")     top of the run
# section("title")    blank line, then a short heading, then a blank line
# note("...")         indented remark under a section

function banner(title)
    println()
    println("="^72)
    println("  ", title)
    println("="^72)
    println()
end

function section(title, notes...)
    println()
    println("--- ", title)
    println()
    for line in notes
        println("  ", line)
    end
    isempty(notes) || println()
end

note(text) = println("  ", text)
