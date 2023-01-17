push!(LOAD_PATH, "../src/")

using CausticsEngineering
using Documenter, DocStringExtensions

#
# HTML docs
#
Documenter.makedocs(
    modules = [CausticsEngineering],
    format = Documenter.HTML(),
    build = "build_html",
    sitename = "Amazing Caustics Back-Engineering!",
    pages = ["Table of Contents" => "table_of_contents.md", "Index" => "index.md"],

    # root = "./doc",
    # source = ".",
    # clean = true,
    # doctest = true,
    # repo = "",
    # highlightsig = true,
    # expand = [],
)


#
# pdf docs
#
Documenter.makedocs(
    modules = [CausticsEngineering],
    format = Documenter.LaTeX(),
    build = "build_pdf",
    sitename = "Amazing Caustics Back-Engineering!",
    pages = ["Table of Contents" => "table_of_contents.md", "Index" => "index.md"],
)
