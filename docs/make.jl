using Documenter
using VibrationalAnalysis

makedocs(
	sitename = "VibrationalAnalysis.jl",
	format = Documenter.HTML(
		prettyurls = (get(ENV, "PRETTY_URLS", nothing) == "true"),
	),
	pages = [
		"Home" => "index.md",
		"Installation" => "installation.md",
		"Examples" => "examples.md",
		"License" => "license.md",
		"API" => "api.md",
	],
)
