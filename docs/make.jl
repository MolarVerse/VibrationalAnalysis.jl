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
		"Quickstart" => "quickstart.md",
		"Input Files" => "input-files.md",
		"CLI" => "cli.md",
		"Examples" => "examples.md",
		"Validation" => "validation.md",
		"Migration" => "migration.md",
		"License" => "license.md",
		"API" => "api.md",
	],
)

deploydocs(
	repo = "github.com/MolarVerse/VibrationalAnalysis.jl.git",
)
