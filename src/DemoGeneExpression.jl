module DemoGeneExpression
using DataFrames
using CSV
using StatsBase
using StatsPlots
using PlotlyJS
using Mux
using Interact

function __axes(x::AbstractVector)
    labels = sort!(unique(x))
    (labels, Dict(zip(labels, 1:length(labels))))
end

"""
    spm(x)

Return the Specificity measure
"""
function spm(x)
    x² = x .^ 2
    s  = sum(x²)
    x² ./ (x .* s) 
end

"""
    tsi(x)

Return the Tissue specific index
"""
function tsi(x)
    s = sum(x)
    x ./ sum(x)
end

"""
    data(func::Function = zscore, name::String = "rna_tissue_hpa_TPM_filtered_stack.csv")

Compute the data needed by the heatmap
"""
function data(func::Function = zscore, name::String = "rna_tissue_hpa_TPM_filtered_stack.csv")
	df = DataFrame(CSV.File(
        isfile(name) ? name : joinpath(@__DIR__, "..", "data", name)
    ))
	labgene,   ordgene   = __axes(df.Gene)
	labtissue, ordtissue = __axes(df.Tissue)

	mat       = Array{Float64, 2}(undef, (length(labtissue), length(labgene)))
	for gene ∈ groupby(df, :Gene)
		igene = ordgene[gene.Gene[begin]]
		for (tissue, y) ∈ zip(gene.Tissue, func(gene.TPM))
			mat[ordtissue[tissue], igene] = y
		end
	end
    (; x = labgene, y = labtissue, z = mat)
end

"""
    heatmap(func::Function = zscore, name::String = "rna_tissue_hpa_TPM_filtered_stack.csv"; heatmapopts...)
    heatmap(data::NamedTuple; heatmapopts...)

Creates a heatmap from the data provided.
"""
function heatmap(data::NamedTuple; kwa...)
    plotlyjs()
    opts = Dict{Symbol, Any}(kwa)
    get!(opts, :xrotation, 90)
    get!(opts, :colorbar,  false)
    get!(opts, :size,      (1000, 800))
    StatsPlots.heatmap(
        data.z;
        xticks    = (1:length(data.x), data.x),
        yticks    = (1:length(data.y), data.y),
        opts...
    )
end

heatmap(args...; kwa...) =  heatmap(data(args...); kwa...)

function ui()
    files = Dict(
        "HPA database" => "hpa",
        "GTEX database" => "gtex",
    )

    scores = Dict(
        "z-score" => "zscore",
        "Tissue-specific index" => "tsi",
        "Specificity measure" => "spm"
    )

    dbinfo = Dict(
        "hpa" => (
            HTML("""<h2 class="subtitle">The Human Protein Atlas (HPA)</h2>"""),
            vskip(1em),
            md"""
            A Swedish-based program started in
            2003 with the aim to map all the human proteins in cells, tissues
            and organs using integration of various omics technologies,
            including antibody-based imaging, mass spectrometry-based
            proteomics, transcriptomics and systems biology. All the data in
            the knowledge resource is open access to allow scientists both in
            academia and industry to freely access the data for exploration of
            the human proteome.""",
            vskip(2em),
        ),
        "gtex" => (
            HTML("""<h2 class="subtitle">The Genotype-Tissue Expression (GTEx)</h2>"""),
            vskip(1em),
            md"""

            This project is an ongoing effort to build a comprehensive public
            resource to study tissue-specific gene expression and regulation.
            Samples were collected from 54 non-diseased tissue sites across
            nearly 1000 individuals, primarily for molecular assays including
            WGS, WES, and RNA-Seq. Remaining samples are available from the
            GTEx Biobank. The GTEx Portal provides open access to data
            including gene expression, QTLs, and histology images.""",
            vskip(2em),
        )
    )

    scoreinfo = Dict(
        "zscore" => (
            HTML("""<h2 class="subtitle">Z-score</h2>"""),
            vskip(1em),
            md"""
            Supposing a distribution is normal, the score is the number of standard deviations to the mean.

            $$µ = \frac{\sum{x_i}}{N}$$

            $$σ = \frac{\sum{(x_i-µ)^2}}{N}$$

            $$zscore(x) = \frac{x-µ}{σ}$$
            """,
            vskip(2em),
        ),
        "tsi" => (
            HTML("""<h2 class="subtitle">Tissue-specific index</h2>"""),
            vskip(1em),
            md"$$tsi(x) = x\sum{x_i}$$",
            vskip(2em),
        )
    )

    precomputed = Dict{Tuple{String, String}, NamedTuple}()

    vbox(
        vskip(1em),
        HTML("<h1 class=\"title\">Gene Expression</h1>"),
        begin
            @manipulate for var"Database" ∈ files, var"Score type" ∈ scores
                db, score = var"Database", var"Score type"
                vbox(
                    hbox(
                        vbox(get(scoreinfo, score, ("",))...),
                        hskip(1em),
                        vbox(get(dbinfo, db, ("",))...),
                    ),
                    vskip(1em),
                    heatmap(
                        get!(precomputed, (db, score)) do
                            data(
                                eval(Symbol(score)),
                                "rna_tissue_$(db)_TPM_filtered_stack.csv"
                            )
                        end
                    )
                )
            end
        end
    )
end

function launch()
  port = rand(8000:9000)

  @info "http://localhost:$port"
  WebIO.webio_serve(page("/", req -> ui()), port) # serve on a random port
end

end # module
