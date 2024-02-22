import pandas as pd
import plotly.express as px
import os

# File with counts per species as .csv
infile = "/Users/simplexdna/Library/CloudStorage/GoogleDrive-christopher.hempel@kaust.edu.sa/.shortcut-targets-by-id/1c91orCwfstL7NmlFbhdJ8Ssz87pDMHx-/SIREN project/Red Sea species list processing - chris - NCBI + BOLD counts/red_sea_species_list_standardized_synonyms_merged_unaccepted_merged.csv"
# Graph output directory
plot_outdir = (
    "/Users/simplexdna/Library/CloudStorage/GoogleDrive-christopher.hempel@kaust.edu.sa/.shortcut-targets-by-id/1c91orCwfstL7NmlFbhdJ8Ssz87pDMHx-/SIREN project/coi_counts_graphs_NCBI_and_BOLD"
)
# Min number of available COI seq for taxon to count as found
min_seqs = 1
# Maximum number of taxa to show in barplots
num_taxa_barplots = 50
# Determine ranks to visualize
# ranks = ["kingdom", "phylum", "class", "order", "family", "genus"]
ranks = ["phylum"]
# Kingdom filter (all listed will be kept)
# kingdoms = ["Animalia", "Archaea", "Bacteria", "Chromista", "Plantae"]
kingdoms = ["Animalia"]

# Make plot directory
os.makedirs(plot_outdir, exist_ok=True)

# Import data
df = pd.read_csv(infile)

# Filter the DataFrame based on the kingdom filter
df = df[df["kingdom"].isin(kingdoms)]

# Turn to p/a based on min_seqs
df["COI sequence available GenBank"] = df["COI_sequences_on_GenBank"].apply(
    lambda x: 1 if x >= min_seqs else 0
)
df["COI sequence available BOLD"] = df["COI_sequences_on_BOLD"].apply(
    lambda x: 1 if x >= min_seqs else 0
)

# For both databases:
for db in ["GenBank", "BOLD"]:

    # Strip plots

    ## Loop over ranks
    for rank in ranks:
        ## Group counts of ranks
        rank_series = df.groupby(rank)[f"COI sequence available {db}"].sum()
        ## Make separate violin plot
        fig = px.strip(rank_series, width=300).update_traces(jitter=1)
        median = rank_series.median()
        if int(median) == median and isinstance(median, float):
            median = round(median)
        mean = round(rank_series.mean(), 1)
        fig.update_xaxes(showticklabels=False)
        fig.update_layout(
            xaxis_title=f"{rank} - median: {median} - mean: {mean}",
            yaxis_title=f"# of taxa with ≥ {min_seqs} COI seqs on {db}",
        )
        fig.show()
        fig.write_image(
            os.path.join(plot_outdir, f"stripplot_{db}_{rank}_{min_seqs}_seqs.png"), width=300
        )

    if "pecies" in ranks[-1]:
        species_series = df.groupby(ranks[-1])[f"COI sequence available {db}"].sum()
        species_df = pd.DataFrame(
            {
                "spacer": ["a", "a"],
                "category": ["without", "with"],
                "count": [(species_series == 0).sum(), (species_series == 1).sum()],
            }
        )
        fig = px.bar(
            species_df,
            text_auto=True,
            x="spacer",
            y="count",
            color="category",
            title=f"Species with and without ≥ {min_seqs} COI seqs",
        )
        fig.update_xaxes(showticklabels=False)
        fig.update_layout(
            xaxis_title="",
            yaxis_title="",
        )
        fig.update_traces(insidetextanchor="middle")
        fig.show()
        fig.write_image(
            os.path.join(plot_outdir, f"stackedbarplot_species_{db}_{min_seqs}_seqs.png")
        )

    # Bar plots
    ## Loop over ranks
    for rank in ranks:
        ### Group counts of ranks and count how many taxa are grouped
        rank_df = df.groupby(rank).agg({f"COI sequence available {db}": ["sum", "count"]})
        rank_df.columns = [
            f"COI sequence available {db}",
            "Taxon count",
        ]
        ### Determine number of taxa with missing COI sequences
        rank_df[f"COI sequence not available {db}"] = (
            rank_df["Taxon count"] - rank_df[f"COI sequence available {db}"]
        )
        ### Determine relative number of counts
        rank_df["proportion"] = (
            rank_df[f"COI sequence available {db}"] / rank_df["Taxon count"] * 100
        )

        # Save df
        rank_df.to_csv(
            f"/Users/simplexdna/Library/CloudStorage/GoogleDrive-christopher.hempel@kaust.edu.sa/.shortcut-targets-by-id/1c91orCwfstL7NmlFbhdJ8Ssz87pDMHx-/SIREN project/Red Sea species list processing - chris - NCBI + BOLD counts/red_sea_species_list_results_{db}_{rank}_{min_seqs}_seqs.csv"
        )

        ### Define barplot title based on number of taxa
        if len(rank_df) >= num_taxa_barplots:
            barplot_title_total = f"Total number of taxa with ≥ {min_seqs} COI seqs available on {db} - {rank} - first {num_taxa_barplots} taxa"
            barplot_title_relative = f"Relative number of taxa with ≥ {min_seqs} COI seqs available on {db} - {rank} - first {num_taxa_barplots} taxa"
        else:
            barplot_title_total = f"Total number of taxa with ≥ {min_seqs} COI seq available on {db} - {rank}"
            barplot_title_relative = f"Relative number of taxa with ≥ {min_seqs} COI seqs available on {db} - {rank}"

        ### Total barplots
        #### Sort and keep top X taxa, where X = num_taxa_barplots
        rank_df_total = rank_df[f"COI sequence available {db}"].sort_values(ascending=False)[
            :num_taxa_barplots
        ]
        #### Make bar plots
        barplot = px.bar(rank_df_total)
        barplot.update_layout(
            xaxis_title="Taxa",
            yaxis_title="Counts",
            title={
                "text": barplot_title_total,
                "font": {"size": 15},
            },
            showlegend=False,
            width=max(300, 25 * len(rank_df_total)),
        )
        barplot.update_xaxes(tickangle=45)
        barplot.show()
        barplot.write_image(
            os.path.join(plot_outdir, f"barplot_total_{db}_{rank}_{min_seqs}_seqs.png")
        )

        ### Relative barplots
        #### Sort and keep top X taxa, where X = num_taxa_barplots
        rank_df_relative = rank_df["proportion"].sort_values(ascending=False)[
            :num_taxa_barplots
        ]
        #### Make bar plots
        barplot = px.bar(rank_df_relative)
        barplot.update_layout(
            xaxis_title="Taxa",
            yaxis_title="%",
            title={
                "text": barplot_title_relative,
                "font": {"size": 15},
            },
            showlegend=False,
            width=max(300, 25 * len(rank_df_relative)),
        )
        barplot.update_xaxes(tickangle=45)
        barplot.update_layout(yaxis_range=[0, 100])
        barplot.show()
        barplot.write_image(
            os.path.join(plot_outdir, f"barplot_relative_{db}_{rank}_{min_seqs}_seqs.png")
        )

        ### Stacked barplot
        stacked_df = rank_df.reset_index().sort_values("Taxon count", ascending=True)
        stacked_bar = px.bar(
            stacked_df,
            y=rank,
            x=[f"COI sequence available {db}", "Taxon count"],
            title=f"Total number of taxa on rank {rank} and number of taxa with ≥ {min_seqs} COI seqs available on {db}",
            barmode="group",
            height=700,
            width=500,
            text_auto=True,
            # log_x=True,
        )

        stacked_bar.update_layout(
            yaxis_title="Taxa",
            xaxis_title="Sequence counts",
            legend=dict(yanchor="bottom", y=0.01, xanchor="right", x=0.99),
        )
        stacked_bar.show()
        stacked_bar.write_image(
            os.path.join(plot_outdir, f"stacked_barplot_{db}_{rank}_{min_seqs}_seqs.png")
        )
        stacked_bar.write_image(
            os.path.join(plot_outdir, f"stacked_barplot_{db}_{rank}_{min_seqs}_seqs.svg")
        )

        ### Bubble plot
        bubble_df = rank_df.reset_index()
        bubble_plot = px.scatter(
            bubble_df,
            x=f"COI sequence available {db}",
            y=f"COI sequence not available {db}",
            size="Taxon count",
            color="proportion",
            hover_name=rank,
            size_max=60,
            color_continuous_scale="Plasma_r",
            width=550,
            height=550,
            title=f"{rank} with at least {min_seqs} available seqs - {db}",
            text=f"{rank}",
        )
        bubble_plot.update_layout(
            yaxis_range=[
                -20,
                1800,
            ],
            xaxis_range=[
                -100,
                1800,
            ],
        )
        bubble_plot.update_yaxes(nticks=4)
        bubble_plot.update_traces(marker=dict(line=dict(width=0)))
        bubble_plot.show()
        bubble_plot.write_image(
            os.path.join(plot_outdir, f"bubbleplot_{db}_{rank}_{min_seqs}_seqs.png")
        )
        bubble_plot.write_image(
            os.path.join(plot_outdir, f"bubbleplot_{db}_{rank}_{min_seqs}_seqs.svg")
        )
