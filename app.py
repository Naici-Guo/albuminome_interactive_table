import pandas as pd
from shiny import render, ui, App, reactive
from shared import data_index
from shared import data_main

other_proteins_options = sorted(data_index['Other_proteins'].dropna().str.split(",").explode().str.strip().unique())
other_proteins_options.remove('Albumin')
albumin_only_options = ["All", "Yes", "No"]

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.h3("About:"),
        ui.p("Albumin can bind to a wide range of proteins and peptides in human blood, leading to the concept of the 'albuminome' as defined by Stanley et al. (2004)."),
        ui.p("In quantitative proteomics studies, albumin depletion is a common preprocessing step. This is typically done using commercial kits that target high-abundance proteins (HAPs), including albumin. However, because albumin binds to many other proteins, this process can unintentionally remove additional proteins bound to albumin or other HAPs."),
        ui.p("As a result, it is often unclear whether co-removed proteins are due to albumin binding or associations with other HAPs. These ambiguities make it difficult to interpret depletion results without further analysis."),
        ui.p("In this app, we include studies of two types:"),
        ui.tags.ul(
            ui.tags.li("Studies investigating direct albumin-binding proteins/peptides (the albuminome)."),
            ui.tags.li("Studies that analyze proteins unintentionally removed alongside albumin and other high-abundance proteins (HAPs), " \
            "where the exact binding source is unclear.")),  
        ui.p("The results from these studies are summarized in an interactive table. Use this app to explore the albuminome and proteins that may be unintentionally removed during depletion protocols."),

        ui.h3("Filters to select:"),
         ui.p("1. Please select the type of studies you want to include based on the proteins of interest. "),
                    ui.tags.ul(
                        ui.tags.li("Select 'All' to include all studies, regardless of whether they focus on albumin or other HAPs."),
                        ui.tags.li("Select 'Yes': to include studies that specifically investigate albumin-binding proteins (the albuminome)."),
                        ui.tags.li("Select 'No' to include studies that investigate proteins co-removed with albumin and other HAPs, other HAPs can be further selected below."),
                    ), 
        ui.input_select("albumin_only", 
                    ui.strong("Albumin only study?"),
                    choices = albumin_only_options,
                    selected="Yes"),
        ui.p("2. If you selected ‘No’ above, you can now choose specific high-abundance proteins (HAPs) that you’re interested in. " \
        "This will allow you to explore proteins that are potentially co-removed along with the selected HAP(s) during depletion."),
        ui.p("If you selected ‘All’ above, all HAPs listed below will be included automatically in the results."),
        ui.input_checkbox_group("other_proteins", 
                            ui.strong("Other HAPs:" ), 
                            choices=other_proteins_options, 
                            selected="All"),
        width="600px",
        bg="#f8f8f8"
    ),
    # this card returns the list of studies that match the selected criteria
    ui.card(
        ui.h3("Details of included studies:"),
        ui.output_table("selected_papers"),
        style="""
            max-height: 400px; 
            overflow-y: auto; 
            overflow-x: auto;
        """,
        class_="left-align-table"
    ),
    # this card returns the aggregated table of proteins that match the selected criteria
    ui.card(
        ui.h3("Summary of proteins identified as albuminome and/or co-removed with selected HAPs:"),
        ui.output_table("aggregated_table"),
        full_screen=True,
        style="""
            max-height: 800px;
            overflow-y: auto; 
            overflow-x: auto;
        """,
        class_="left-align-table"
    ),
    ui.include_css("style.css"),
    title = ui.div(
        ui.h3("Comprehensive Profiling of Albuminome and High Abundance Protein Interactions in Human Blood",
            style="font-weight: 400; border-bottom: 1px solid #eee; padding-bottom: 0.5rem;"),
        ui.p("Exploring proteins co-removed with high abundance proteins (HAPs) in depletion-based proteomics studies",
            style="font-style: italic; color: #666; font-size: 0.95rem; margin-top: 0.5rem;"),
        style="text-align: center;")
)

# Define Server Logic
def server(input, output, session):

    @reactive.Calc
    def filtered_index():
        """Filter the data_index DataFrame based on user inputs."""
        filtered = data_index.copy()

        # Filter by albumin_only
        if input.albumin_only() != "All":
            filtered = filtered[filtered["Albumin_only"] == input.albumin_only()]

        # Filter by other_proteins if albumin_only is "No"
        if input.albumin_only() == "No" and input.other_proteins():
            selected_proteins = input.other_proteins()
            mask = filtered["Other_proteins"].dropna().str.split(",").apply(
                lambda proteins: any(p.strip() in selected_proteins for p in proteins)
            )
            filtered = filtered[mask]

        return filtered

    @output
    @render.table
    def selected_papers():
        """Display detailed study information in a table format."""
        studies = filtered_index()

        if studies.empty:
            # Return an empty table with a message
            return pd.DataFrame({"Message": ["No studies match the selected criteria."]})

        # Select relevant columns for display
        selected_columns = [
            "Paper",  # Study reference
            "Other_proteins",  # Full list of other proteins
            "Sample",  # Biological sample description
            "Separation_method",  # Separation method used
            "Detection_method",  # Detection method used
        ]

        # Return the selected data as a DataFrame
        return studies[selected_columns].reset_index(drop=True)

    @reactive.Calc
    def filtered_main():
        """Filter the data_main DataFrame based on selected studies in data_index."""
        selected_labels = filtered_index()["Label"].astype(str).tolist()

        # Subset data_main to include only relevant columns
        relevant_columns = ["Protein name", "Protein Uniprot ID"] + selected_labels
        return data_main.loc[:, data_main.columns.isin(relevant_columns)]

    @output
    @render.table
    def aggregated_table():
        """Aggregate protein data and rank by frequency of mentions."""
        # Get the selected paper labels
        selected_labels = filtered_index()["Label"].astype(str).tolist()

        if not selected_labels:
            # Return an empty table with a message if no papers are selected
            return pd.DataFrame({"Message": ["No proteins match the selected criteria."]})

        # Filter `data_main` to include only relevant columns
        relevant_columns = ["Protein", "Uniprot ID", "Protein Name"] + selected_labels
        filtered_main = data_main[relevant_columns]

        # Melt the data to create a long format for counting
        melted = filtered_main.melt(
            id_vars=["Protein", "Uniprot ID", "Protein Name"] ,
            value_vars=selected_labels,
            var_name="Paper",
            value_name="Mentioned",
        )

        # Filter for proteins that are mentioned (value = 1)
        mentioned_proteins = melted[melted["Mentioned"] == 1]

        # Group by Protein name and Uniprot ID to count occurrences
        aggregated = (
            mentioned_proteins.groupby(["Protein", "Uniprot ID", "Protein Name"] )
            .size()
            .reset_index(name="Count")
        )

        # Sort by count in descending order
        aggregated = aggregated.sort_values(by="Count", ascending=False)

        return aggregated
    
app = App(app_ui, server)