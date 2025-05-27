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
        "Albumin can bind to a wide quantity of proteins and peptides in humnan blood, which leads to the term 'albuminome'" \
        "defined by Stanley et al. (2004). In quantitative proteomics studies, albumin-depletion, which is the standard procedure" \
        "can also remove other proteins due to the binding of albumin and other proteins. " \
        "Removal of albumin in proteomics studies is usually done by commercialised kits that target at other high abundant " \
        "proteins (HAPs), so that the co-removed proteins can possibly due to co-removal with albumin or other HAPs that cannot and have never" \
        "be distinguished. " \
        "We include studies that are of two types: (1) studies that are studying albumin-binding proteins only, and "
        "(2) studies that are studying co-removed proteins with albumin or other HAPs. And summerised the results in this interactive table. " \
        "This app allows you to explore the albuminome and/or co-removed proteins ",
        ui.strong("Select HAPs of interest:"),
        #ui.h3("Filters to select"),
        ui.input_select("albumin_only", 
                    ui.strong("If you are only interested in proteins that potentially bind with albumin, i.e., included studies are " \
                    "studying albumin-binding only, and please select 'Yes',:",\
                    "if you are interested in all proteins, select 'All'." \
                    " If you are interested in proteins that are co-removed with albumin, select 'No'."), 
                    choices = albumin_only_options,
                    selected="Yes"),
        ui.input_checkbox_group("other_proteins", 
                            ui.strong("If selected 'No' above, you can select cofounding other HAPs proteins that are usually co-removed with albumin and bind with. " \
                                      "other proteins" ), 
                            choices=other_proteins_options, 
                            selected="All"),
        width="600px",
        bg="#f8f8f8"
    ),
    # this card returns the list of studies that match the selected criteria
    ui.card(
        ui.h3("Details about selected studies"),
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
        ui.h3("Albuminome and/or co-removed proteins with HAPs"),
        ui.output_table("aggregated_table"),
        full_screen=True,
        style="""
            overflow-y: auto; 
            overflow-x: auto;
        """,
        class_="left-align-table"
    ),
    ui.include_css("style.css"),
    #title= "Albuminome and/or co-removed proteins with high abundance proteins in human blood based studies",
    #title = ui.div(
        #ui.h2("Albuminome Research Explorer", style="margin-bottom: 0.2rem;"),
        #ui.p("Studying albumin-binding and co-removed proteins ", ui.strong("(HAPs)"), " in human blood",
        #     style="font-size: 1.1rem; margin-top: 0; color: #555;"),
        #style="text-align: center; padding: 0.5rem;")
    title = ui.div(
        ui.h3("Characterization of Albuminome and High-Abundance Protein Interactions",
            style="font-weight: 400; border-bottom: 1px solid #eee; padding-bottom: 0.5rem;"),
        ui.p("Analysis of albumin-binding proteins and co-removed HAPs in human serum studies",
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

"""
with ui.sidebar(bg="#f8f8f8", width="500px"):
    ui.h3("About:")
    "This app allows you to explore the albuminome and/or co-removed proteins "
    "with high abundance proteins (HAPs) in human blood-based studies. You can select specific proteins of interest, "
    "view their abundance, and filter by species."
    ui.strong("Select HAPs of interest:") 
    ui.input_select("choose_albumin_only", 
                    ui.strong("Albumin only?"), 
                    choices=["Yes", "No"])
    ui.input_checkbox_group("choose_other_proteins", 
                            ui.strong("If select No above, proteins co-removed other HAPs:"), 
                            choices=other_proteins_options, 
                            selected="All")
"""