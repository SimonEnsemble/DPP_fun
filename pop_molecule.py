import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem.Draw import rdDepictor
    from marimo_chem_utils import (
        add_fingerprint_column,
        add_image_column,
        add_inchi_key_column,
        add_tsne_columns,
        interactive_chart
    )
    from rdkit.Chem import PandasTools
    from rdkit.Chem import Descriptors
    import useful_rdkit_utils as uru
    import altair as alt

    return PandasTools, alt, mo, np, pd, rdDepictor, uru


@app.cell
def _(pd):
    df = pd.read_csv('pca_data.csv')
    df
    return (df,)


@app.cell
def _(PandasTools, df, rdDepictor, uru):
    data = df.copy()
    PandasTools.AddMoleculeColumnToFrame(data, smilesCol="molecule")
    data.ROMol.apply(rdDepictor.Compute2DCoords)
    data["image"] = data.ROMol.apply(uru.mol_to_base64_image, target="altair")
    data = data.drop(columns=["ROMol"])
    # descs = [Descriptors.CalcMolDescriptors(m) for m in data.ROMol]
    # data = df.join(pd.DataFrame(descs))
    return (data,)


@app.cell
def _(data, mo):
    x = mo.ui.dropdown(options=data.columns.drop("image"), value="pca_1")
    x
    return (x,)


@app.cell
def _(data, mo):
    y = mo.ui.dropdown(options=data.columns.drop("image"), value="pca_2")
    y
    return (y,)


@app.cell
def _():
    pca_range = [-1, 1]
    return (pca_range,)


@app.cell
def _(alt, data, mo, pca_range, x, y):
    mo.ui.altair_chart(alt.Chart(data).mark_point().encode(
        x=alt.X(x.value, scale=alt.Scale(domain=pca_range)), 
        y=alt.Y(y.value, scale=alt.Scale(domain=pca_range)),
        tooltip=alt.Tooltip(['image']))
                      )
    return


@app.cell
def _(data):
    smell_list = data.columns.drop("image").drop("molecule").drop("pca_1").drop("pca_2")
    return


@app.cell
def _():
    interesting_smells = ["vanilla", "mint", "fish", "wine", "chamomile", "popcorn", "currant"]
    return (interesting_smells,)


@app.cell
def _(interesting_smells, mo):
    smell = mo.ui.dropdown(options=interesting_smells, value="wine") #  np.sort(smell_list)
    smell
    return (smell,)


@app.cell
def _(alt, data, mo, smell, x, y):
    other_points = alt.Chart(data[data[smell.value] == 0]).mark_point(
        color="lightgray",
        opacity=0.8
    ).encode(
        x=alt.X(x.value, scale=alt.Scale(domain=[-1, 1])),
        y=alt.Y(y.value, scale=alt.Scale(domain=[-1, 1])),
        tooltip=["image", smell.value]
    )

    smell_points = alt.Chart(data[data[smell.value] == 1]).mark_point(
        color="red",
        opacity=1.0
    ).encode(
        x=alt.X(x.value, scale=alt.Scale(domain=[-1, 1])),
        y=alt.Y(y.value, scale=alt.Scale(domain=[-1, 1])),
        tooltip=["image", smell.value]
    )

    chart = (other_points + smell_points).properties(
        width=400,
        height=400
    )

    mo.ui.altair_chart(chart)
    return


@app.cell
def _(alt, mo, pca_range, smell, x, y):
    def selected_molecules(data, selected_ids):
        plot_data = data.copy()
        plot_data["selected"] = 0
        plot_data.loc[selected_ids, "selected"] = 1

        plot_x = alt.X(x.value, scale=alt.Scale(domain=pca_range))
        plot_y = alt.Y(y.value, scale=alt.Scale(domain=pca_range))

        other_molecules = alt.Chart(plot_data[plot_data["selected"] == 0]).mark_point(
            color="lightgray",
            opacity=0.8
        ).encode(
            x=plot_x,
            y=plot_y,
            tooltip=["image"]
        )

        selected_data = plot_data[plot_data["selected"] == 1]
        selected_molecules = alt.Chart(selected_data).mark_point(
            color="blue",
            opacity=1.0
        ).encode(
            x=plot_x,
            y=plot_y,
            tooltip=["image"]
        )

        smell_points = alt.Chart(selected_data[selected_data[smell.value] == 1]).mark_point(
            color="red",
            opacity=1.0
        ).encode(
            x=plot_x,
            y=plot_y,
            tooltip=["image", smell.value]
        )

        return mo.ui.altair_chart(other_molecules + selected_molecules + smell_points)

    return (selected_molecules,)


@app.cell
def _(np):
    with open("ids_molecule_dpp.txt", "r") as f:
        ids_molecule_dpp = list(map(int, f.read().split()))
    ids_molecule_dpp = np.array(ids_molecule_dpp) - 1 
    return (ids_molecule_dpp,)


@app.cell
def _(data, ids_molecule_dpp, selected_molecules):
    selected_molecules(data, ids_molecule_dpp)
    return


@app.cell
def _(np):
    with open("ids_molecule_uniform.txt", "r") as f_uniform:
        ids_molecule_uniform = list(map(int, f_uniform.read().split()))
    ids_molecule_uniform = np.array(ids_molecule_uniform) - 1 
    return (ids_molecule_uniform,)


@app.cell
def _(data, ids_molecule_uniform, selected_molecules):
    selected_molecules(data, ids_molecule_uniform)
    return


@app.cell
def _(data, ids_molecule_dpp, ids_molecule_uniform, interesting_smells, pd):
    percent_smell = pd.DataFrame({
        "all": data[interesting_smells].sum(axis=0) / len(data),
        "DPP": data.loc[ids_molecule_dpp, interesting_smells].sum(axis=0) / len(ids_molecule_dpp),
        "uniform": data.loc[ids_molecule_uniform, interesting_smells].sum(axis=0) / len(ids_molecule_uniform),
    })
    percent_smell
    return (percent_smell,)


@app.cell
def _(percent_smell):
    (percent_smell > 0).sum(axis=0)
    return


if __name__ == "__main__":
    app.run()
