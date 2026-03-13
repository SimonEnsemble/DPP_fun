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

    return PandasTools, alt, mo, pd, rdDepictor, uru


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
    data = data[["pca_1", "pca_2", "image"]]
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
def _(alt, data, mo, x, y):
    mo.ui.altair_chart(alt.Chart(data).mark_point().encode(
        x=x.value,
        y=y.value,
        tooltip=alt.Tooltip(['image']))
                      )
    return


if __name__ == "__main__":
    app.run()
