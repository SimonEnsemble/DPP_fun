import marimo

__generated_with = "0.23.1"
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
    import matplotlib.pyplot as plt
    import seaborn as sns

    from sklearn.model_selection import train_test_split
    from sklearn.kernel_ridge import KernelRidge
    from sklearn.metrics import mean_squared_error

    from sklearn.svm import SVC
    from sklearn.multiclass import OneVsRestClassifier
    from sklearn.metrics import f1_score

    return (
        KernelRidge,
        OneVsRestClassifier,
        PandasTools,
        SVC,
        alt,
        f1_score,
        mean_squared_error,
        mo,
        np,
        pd,
        plt,
        rdDepictor,
        uru,
    )


@app.cell
def _(pd):
    df = pd.read_csv('pca_data.csv')
    df
    return (df,)


@app.cell
def _(PandasTools, df, rdDepictor, uru):
    data = df.copy()
    PandasTools.AddMoleculeColumnToFrame(data, smilesCol="molecule")#molecule
    data.ROMol.apply(rdDepictor.Compute2DCoords)
    data["image"] = data.ROMol.apply(uru.mol_to_base64_image, target="altair")
    data = data.drop(columns=["ROMol"])
    return (data,)


@app.cell
def _(data):
    data
    return


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
def _(data):
    smell_list = data.columns.drop("image").drop("molecule").drop("pca_1").drop("pca_2").drop("importance_pca_1").drop("importance_pca_2")
    return (smell_list,)


@app.cell
def _(smell_list):
    smell_list
    return


@app.cell
def _():
    interesting_smells = ["vanilla", "mint", "wine", "chamomile", "currant", "fish", "popcorn"]
    return (interesting_smells,)


@app.cell
def _(mo, np, smell_list):
    smell = mo.ui.dropdown(options=np.sort(smell_list), value="wine") #  np.sort(smell_list)
    smell
    return (smell,)


@app.cell
def _(alt, data, mo, smell, x, y):
    other_points = alt.Chart(data[data[smell.value] == 0]).mark_point(
        color="lightgray",
        opacity=0.8
    ).encode(
        x=alt.X(x.value, scale=alt.Scale(domain=[-1, 1]), title=f"{x.value} ({data['importance_' + x.value][0]:.2f})"),
        y=alt.Y(y.value, scale=alt.Scale(domain=[-1, 1]), title=f"{y.value} ({data['importance_' + y.value][0]:.2f})"),
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
    def selected_molecules(data, selected_ids, smell=smell.value):
        plot_data = data.copy()
        plot_data["selected"] = 0
        plot_data.loc[selected_ids, "selected"] = 1

        plot_x = alt.X(x.value, scale=alt.Scale(domain=pca_range), title=f"{x.value} ({data['importance_' + x.value][0]:.2f})")
        plot_y = alt.Y(y.value, scale=alt.Scale(domain=pca_range), title=f"{y.value} ({data['importance_' + y.value][0]:.2f})")

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

        # smell_points = alt.Chart(selected_data[selected_data[smell] == 1]).mark_point(
        #     color="red",
        #     opacity=1.0
        # ).encode(
        #     x=plot_x,
        #     y=plot_y,
        #     tooltip=["image", smell]
        # )

        return mo.ui.altair_chart(other_molecules + selected_molecules) #+ smell_points)

    return (selected_molecules,)


@app.cell
def _(pd):
    def get_selected_data(file):
        ids_molecule = pd.read_csv(file)
        return ids_molecule - 1

    return (get_selected_data,)


@app.cell
def _(ids_molecule_mcmc_dpp):
    n_runs = ids_molecule_mcmc_dpp.shape[1]
    return (n_runs,)


@app.cell
def _(mo, n_runs):
    slider = mo.ui.slider(1, n_runs, label="number of runs", value=1)
    return (slider,)


@app.cell
def _(mo, slider):
    mo.hstack([slider, mo.md(f"number of runs: {slider.value}")])
    return


@app.cell
def _(get_selected_data):
    ids_molecule_mcmc_dpp = get_selected_data('ids_molecule_mcmc_dpp_20.csv')
    return (ids_molecule_mcmc_dpp,)


@app.cell
def _(data, ids_molecule_mcmc_dpp, selected_molecules, slider):
    selected_molecules(data, ids_molecule_mcmc_dpp.iloc[:, slider.value-1])
    return


@app.cell
def _(get_selected_data):
    ids_molecule_uniform = get_selected_data('ids_molecule_uniform_20.csv')
    return (ids_molecule_uniform,)


@app.cell
def _(data, ids_molecule_uniform, selected_molecules, slider):
    selected_molecules(data, ids_molecule_uniform.iloc[:, slider.value-1])
    return


@app.cell
def _(data, np):
    def num_smells(ids_molecule, n_runs, smell):
        n_smells = np.mean([
            (data.loc[ids_molecule.iloc[:, i], smell].sum(axis=0) > 0).sum()
            for i in range(n_runs)
        ])
        return [n_smells]

    return (num_smells,)


@app.cell
def _(
    data,
    ids_molecule_mcmc_dpp,
    ids_molecule_uniform,
    n_runs,
    num_smells,
    pd,
    smell_list,
):
    n_smell = pd.DataFrame({
        "all": (data[smell_list].sum(axis=0) > 0).sum(),
        "mcmc_DPP": num_smells(ids_molecule_mcmc_dpp, n_runs, smell_list),
        "uniform": num_smells(ids_molecule_uniform, n_runs, smell_list),
    })
    n_smell
    return


@app.cell
def _(pd):
    def proportion_smell(data, ids, smells, n_runs):
        smell = pd.DataFrame(
            [data.loc[ids.iloc[:, i], smells].sum(axis=0) for i in range(n_runs)]
        ).mean()
        smell = smell / smell.sum()
        return smell

    return (proportion_smell,)


@app.cell
def _(data):
    data
    return


@app.cell
def _(ids_molecule_mcmc_dpp):
    ids_molecule_mcmc_dpp
    return


@app.cell
def _(dpp_smell):
    dpp_smell
    return


@app.cell
def _(
    data,
    ids_molecule_mcmc_dpp,
    ids_molecule_uniform,
    n_runs,
    pd,
    plt,
    proportion_smell,
    smell_list,
):
    dpp_smell = proportion_smell(data, ids_molecule_mcmc_dpp, smell_list, n_runs)
    uniform_smell = proportion_smell(data, ids_molecule_uniform, smell_list, n_runs)

    all_smell = data[smell_list].sum()
    all_smell = all_smell / all_smell.sum()

    pd.concat(
        [dpp_smell.rename("dpp"),
         uniform_smell.rename("uniform"),
         all_smell.rename("all data")],
        axis=1
    ).sort_values("all data").plot.bar()

    plt.xlabel("smell")
    plt.ylabel("proportion")
    plt.tight_layout()
    plt.show()
    return (dpp_smell,)


@app.cell
def _(np):
    def entropy(p):
        p = np.array(p)
        p = p[p > 0]  
        return -np.sum(p * np.log(p))

    return (entropy,)


@app.cell
def _(
    data,
    entropy,
    ids_molecule_mcmc_dpp,
    ids_molecule_uniform,
    n_runs,
    pd,
    proportion_smell,
    smell_list,
):
    all_data_smell = data[smell_list].sum()
    all_data_smell = all_data_smell / all_data_smell.sum()
    pd.DataFrame({
    "entropy": [
        entropy(proportion_smell(data, ids_molecule_mcmc_dpp, smell_list, n_runs)),
        entropy(proportion_smell(data, ids_molecule_uniform, smell_list, n_runs)),
        entropy(all_data_smell),
    ]}, index=["dpp", "uniform", "all data"])
    return


@app.cell
def _(OneVsRestClassifier, SVC, f1_score, np):
    def compute_micro_f1(train_idx, other_idx, K, target):
        all_idx = np.arange(len(K))
    
        combine_idx = np.unique(np.concatenate([train_idx, other_idx]))
        test_idx = np.setdiff1d(all_idx, combine_idx)

        K_train = K[np.ix_(train_idx, train_idx)]
        K_test = K[np.ix_(test_idx, train_idx)]

        target_train = target[train_idx]
        target_test = target[test_idx]

        model = OneVsRestClassifier(
            SVC(kernel="precomputed", class_weight="balanced")
        )
        model.fit(K_train, target_train)
        target_pred = model.predict(K_test)

        return f1_score(target_test, target_pred, average="macro", zero_division=0)

    return (compute_micro_f1,)


@app.cell
def _(np, pd):
    K = np.array(pd.read_csv('Gram_matrix.csv'))
    K
    return (K,)


@app.cell
def _():
    n_s_data_list = [5, 10, 20]
    n_s_data = len(n_s_data_list)
    return n_s_data, n_s_data_list


@app.cell
def _(
    K,
    compute_micro_f1,
    data,
    get_selected_data,
    interesting_smells,
    n_runs,
    n_s_data,
    n_s_data_list,
    np,
):
    dpp_f1 = np.zeros(n_s_data)
    dpp_f1_std = np.zeros(n_s_data)
    uniform_f1 = np.zeros(n_s_data)
    uniform_f1_std = np.zeros(n_s_data)

    Y = data[interesting_smells].to_numpy()

    for n_s in range(n_s_data):
        dpp_s = get_selected_data(f'ids_molecule_mcmc_dpp_{n_s_data_list[n_s]}.csv')
        uniform_s = get_selected_data(f'ids_molecule_uniform_{n_s_data_list[n_s]}.csv')

        dpp_s_list = [
            compute_micro_f1(np.array(dpp_s.iloc[:, n_run]), np.array(uniform_s.iloc[:, n_run]), K, Y)
            for n_run in range(n_runs)
        ]
        dpp_f1[n_s] = np.mean(dpp_s_list)
        dpp_f1_std[n_s] = np.std(dpp_s_list) / np.sqrt(n_runs)

        uniform_s_list = [
            compute_micro_f1(np.array(uniform_s.iloc[:, n_run]), np.array(dpp_s.iloc[:, n_run]), K, Y)
            for n_run in range(n_runs)
        ]
        uniform_f1[n_s] = np.mean(uniform_s_list)
        uniform_f1_std[n_s] = np.std(uniform_s_list) / np.sqrt(n_runs)
    return dpp_f1, dpp_f1_std, uniform_f1, uniform_f1_std


@app.cell
def _(dpp_f1, dpp_f1_std, n_s_data_list, plt, uniform_f1, uniform_f1_std):
    plt.errorbar(n_s_data_list, dpp_f1, yerr=dpp_f1_std, capsize=3, marker="o",label="dpp")
    plt.errorbar(n_s_data_list, uniform_f1, yerr=uniform_f1_std, capsize=3, marker="^",label="uniform")
    plt.legend()
    plt.xlabel("# of training data")
    plt.ylabel("f1")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #SURFACTANT
    """)
    return


@app.cell
def _(pd):
    s_df = pd.read_csv('s_pca_data.csv')
    s_df
    return (s_df,)


@app.cell
def _(PandasTools, rdDepictor, s_df, uru):
    s_data = s_df.copy()
    PandasTools.AddMoleculeColumnToFrame(s_data, smilesCol="SMILES")
    s_data.ROMol.apply(rdDepictor.Compute2DCoords)
    s_data["image"] = s_data.ROMol.apply(uru.mol_to_base64_image, target="altair")
    s_data = s_data.drop(columns=["ROMol"])
    return (s_data,)


@app.cell
def _(np, s_data):
    plot_data = s_data.copy()
    plot_data["log_CMC"] = np.log10(plot_data["CMC"])
    return (plot_data,)


@app.cell
def _(alt, mo, plot_data, x, y):
    other_points1 = alt.Chart(plot_data).mark_point(
        color="lightgray",
        opacity=0.8
    ).encode(
        x=alt.X(x.value, scale=alt.Scale(domain=[-1, 1]), title=f"{x.value} ({plot_data['importance_' + x.value][0]:.2f})"),
        y=alt.Y(y.value, scale=alt.Scale(domain=[-1, 1]), title=f"{y.value} ({plot_data['importance_' + y.value][0]:.2f})"),
        tooltip=["image"]
    )

    smell_points1 = alt.Chart(plot_data).mark_point(
        opacity=1.0
    ).encode(
        x=alt.X(x.value, scale=alt.Scale(domain=[-1, 1])),
        y=alt.Y(y.value, scale=alt.Scale(domain=[-1, 1])),
        color=alt.Color("log_CMC:Q", scale=alt.Scale(scheme="viridis")),
        tooltip=["image"]
    )


    chart1 = (other_points1 + smell_points1).properties(
        width=400,
        height=400
    )

    mo.ui.altair_chart(chart1)
    return


@app.cell
def _(get_selected_data):
    s_ids_molecule_uniform = get_selected_data('ss_ids_molecule_uniform_50.csv')
    return (s_ids_molecule_uniform,)


@app.cell
def _(s_data, s_ids_molecule_uniform, selected_molecules, slider):
    selected_molecules(s_data, s_ids_molecule_uniform.iloc[:, slider.value-1])
    return


@app.cell
def _(get_selected_data):
    s_ids_molecule_mcmc = get_selected_data('ss_ids_molecule_mcmc_dpp_50.csv')
    return (s_ids_molecule_mcmc,)


@app.cell
def _(s_data, s_ids_molecule_mcmc, selected_molecules, slider):
    selected_molecules(s_data, s_ids_molecule_mcmc.iloc[:, slider.value-1])
    return


@app.cell
def _(s_ids_molecule_mcmc):
    n_s_runs = s_ids_molecule_mcmc.shape[1]
    return (n_s_runs,)


@app.cell
def _(np, plt, s_data, s_ids_molecule_mcmc, s_ids_molecule_uniform):
    log_cmc_all = np.log10(s_data["CMC"])

    bins = np.linspace(log_cmc_all.min(), log_cmc_all.max(), 30)

    def average_histogram(id_df, s_data, bins):
        densities = []

        for col in id_df.columns:
            ids = id_df[col].astype(int).to_numpy()
            values = np.log10(s_data.loc[ids, "CMC"])

            density, _ = np.histogram(values, bins=bins, density=True)
            densities.append(density)

        densities = np.array(densities)

        mean_density = densities.mean(axis=0)
        std_density = densities.std(axis=0, ddof=1)

        return mean_density, std_density

    mcmc_mean, mcmc_std = average_histogram(s_ids_molecule_mcmc, s_data, bins)
    uniform_mean, uniform_std = average_histogram(s_ids_molecule_uniform, s_data, bins)

    all_density, _ = np.histogram(log_cmc_all, bins=bins, density=True)

    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    bin_width = bins[1] - bins[0]

    plt.plot(bin_centers, all_density, label="all data")
    plt.plot(bin_centers, mcmc_mean, label="mcmc mean")
    plt.plot(bin_centers, uniform_mean, label="uniform mean")

    plt.fill_between(
        bin_centers,
        mcmc_mean - mcmc_std,
        mcmc_mean + mcmc_std,
        alpha=0.2
    )

    plt.fill_between(
        bin_centers,
        uniform_mean - uniform_std,
        uniform_mean + uniform_std,
        alpha=0.2
    )

    plt.ylabel("density")
    plt.xlabel("log[CMC]")
    plt.legend()
    plt.show()
    return all_density, mcmc_mean, uniform_mean


@app.cell
def _(np, plt, s_data, s_ids_molecule_mcmc, s_ids_molecule_uniform):
    plt.hist(np.log10(s_data.loc[s_ids_molecule_mcmc.iloc[:, 5], "CMC"]), label="dpp")
    plt.hist(np.log10(s_data.loc[s_ids_molecule_uniform.iloc[:, 5], "CMC"]), label="uniform")
    plt.legend()
    return


@app.cell
def _(all_density, entropy, mcmc_mean, pd, uniform_mean):
    pd.DataFrame({
    "entropy": [
        entropy(mcmc_mean),
        entropy(uniform_mean),
        entropy(all_density),
    ]}, index=["dpp", "uniform", "all data"])
    return


@app.cell
def _(np, pd):
    Ks = np.array(pd.read_csv('s_Gram_matrix.csv'))
    Ks
    return (Ks,)


@app.cell
def _(KernelRidge, mean_squared_error, np):
    def compute_mse(train_idx, other_idx, K, target): 
        n = len(K)
        all_idx = np.arange(n)
    
        combine_idx = np.unique(np.concatenate([train_idx, other_idx]))
        test_idx = np.setdiff1d(all_idx, combine_idx)
    
        K_train = K[np.ix_(train_idx, train_idx)]
        K_test = K[np.ix_(test_idx, train_idx)]
    
        target_train = target[train_idx]
        target_test = target[test_idx]
    
        model = KernelRidge(alpha=1e-3, kernel="precomputed")
        model.fit(K_train, target_train)
    
        target_pred = model.predict(K_test)
    
        mse = mean_squared_error(target_test, target_pred)
    
        return mse

    return (compute_mse,)


@app.cell
def _():
    n_data_list = [10, 20, 50, 100]# , 60, 80, 100, 150, 200
    n_data = len(n_data_list)
    return n_data, n_data_list


@app.cell
def _(
    Ks,
    compute_mse,
    get_selected_data,
    n_data,
    n_data_list,
    n_s_runs,
    np,
    s_df,
):
    dpp_mse = np.zeros(n_data)
    dpp_mse_std = np.zeros(n_data)
    uniform_mse = np.zeros(n_data)
    uniform_mse_std = np.zeros(n_data)
    for n in range(n_data):
        dpp = get_selected_data(f'ss_ids_molecule_mcmc_dpp_{n_data_list[n]}.csv')
        uniform = get_selected_data(f'ss_ids_molecule_uniform_{n_data_list[n]}.csv')
    
        dpp_list = [
            compute_mse(np.array(dpp.iloc[:, n_run]), np.array(uniform.iloc[:, n_run]), Ks, np.log10(np.array(s_df["CMC"]))) 
            for n_run in range(n_s_runs)
        ]
        dpp_mse[n] = np.mean(dpp_list)
        dpp_mse_std[n] = np.std(dpp_list) / np.sqrt(n_s_runs)
    
        uniform_list = [
            compute_mse(np.array(uniform.iloc[:, n_run]), np.array(dpp.iloc[:, n_run]), Ks, np.log10(np.array(s_df["CMC"]))) 
            for n_run in range(n_s_runs)]
    
        uniform_mse[n] = np.mean(uniform_list)
        uniform_mse_std[n] = np.std(uniform_list) / np.sqrt(n_s_runs)
    return dpp_mse, dpp_mse_std, uniform_mse, uniform_mse_std


@app.cell
def _(dpp_mse, dpp_mse_std, n_data_list, plt, uniform_mse, uniform_mse_std):
    plt.errorbar(n_data_list, dpp_mse, yerr=dpp_mse_std, capsize=3, marker="o",label="dpp")
    plt.errorbar(n_data_list, uniform_mse, yerr=uniform_mse_std, capsize=3, marker="^",label="uniform")
    plt.legend()
    plt.xlabel("# of training data")
    plt.ylabel("MSE of log[CMC]")
    return


if __name__ == "__main__":
    app.run()
