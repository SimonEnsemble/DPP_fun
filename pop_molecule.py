import marimo

__generated_with = "0.23.8"
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
    from sklearn.metrics import root_mean_squared_error

    from sklearn.svm import SVC
    from sklearn.multiclass import OneVsRestClassifier
    from sklearn.metrics import f1_score
    from sklearn.dummy import DummyClassifier

    return (
        Chem,
        DummyClassifier,
        KernelRidge,
        OneVsRestClassifier,
        PandasTools,
        SVC,
        alt,
        f1_score,
        mo,
        np,
        pd,
        plt,
        rdDepictor,
        root_mean_squared_error,
        uru,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #read in the data
    """)
    return


@app.cell
def _(mo):
    dropdown = mo.ui.dropdown(
        options=["smells", "bees", "surfactants"], value="bees", label="dataset"
    )
    dropdown
    return (dropdown,)


@app.cell
def _(dropdown):
    datapath = f"data_from_{dropdown.value}"
    return (datapath,)


@app.cell
def _(datapath, dropdown, np, pd):
    df = pd.read_csv(datapath + '/PCA_perfect.csv')
    if dropdown.value == "surfactants":
        df["log_CMC"] = np.log10(df["CMC"])
    df
    return (df,)


@app.cell
def _(PandasTools, df, rdDepictor, uru):
    data = df.copy()
    PandasTools.AddMoleculeColumnToFrame(data, smilesCol="SMILES")#molecule
    data.ROMol.apply(rdDepictor.Compute2DCoords)
    data["image"] = data.ROMol.apply(uru.mol_to_base64_image, target="altair")
    data = data.drop(columns=["ROMol"])
    return (data,)


@app.cell
def _(data):
    data
    return


@app.cell
def _(datapath, pd):
    PCA_importance = pd.read_csv(datapath + "/PCA_importance.csv")
    PCA_importance
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #viz PCA
    """)
    return


@app.cell
def _(data, mo):
    x = mo.ui.dropdown(options=data.columns.drop("image"), value="x1")
    x
    return (x,)


@app.cell
def _(data, mo):
    y = mo.ui.dropdown(options=data.columns.drop("image"), value="x2")
    y
    return (y,)


@app.cell
def _():
    pca_range = [-1, 1]
    return (pca_range,)


@app.cell
def _(data, dropdown, mo, np):
    if dropdown.value == "smells":
        smell_list = data.columns.drop(["image", "SMILES", "x1", "x2"])
        interesting_smells = ["vanilla", "mint", "wine", "chamomile", "currant", "fish", "popcorn"]
        smell = mo.ui.dropdown(options=np.sort(smell_list), value="wine") #  np.sort(smell_list)
    return interesting_smells, smell, smell_list


@app.cell
def _(dropdown, mo, smell):
    mo.vstack([
        smell
    ]) if dropdown.value == "smells" else mo.md("")
    return


@app.cell
def _(alt, data, dropdown, mo, smell, x, y):
    if dropdown.value == "smells":
        other_points = alt.Chart(data[data[smell.value] == 0]).mark_point(
            color="lightgray",
            opacity=0.8
        ).encode(
            x=x.value,#alt.X(x.value, scale=alt.Scale(domain=[-1, 1]), title=f"{x.value} ({PCA_importance['imp'][0]:.2f})"),
            y=y.value,#alt.Y(y.value, scale=alt.Scale(domain=[-1, 1]), title=f"{y.value} ({PCA_importance['imp'][1]:.2f})"),
            tooltip=["image", smell.value]
        )

        smell_points = alt.Chart(data[data[smell.value] == 1]).mark_point(
            color="red",
            opacity=1.0
        ).encode(
            x=x.value,#alt.X(x.value, scale=alt.Scale(domain=[-1, 1])),
            y=y.value,#alt.Y(y.value, scale=alt.Scale(domain=[-1, 1])),
            tooltip=["image", smell.value]
        )

        chart = (other_points + smell_points).properties(
            width=400,
            height=400
        )

        out = mo.ui.altair_chart(chart)

    elif dropdown.value == "surfactants":
        other_points = alt.Chart(data).mark_point(
            color="lightgray",
            opacity=0.8
        ).encode(
            x=x.value,#alt.X(x.value, scale=alt.Scale(domain=[-1, 1]), title=f"{x.value} ({PCA_importance['imp'][0]:.2f})"),
            y=y.value,#alt.Y(y.value, scale=alt.Scale(domain=[-1, 1]), title=f"{y.value} ({PCA_importance['imp'][1]:.2f})"),
            tooltip=["image"]
        )

        smell_points = alt.Chart(data).mark_point(
            opacity=1.0
        ).encode(
            x=x.value,#alt.X(x.value, scale=alt.Scale(domain=[-1, 1])),
            y=y.value,#alt.Y(y.value, scale=alt.Scale(domain=[-1, 1])),
            color=alt.Color("log_CMC:Q", scale=alt.Scale(scheme="viridis")),
            tooltip=["image"]
        )


        chart = (other_points + smell_points).properties(
            width=400,
            height=400
        )

        out = mo.ui.altair_chart(chart)

    elif dropdown.value == "bees":
        non_toxic_points = alt.Chart(data[data["toxic"] == 0]).mark_point(
            color="lightgray",
            opacity=0.8
        ).encode(
            x=x.value,#=alt.X(x.value, scale=alt.Scale(domain=[-1, 1]), title=f"{x.value} ({PCA_importance['imp'][0]:.2f})"),
            y=y.value,#=alt.Y(y.value, scale=alt.Scale(domain=[-1, 1]), title=f"{y.value} ({PCA_importance['imp'][1]:.2f})"),
            tooltip=["image", "toxic"]
        )

        toxic_points = alt.Chart(data[data["toxic"] == 1]).mark_point(
            color="red",
            opacity=1.0
        ).encode(
            x=x.value,#=alt.X(x.value, scale=alt.Scale(domain=[-1, 1])),
            y=y.value,#=alt.Y(y.value, scale=alt.Scale(domain=[-1, 1])),
            tooltip=["image", "toxic"]
        )

        bee_chart = (non_toxic_points + toxic_points).properties(
            width=400,
            height=400
        )

        out = mo.ui.altair_chart(bee_chart)

    out
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # get selected data
    """)
    return


@app.cell
def _(pd):
    def get_selected_data(file):
        ids_molecule = pd.read_csv(file)
        return ids_molecule - 1

    return (get_selected_data,)


@app.cell
def _(alt, datapath, get_selected_data, mo, np, pca_range, x, y):
    def selected_molecules(data, selected_ids, pca_range=pca_range):
        plot_data = data.copy()
        plot_data["selected"] = 0

        train_idx = np.array(get_selected_data(datapath + '/train_data.csv')).ravel()
        train_data = plot_data.loc[train_idx]
        plot_data.loc[train_idx[selected_ids], "selected"] = 1

        plot_x = x.value#alt.X(x.value, scale=alt.Scale(domain=pca_range), title=f"{x.value} ({PCA_importance['imp'][0]:.2f})")
        plot_y = y.value#alt.Y(y.value, scale=alt.Scale(domain=pca_range), title=f"{y.value} ({PCA_importance['imp'][1]:.2f})")

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

        return mo.ui.altair_chart(other_molecules + selected_molecules) #+ smell_points)

    return (selected_molecules,)


@app.cell
def _(data, datapath, get_selected_data, np):
    train_data = data.loc[np.array(get_selected_data(datapath + '/train_data.csv')).ravel()]
    test_data = data.loc[np.array(get_selected_data(datapath + '/test_data.csv')).ravel()]
    return


@app.cell
def _(datapath, get_selected_data):
    ids_molecule_mcmc_dpp = get_selected_data(datapath + '/ids_dpp_10.csv')
    return (ids_molecule_mcmc_dpp,)


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
def _(data, ids_molecule_mcmc_dpp, np, selected_molecules, slider):
    selected_molecules(data, np.array(ids_molecule_mcmc_dpp.iloc[:, slider.value-1]).ravel())
    return


@app.cell
def _(datapath, get_selected_data):
    ids_molecule_uniform = get_selected_data(datapath + '/ids_uniform_10.csv')
    return (ids_molecule_uniform,)


@app.cell
def _(data, ids_molecule_uniform, np, selected_molecules, slider):
    selected_molecules(data, np.array(ids_molecule_uniform.iloc[:, slider.value-1]).ravel())
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # estimate uniform and DPP sampling
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## regression
    """)
    return


@app.cell
def _(KernelRidge, datapath, get_selected_data, np, root_mean_squared_error):
    def compute_mse(selected_ids, K, target, alpha=0.1): 
        train_idx = np.array(get_selected_data(datapath + '/train_data.csv')).ravel()
        train_ids = train_idx[selected_ids]

        test_idx = np.array(get_selected_data(datapath + '/test_data.csv')).ravel()

        K_train = K[np.ix_(train_ids, train_ids)]
        K_test = K[np.ix_(test_idx, train_ids)]

        target_train = target[train_ids]
        target_test = target[test_idx]

        model = KernelRidge(alpha=alpha, kernel="precomputed")
        model.fit(K_train, target_train)

        target_pred = model.predict(K_test)

        mse = root_mean_squared_error(target_test, target_pred)

        return mse
        # return model.score(K_test, target_test)
    return (compute_mse,)


@app.cell
def _(Ks, compute_mse, data, dropdown, ids_molecule_mcmc_dpp, n_runs, np):
    if dropdown.value == "surfactants":
        alpha_list = [0.0001, 0.001, 0.01, 0.1, 1.0, 10.0]
        [
            [compute_mse(np.array(ids_molecule_mcmc_dpp.iloc[:, n_run]), Ks, data["log_CMC"].to_numpy(), alpha) for alpha in alpha_list] 
            for n_run in range(n_runs)
        ]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## classification
    """)
    return


@app.cell
def _(OneVsRestClassifier, SVC, datapath, f1_score, get_selected_data, np):
    def compute_f1(selected_ids, K, target, C=10, method="macro"):
        train_idx = np.array(get_selected_data(datapath + '/train_data.csv')).ravel()
        train_ids = train_idx[selected_ids]

        test_idx = np.array(get_selected_data(datapath + '/test_data.csv')).ravel()

        K_train = K[np.ix_(train_ids, train_ids)]
        K_test = K[np.ix_(test_idx, train_ids)]

        target_train = target[train_ids]
        target_test = target[test_idx]

        model = OneVsRestClassifier(
            SVC(C=C, kernel="precomputed", class_weight="balanced")
        )
        model.fit(K_train, target_train)
        target_pred = model.predict(K_test)

        return f1_score(target_test, target_pred, average=method, zero_division=0)
        # return model.score(K_test, target_test)
    return (compute_f1,)


@app.cell
def _(Ks, compute_f1, data, dropdown, ids_molecule_mcmc_dpp, n_runs, np):
    if dropdown.value == "bees":
        C_grid = [0.001, 0.01, 0.1, 1.0, 10.0, 100]
        [
            [compute_f1(np.array(ids_molecule_mcmc_dpp.iloc[:, n_run]), Ks, data["toxic"].to_numpy(), C, method="macro") for C in C_grid] 
            for n_run in range(n_runs)
        ]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### kernel training
    """)
    return


@app.cell
def _(datapath, np, pd):
    Ks = np.array(pd.read_csv(datapath + '/Gram_matrix_perfect.csv'))
    Ks
    return (Ks,)


@app.cell
def _(
    DummyClassifier,
    Ks,
    compute_f1,
    compute_mse,
    data,
    datapath,
    df,
    dropdown,
    f1_score,
    get_selected_data,
    interesting_smells,
    n_runs,
    np,
    plt,
    root_mean_squared_error,
):
    n_data_list = [10, 20, 50, 80, 100]
    n_data = len(n_data_list)
    train_idx = np.array(get_selected_data(datapath + '/train_data.csv')).ravel()
    test_idx = np.array(get_selected_data(datapath + '/test_data.csv')).ravel()

    if dropdown.value == "smells":
        dpp_f1 = np.zeros(n_data)
        dpp_f1_std = np.zeros(n_data)
        uniform_f1 = np.zeros(n_data)
        uniform_f1_std = np.zeros(n_data)

        Y = data[interesting_smells].to_numpy()

        for n in range(n_data):
            dpp_s = get_selected_data(datapath + f'/ids_dpp_{n_data_list[n]}.csv')
            uniform_s = get_selected_data(datapath + f'/ids_dpp_{n_data_list[n]}.csv')

            dpp_s_list = [
                compute_f1(np.array(dpp_s.iloc[:, n_run]), Ks, Y)
                for n_run in range(n_runs)
            ]
            dpp_f1[n] = np.mean(dpp_s_list)
            dpp_f1_std[n] = np.std(dpp_s_list) / np.sqrt(n_runs)

            uniform_s_list = [
                compute_f1(np.array(uniform_s.iloc[:, n_run]), Ks, Y)
                for n_run in range(n_runs)
            ]
            uniform_f1[n] = np.mean(uniform_s_list)
            uniform_f1_std[n] = np.std(uniform_s_list) / np.sqrt(n_runs)

        baseline = DummyClassifier(strategy="most_frequent")
        K_train = Ks[np.ix_(train_idx, train_idx)]
        K_test = Ks[np.ix_(test_idx, train_idx)]
        target_train = Y[train_idx]
        target_test = Y[test_idx]
        baseline.fit(K_train, target_train)
        target_pred = baseline.predict(K_test)
        f1 = f1_score(target_test, target_pred, average="macro", zero_division=0)

        plt.errorbar(n_data_list, dpp_f1, yerr=dpp_f1_std, capsize=3, marker="o",label="dpp")
        plt.errorbar(n_data_list, uniform_f1, yerr=uniform_f1_std, capsize=3, marker="^",label="uniform")
        plt.plot(n_data_list, np.ones(n_data)*f1, label="baseline")
        plt.legend()
        plt.xlabel("# of training data")
        plt.ylabel("f1")
        plt.show()

    elif dropdown.value == "surfactants":
        dpp_mse = np.zeros(n_data)
        dpp_mse_std = np.zeros(n_data)
        uniform_mse = np.zeros(n_data)
        uniform_mse_std = np.zeros(n_data)
        for n in range(n_data):
            dpp = get_selected_data(datapath + f'/ids_dpp_{n_data_list[n]}.csv')
            uniform = get_selected_data(datapath + f'/ids_uniform_{n_data_list[n]}.csv')

            dpp_list = [
                compute_mse(np.array(dpp.iloc[:, n_run]), Ks, np.array(df["log_CMC"])) 
                for n_run in range(n_runs)
            ]
            dpp_mse[n] = np.mean(dpp_list)
            dpp_mse_std[n] = np.std(dpp_list) / np.sqrt(n_runs)

            uniform_list = [
                compute_mse(np.array(uniform.iloc[:, n_run]), Ks, np.array(df["log_CMC"])) 
                for n_run in range(n_runs)]

            uniform_mse[n] = np.mean(uniform_list)
            uniform_mse_std[n] = np.std(uniform_list) / np.sqrt(n_runs)

        baseline = root_mean_squared_error(data.loc[test_idx, "log_CMC"], np.ones(len(test_idx))*np.mean(data.loc[train_idx, "log_CMC"]))

        plt.errorbar(n_data_list, dpp_mse, yerr=dpp_mse_std, capsize=3, marker="o",label="dpp")
        plt.errorbar(n_data_list, uniform_mse, yerr=uniform_mse_std, capsize=3, marker="^",label="uniform")
        plt.plot(n_data_list, np.ones(n_data)*baseline, label="baseline")
        plt.legend()
        plt.ylim(ymin=0)
        plt.xlabel("# of training data")
        plt.ylabel("MSE of log[CMC]")
        plt.show()

    elif dropdown.value == "bees":
        bee_dpp_f1 = np.zeros(n_data)
        bee_dpp_f1_std = np.zeros(n_data)
        bee_uniform_f1 = np.zeros(n_data)
        bee_uniform_f1_std = np.zeros(n_data)

        bee_Y = data["toxic"].to_numpy()

        for n in range(n_data):
            bee_ids_molecule_mcmc = get_selected_data(datapath + f'/ids_dpp_{n_data_list[n]}.csv')
            bee_ids_molecule_uniform = get_selected_data(datapath + f'/ids_uniform_{n_data_list[n]}.csv')

            bee_dpp_s_list = [
                compute_f1(
                    np.array(bee_ids_molecule_mcmc.iloc[:, n_run]).ravel(), Ks, bee_Y, method="macro"
                )
                for n_run in range(n_runs)
            ]
            bee_dpp_f1[n] = np.mean(bee_dpp_s_list)
            bee_dpp_f1_std[n] = np.std(bee_dpp_s_list) / np.sqrt(n_runs)

            bee_uniform_s_list = [
                compute_f1(
                    np.array(bee_ids_molecule_uniform.iloc[:, n_run]).ravel(), Ks, bee_Y, method="macro"
                )
                for n_run in range(n_runs)
            ]
            bee_uniform_f1[n] = np.mean(bee_uniform_s_list)
            bee_uniform_f1_std[n] = np.std(bee_uniform_s_list) / np.sqrt(n_runs)

        baseline = DummyClassifier(strategy="most_frequent")
        K_train = Ks[np.ix_(train_idx, train_idx)]
        K_test = Ks[np.ix_(test_idx, train_idx)]
        target_train = bee_Y[train_idx]
        target_test = bee_Y[test_idx]
        baseline.fit(K_train, target_train)
        target_pred = baseline.predict(K_test)
        f1 = f1_score(target_test, target_pred, average="macro", zero_division=0)

        plt.errorbar(n_data_list, bee_dpp_f1, yerr=bee_dpp_f1_std, capsize=3, marker="o",label="dpp")
        plt.errorbar(n_data_list, bee_uniform_f1, yerr=bee_uniform_f1_std, capsize=3, marker="^",label="uniform")
        plt.plot(n_data_list, np.ones(n_data)*f1, label="baseline")
        plt.legend()
        plt.ylim(0, 1)
        plt.xlabel("# of training data")
        plt.ylabel("f1")
        plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## estimation
    """)
    return


@app.cell
def _(np):
    def entropy(p):
        p = np.array(p)
        p = p[p > 0]  
        return -np.sum(p * np.log(p))

    return (entropy,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ##viz similarity
    """)
    return


@app.cell
def _(np):
    def compute_pairwise_kernel_and_jaccard(K, smell_df, upper_only=True):

        K = np.asarray(K)
        Y = smell_df.to_numpy(dtype=np.int8)

        n = K.shape[0]

        # intersection: number of labels both molecules have
        intersection = Y @ Y.T  # shape (n, n)

        # row sums = number of positive labels per molecule
        row_sum = Y.sum(axis=1, keepdims=True)  # shape (n, 1)

        # union = labels in i + labels in j - intersection
        union = row_sum + row_sum.T - intersection

        # Jaccard similarity
        jaccard = np.divide(
            intersection,
            union,
            out=np.zeros_like(intersection, dtype=float),
            where=(union > 0)
        )

        if upper_only:
            pair_idx = np.triu_indices(n, k=1)
        else:
            pair_idx = np.indices((n, n)).reshape(2, -1)

        kernel_sims = K[pair_idx]
        jaccard_sims = jaccard[pair_idx]

        return kernel_sims, jaccard_sims

    return (compute_pairwise_kernel_and_jaccard,)


@app.cell
def _(plt):
    def plot_similarity(kernel_sims, jaccard_sims, ylabel, gridsize=40):
        plt.figure(figsize=(6, 5))

        hb = plt.hexbin(
            kernel_sims,
            jaccard_sims,
            gridsize=gridsize,
            mincnt=1,
            bins="log",
            cmap="viridis"
        )

        plt.xlabel("kernel similarity")
        plt.ylabel(ylabel)

        cbar = plt.colorbar(hb)
        cbar.set_label("log(count)")

        plt.tight_layout()
        plt.show()

    return (plot_similarity,)


@app.cell
def _(np):
    def compute_pairwise_kernel_and_logcmc_diff(K, cmc_df, target_col="log_CMC", upper_only=True):

        K = np.asarray(K)
        y = cmc_df[target_col].to_numpy(dtype=float)

        n = K.shape[0]

        if upper_only:
            pair_idx = np.triu_indices(n, k=1)
        else:
            pair_idx = np.indices((n, n)).reshape(2, -1)

        i_idx, j_idx = pair_idx
        kernel_sims = K[i_idx, j_idx]
        logcmc_diffs = np.abs(y[i_idx] - y[j_idx])

        return kernel_sims, logcmc_diffs

    return (compute_pairwise_kernel_and_logcmc_diff,)


@app.cell
def _(
    Ks,
    compute_pairwise_kernel_and_jaccard,
    compute_pairwise_kernel_and_logcmc_diff,
    data,
    dropdown,
    plot_similarity,
    smell_list,
):
    if dropdown.value == "smells":
        kernel_sims, jaccard_sims = compute_pairwise_kernel_and_jaccard(Ks, data[smell_list], upper_only=True)
        plot_similarity(kernel_sims, jaccard_sims, "property similarity")
    elif dropdown.value == "surfactants":
        s_kernel_sims, s_logcmc_diffs = compute_pairwise_kernel_and_logcmc_diff(Ks, data)
        plot_similarity(s_kernel_sims, s_logcmc_diffs, "log[CMC] differences")
    elif dropdown.value == "bees":
        bee_kernel_sims, bee_jaccard_sims = compute_pairwise_kernel_and_jaccard(Ks, data[["toxic"]], upper_only=True)
        plot_similarity(bee_kernel_sims, bee_jaccard_sims, "property similarity")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #still need to work
    """)
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
    interesting_smells,
    n_runs,
    num_smells,
    pd,
):
    n_smell = pd.DataFrame({
        "all": (data[interesting_smells].sum(axis=0) > 0).sum(),
        "mcmc_DPP": num_smells(ids_molecule_mcmc_dpp, n_runs, interesting_smells),
        "uniform": num_smells(ids_molecule_uniform, n_runs, interesting_smells),
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
def _(
    data,
    ids_molecule_mcmc_dpp,
    ids_molecule_uniform,
    interesting_smells,
    n_runs,
    pd,
    plt,
    proportion_smell,
):
    dpp_smell = proportion_smell(data, ids_molecule_mcmc_dpp, interesting_smells, n_runs)
    uniform_smell = proportion_smell(data, ids_molecule_uniform, interesting_smells, n_runs)

    all_smell = data[interesting_smells].sum()
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
    return


@app.cell
def _(
    data,
    entropy,
    ids_molecule_mcmc_dpp,
    ids_molecule_uniform,
    interesting_smells,
    n_runs,
    pd,
    proportion_smell,
):
    all_data_smell = data[interesting_smells].sum()
    all_data_smell = all_data_smell / all_data_smell.sum()
    pd.DataFrame({
    "entropy": [
        entropy(proportion_smell(data, ids_molecule_mcmc_dpp, interesting_smells, n_runs)),
        entropy(proportion_smell(data, ids_molecule_uniform, interesting_smells, n_runs)),
        entropy(all_data_smell),
    ]}, index=["dpp", "uniform", "all data"])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #SURFACTANT
    """)
    return


@app.cell
def _(data, ids_molecule_mcmc, ids_molecule_uniform, np, plt, s_data):
    bins = np.linspace(np.log10(s_data["CMC"]).min(), np.log10(s_data["CMC"]).max(), 30)

    def average_histogram(ids_molecule, s_data, bins):
        densities = []

        for col in ids_molecule.columns:
            ids = ids_molecule[col].astype(int).to_numpy()
            values = np.log10(s_data.loc[ids, "CMC"])

            density, _ = np.histogram(values, bins=bins, density=True)
            densities.append(density)

        densities = np.array(densities)

        mean_density = densities.mean(axis=0)
        std_density = densities.std(axis=0, ddof=1)

        return mean_density, std_density

    mcmc_mean, mcmc_std = average_histogram(ids_molecule_mcmc, data, bins)
    uniform_mean, uniform_std = average_histogram(ids_molecule_uniform, data, bins)

    all_density, _ = np.histogram(np.log10(data["CMC"]), bins=bins, density=True)

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
    row=8
    plt.hist(np.log10(s_data.loc[s_ids_molecule_mcmc.iloc[:, row], "CMC"]), label="dpp", density=True)
    plt.hist(np.log10(s_data.loc[s_ids_molecule_uniform.iloc[:, row], "CMC"]), label="uniform", density=True)
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
def _():
    # tune alpha
    return


@app.cell
def _():
    # use the mean of the train set as the prediction
    return


@app.cell
def _(Chem):
    def atom_labels_with_chirality(smiles):
        mol = Chem.MolFromSmiles(smiles)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

        labels = []
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if atom.HasProp("_CIPCode"):
                labels.append(f"{symbol}_{atom.GetProp('_CIPCode')}")
            else:
                labels.append(symbol)
        return labels

    smiles = "C[C@H](O)Cl"
    labels = atom_labels_with_chirality(smiles)

    print(labels)
    return (smiles,)


@app.cell
def _(Chem, smiles):

    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D


    def draw_molecule_with_rs_labels(smiles, output_file="molecule_rs.svg"):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")

        # Assign R/S stereochemistry
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

        # Make 2D coordinates for drawing
        Chem.rdDepictor.Compute2DCoords(mol)

        # Create SVG drawer
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 350)
        opts = drawer.drawOptions()

        # Replace atom labels with C_R / C_S when available
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()

            if atom.HasProp("_CIPCode"):
                cip = atom.GetProp("_CIPCode")  # "R" or "S"
                opts.atomLabels[idx] = f"{symbol}_{cip}"
            else:
                opts.atomLabels[idx] = symbol

        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()

        svg = drawer.GetDrawingText()

        with open(output_file, "w") as f:
            f.write(svg)

        return svg

    draw_molecule_with_rs_labels(smiles, "molecule_ss.svg")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
