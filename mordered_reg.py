import marimo

__generated_with = "0.23.9"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    from sklearn.ensemble import ExtraTreesRegressor
    from sklearn.metrics import root_mean_squared_error
    from dppy.finite_dpps import FiniteDPP
    from sklearn.decomposition import PCA
    from sklearn.metrics.pairwise import pairwise_kernels
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.preprocessing import StandardScaler
    import seaborn as sns
    from rdkit import Chem
    # use mordred community https://github.com/jacksonburns/mordred-community
    # uv pip install "mordredcommunity[full]==2.0.6"
    from mordred import Calculator, descriptors
    from deepchem.molnet import load_freesolv
    import deepchem as dc

    return (
        ExtraTreesRegressor,
        FiniteDPP,
        StandardScaler,
        dc,
        load_freesolv,
        mo,
        np,
        pairwise_kernels,
        pd,
        plt,
        root_mean_squared_error,
        sns,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # read in and split data
    """)
    return


@app.cell
def _(load_freesolv):
    tasks, datasets, transformers = load_freesolv(splitter=None)
    dataset = datasets[0]
    return (dataset,)


@app.cell
def _(dataset, dc):
    # splitter = dc.splits.ButinaSplitter()
    splitter = dc.splits.RandomSplitter()

    train_dataset, test_dataset = splitter.train_test_split(dataset)
    return test_dataset, train_dataset


@app.cell
def _(test_dataset, train_dataset):
    print("# train: ", len(train_dataset))
    print("# test: ", len(test_dataset))
    return


@app.cell
def _(train_dataset):
    train_dataset
    return


@app.cell
def _(test_dataset, train_dataset):
    y_train = train_dataset.y
    y_test = test_dataset.y
    return y_test, y_train


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # descriptors
    """)
    return


@app.cell
def _(dc):
    featurizer = dc.feat.MordredDescriptors(ignore_3D=True)
    return (featurizer,)


@app.cell
def _():
    # old, using mordred-community
    # calc = Calculator(descriptors, ignore_3D=True)
    # len(calc.descriptors)
    # train_mols = [
    #     Chem.MolFromSmiles(smi) for smi in train_dataset.ids
    # ]
    # train_descriptors = calc.pandas(train_mols)
    return


@app.cell
def _(featurizer, train_dataset):
    X_train = featurizer.featurize(train_dataset.ids)
    return (X_train,)


@app.cell
def _(featurizer, test_dataset):
    X_test = featurizer.featurize(test_dataset.ids)
    return (X_test,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # L-matrix for DPP
    """)
    return


@app.cell
def _(StandardScaler, X_train, pairwise_kernels):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_train)
    gamma = 1.0 / (X_scaled.shape[1] * X_scaled.var())
    L = pairwise_kernels(
        X_scaled, metric='rbf', gamma=gamma
    )
    gamma
    return (L,)


@app.cell
def _(L, plt):
    plt.figure(figsize=(10, 8))
    plt.imshow(L, cmap="coolwarm")
    plt.colorbar()
    plt.tight_layout()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # sampler of training data
    """)
    return


@app.cell
def _(np):
    rng = np.random.RandomState(123)
    return (rng,)


@app.cell
def _(FiniteDPP, np, rng):
    def grab_ids_train(X_train, k, L, sample_method):
        n_train = X_train.shape[0]

        if sample_method == "uniform":
            ids_train = np.random.choice(
                np.arange(n_train), size=k, replace=False
            )
        if sample_method == "DPP":
            DPP = FiniteDPP('likelihood', **{'L': L})
            DPP.sample_mcmc_k_dpp(
                size=k, 
                random_state=rng, 
                s_init=grab_ids_train(X_train, k, L, "uniform"),
                nb_iter=1000
            )
            ids_train = DPP.list_of_samples[-1][-1]

        return ids_train

    return (grab_ids_train,)


@app.cell
def _(L, X_train, grab_ids_train):
    _ids_train = grab_ids_train(X_train, 10, L, "DPP")
    X_train[_ids_train].shape
    _ids_train
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # do train/test run
    """)
    return


@app.cell
def _(ExtraTreesRegressor, grab_ids_train, root_mean_squared_error, y_train):
    def train_test(X_train, X_test, y_test, k, L, sample_method):
        ids_train = grab_ids_train(X_train, k, L, sample_method)

        erts = ExtraTreesRegressor()
        erts.fit(X_train[ids_train], y_train[ids_train])

        y_test_pred = erts.predict(X_test)
        return root_mean_squared_error(y_test, y_test_pred)

    return (train_test,)


@app.cell
def _(X_test):
    X_test
    return


@app.cell
def _(L, X_test, X_train, train_test, y_test):
    train_test(X_train, X_test, y_test, 10, L, "DPP")
    return


@app.cell
def _(y_train):
    len(y_train)
    return


@app.cell
def _(L, X_test, X_train, pd, train_test, y_test):
    ks = [100, 200, 300]
    n_runs = 10

    rows = []
    ml_data = pd.DataFrame({"sample method": []})
    for sample_method in ["DPP", "uniform"]:
        for k in ks:
            for r in range(n_runs):
                rmse = train_test(X_train, X_test, y_test, k, L, sample_method)
                rows.append(
                    {"sample_method": sample_method, "k": k, "run": r, "rmse": rmse}
                )
    ml_data = pd.DataFrame(rows)
    ml_data
    return (ml_data,)


@app.cell
def _(ml_data, sns):
    sns.boxplot(data=ml_data, x="k", y="rmse", hue="sample_method", dodge=True)
    return


if __name__ == "__main__":
    app.run()
