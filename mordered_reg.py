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
    # read in data
    """)
    return


@app.cell
def _(load_freesolv):
    tasks, datasets, transformers = load_freesolv(splitter=None)
    dataset = datasets[0]
    return (dataset,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # featurize the molecules
    """)
    return


@app.cell
def _(dc):
    featurizer = dc.feat.MordredDescriptors(ignore_3D=True)
    return (featurizer,)


@app.cell
def _(dataset, featurizer):
    _features = featurizer.featurize(dataset.ids)
    smiles_to_features = dict(
        zip(dataset.ids, _features)
    )
    smiles_to_features
    return (smiles_to_features,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # test/train split
    """)
    return


@app.cell
def _(dc, np):
    def train_test_split(dataset, smiles_to_features, verbose=True):
        splitter = dc.splits.RandomSplitter()
        train_dataset, test_dataset = splitter.train_test_split(dataset)
    
        if verbose:
            print("# train: ", len(train_dataset))
            print("# test: ", len(test_dataset))

        X_train = np.array([smiles_to_features[sm] for sm in train_dataset.ids])
        X_test  = np.array([smiles_to_features[sm] for sm in  test_dataset.ids])

        y_train = train_dataset.y.ravel()
        y_test = test_dataset.y.ravel()
    
        return X_train, y_train, X_test, y_test

    return (train_test_split,)


@app.cell
def _(dataset, smiles_to_features, train_test_split):
    X_train, y_train, X_test, y_test = train_test_split(dataset, smiles_to_features)
    return X_test, X_train, y_test, y_train


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # L-matrix for DPP
    """)
    return


@app.cell
def _(StandardScaler, pairwise_kernels):
    def build_L(X_train, verbose=False):
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_train)
        gamma = 1.0 / (X_scaled.shape[1] * X_scaled.var())
        if verbose:
            print("gamma: ", gamma)
        
        L = pairwise_kernels(
            X_scaled, metric='rbf', gamma=gamma
        )
        return L

    return (build_L,)


@app.cell
def _(X_train, build_L, plt):
    L = build_L(X_train, verbose=True)

    plt.figure(figsize=(10, 8))
    plt.imshow(L, cmap="coolwarm")
    plt.colorbar()
    plt.tight_layout()
    plt.show()
    return (L,)


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
def _(FiniteDPP, build_L, np, rng):
    def grab_ids_train(X_train, k, sample_method):
        n_train = X_train.shape[0]

        if sample_method == "uniform":
            ids_train = np.random.choice(
                np.arange(n_train), size=k, replace=False
            )
        if sample_method == "DPP":
            L = build_L(X_train)
        
            DPP = FiniteDPP('likelihood', **{'L': L})
            DPP.sample_mcmc_k_dpp(
                size=k, 
                random_state=rng, 
                s_init=grab_ids_train(X_train, k, "uniform"),
                nb_iter=1000
            )
            ids_train = DPP.list_of_samples[-1][-1]

        return ids_train

    return (grab_ids_train,)


@app.cell
def _(X_train, grab_ids_train):
    _ids_train = grab_ids_train(X_train, 10, "uniform")
    X_train[_ids_train]#.shape
    # _ids_train
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # do train/test run
    """)
    return


@app.cell
def _(
    ExtraTreesRegressor,
    grab_ids_train,
    root_mean_squared_error,
    train_test_split,
):
    def train_test_run(
        dataset, smiles_to_features, k, verbose=False
    ):
        X_train, y_train, X_test, y_test = train_test_split(
            dataset, smiles_to_features, verbose=verbose
        )

        res = {}
        for sample_method in ["uniform", "DPP"]:
            ids_train = grab_ids_train(X_train, k, sample_method)
    
            erts = ExtraTreesRegressor()
            erts.fit(X_train[ids_train], y_train[ids_train])
    
            y_test_pred = erts.predict(X_test)
        
            res[sample_method] = root_mean_squared_error(y_test, y_test_pred)
        
        return res

    return (train_test_run,)


@app.cell
def _(dataset, smiles_to_features, train_test_run):
    train_test_run(dataset, smiles_to_features, 10)
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
def _(L, X_test, X_train, ks, n_runs, pd, train_test, y_test):
    # ks = [100, 200, 300]
    # n_runs = 10

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


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # diverse smells
    """)
    return


if __name__ == "__main__":
    app.run()
