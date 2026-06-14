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
    from scipy.stats import entropy
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
        entropy,
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
    return (X_train,)


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
def _(plt):
    def viz_L(L):
        plt.figure(figsize=(10, 8))
        plt.imshow(L, cmap="coolwarm")
        plt.colorbar()
        plt.tight_layout()
        plt.show()

    return (viz_L,)


@app.cell
def _(X_train, build_L):
    L = build_L(X_train, verbose=True)
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
def _(FiniteDPP, build_L, np, rng):
    def grab_ids_train(X_train, k, sample_method, L=None):
        n_train = X_train.shape[0]

        if sample_method == "uniform":
            ids_train = np.random.choice(
                np.arange(n_train), size=k, replace=False
            )
        if sample_method == "DPP":
            if L is None:
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

        rmses = {}
        for sample_method in ["uniform", "DPP"]:
            ids_train = grab_ids_train(X_train, k, sample_method)
    
            erts = ExtraTreesRegressor()
            erts.fit(X_train[ids_train], y_train[ids_train])
    
            y_test_pred = erts.predict(X_test)
        
            rmses[sample_method] = root_mean_squared_error(y_test, y_test_pred)
        
        return rmses

    return (train_test_run,)


@app.cell
def _(dataset, smiles_to_features, train_test_run):
    train_test_run(dataset, smiles_to_features, 10)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # multiple runs over different train set sizes
    """)
    return


@app.cell
def _(dataset, pd, smiles_to_features, train_test_run):
    ks = [100, 150, 200, 250, 300]
    n_runs = 0

    rows = []
    for k in ks:
        for r in range(n_runs):
            rmses = train_test_run(dataset, smiles_to_features, k)
            for sample_method in ["DPP", "uniform"]:
                rows.append(
                    {
                        "sample_method": sample_method, 
                        "k": k, 
                        "run": r, 
                        "rmse": rmses[sample_method]
                    }
                )
    ml_data = pd.DataFrame(rows)
    ml_data
    return (ml_data,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    for typical performance see [here](https://moleculenet.org/full-results).
    """)
    return


@app.cell
def _(ml_data, sns):
    sns.boxplot(data=ml_data, x="k", y="rmse", hue="sample_method", dodge=True)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # diverse smells

    ## read in data
    """)
    return


@app.cell
def _(pd):
    data_smells = pd.read_csv("https://raw.githubusercontent.com/SimonEnsemble/Leffingwell_Goodscents_combine/refs/heads/main/pyrfume.csv")
    return (data_smells,)


@app.cell
def _(data_smells):
    unique_smells = data_smells.columns[1:].values
    unique_smells
    return (unique_smells,)


@app.cell
def _(np, unique_smells):
    n_total_smells = len(np.unique(unique_smells))
    n_total_smells
    return (n_total_smells,)


@app.cell
def _(data_smells, unique_smells):
    y_smells = data_smells[unique_smells].apply(
        lambda row: row[row == 1].index.tolist(), axis=1
    ).tolist()
    y_smells
    return (y_smells,)


@app.cell
def _(data_smells):
    data_smells
    return


@app.cell
def _(data_smells, featurizer):
    X_smells = featurizer.featurize(data_smells["molecule"].values)
    X_smells
    return (X_smells,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## L-matrix for DPP
    """)
    return


@app.cell
def _(X_smells, build_L):
    L_smells = build_L(X_smells, verbose=True)
    return (L_smells,)


@app.cell
def _(L_smells, viz_L):
    viz_L(L_smells)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## sample randoms sets, look at label dist'n
    """)
    return


@app.cell
def _(grab_ids_train, np):
    def sample_label_dist(
        k, sample_method, X_smells, y_smells, L_smells
    ):
        ids = grab_ids_train(X_smells, k, sample_method, L=L_smells)
    
        smells = np.concatenate([y_smells[id] for id in ids])
    
        unique, counts = np.unique(smells, return_counts=True)
    
        smell_distn = dict(zip(unique, counts))

        return smell_distn

    return (sample_label_dist,)


@app.cell
def _(entropy, np):
    def score_entropy(smell_distn, n_total_smells):
        full_counts = np.zeros(n_total_smells)
        for i, (smell, count) in enumerate(smell_distn.items()):
            full_counts[i] = count
        p = full_counts / full_counts.sum()
        return entropy(p)  # low if smells missing OR if distribution is skewed

    return (score_entropy,)


@app.cell
def _(L_smells, X_smells, sample_label_dist, y_smells):
    smell_distn = sample_label_dist(
        10, "uniform", X_smells, y_smells, L_smells
    )
    smell_distn
    return (smell_distn,)


@app.cell
def _(n_total_smells, score_entropy):
    score_entropy({'fruity': 1}, n_total_smells)
    return


@app.cell
def _(n_total_smells, np, score_entropy, unique_smells):
    max_ent = score_entropy(
        dict(zip(unique_smells, np.ones(n_total_smells))), n_total_smells
    )
    max_ent
    return (max_ent,)


@app.cell
def _(n_total_smells, score_entropy, smell_distn):
    score_entropy(smell_distn, n_total_smells)
    return


@app.cell
def _(pd, sample_label_dist, score_entropy):
    def sample_smell_label_entropy(
        k, X_smells, y_smells, L_smells, n_total_smells, n_runs=100
    ):
        rows = []
        for r in range(n_runs):
            for sample_method in ["DPP", "uniform"]:
                smell_distn = sample_label_dist(
                    k, sample_method, X_smells, y_smells, L_smells
                )
            
                s = score_entropy(smell_distn, n_total_smells)

                rows.append(
                    {
                        'sample method': sample_method, "k": k, "entropy": s
                    }
                )
        return pd.DataFrame(rows)

    return (sample_smell_label_entropy,)


@app.cell
def _(
    L_smells,
    X_smells,
    n_total_smells,
    sample_smell_label_entropy,
    y_smells,
):
    entropy_data = sample_smell_label_entropy(
        250, X_smells, y_smells, L_smells, n_total_smells, n_runs=100
    )
    entropy_data
    return (entropy_data,)


@app.cell
def _(entropy_data, max_ent, plt, sns):
    _p = sns.histplot(entropy_data, x="entropy", hue="sample method", element="step")
    plt.xlim([0, max_ent])
    _p
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
