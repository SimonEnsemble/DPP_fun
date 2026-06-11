import marimo

__generated_with = "0.23.9"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    import importlib.util, urllib.request, sys, os
    import pandas as pd
    from sklearn.ensemble import ExtraTreesClassifier
    from sklearn.metrics import f1_score
    from dppy.finite_dpps import FiniteDPP
    from sklearn.decomposition import PCA
    from sklearn.metrics.pairwise import pairwise_kernels
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.preprocessing import StandardScaler
    import seaborn as sns

    # install:
    # uv venv --python 3.12 dpp
    # uv pip install chemprop
    # url = "https://raw.githubusercontent.com/JacksonBurns/chemeleon/main/chemeleon_fingerprint.py"
    # urllib.request.urlretrieve(url, "chemeleon_fingerprint.py")

    path = os.path.join(os.path.dirname(__file__), "chemeleon_fingerprint.py")
    spec = importlib.util.spec_from_file_location("chemeleon_fingerprint", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    CheMeleonFingerprint = mod.CheMeleonFingerprint
    return (
        CheMeleonFingerprint,
        ExtraTreesClassifier,
        FiniteDPP,
        StandardScaler,
        f1_score,
        np,
        pairwise_kernels,
        pd,
        plt,
        sns,
    )


@app.cell
def _(pd):
    data_train = pd.read_csv("data_from_bees/time_train.csv")
    data_train
    return (data_train,)


@app.cell
def _(pd):
    data_test = pd.read_csv("data_from_bees/time_test.csv")
    data_test
    return (data_test,)


@app.cell
def _(CheMeleonFingerprint):
    fp = CheMeleonFingerprint()
    return (fp,)


@app.cell
def _(data_train, fp):
    X_train = fp(data_train["SMILES"].values)
    X_train
    return (X_train,)


@app.cell
def _(X_train, plt):
    plt.figure(figsize=(10, 8))
    plt.imshow(X_train, cmap="viridis")
    plt.colorbar()
    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(L):
    L.max()
    return


@app.cell
def _(L):
    L.min()
    return


@app.cell
def _(X_scaled):
    X_scaled.mean(axis=0)
    return


@app.cell
def _(X_scaled):
    X_scaled.std(axis=0)
    return


@app.cell
def _(X_scaled):
    X_scaled
    return


@app.cell
def _(n_train):
    n_train
    return


@app.cell
def _(L, n_train):
    assert L.shape[0] == n_train
    return


@app.cell
def _(StandardScaler, X_train):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_train)
    return (X_scaled,)


@app.cell
def _(X_train):
    gamma = 1.0/(X_train.shape[1] * X_train.var())
    gamma
    return (gamma,)


@app.cell
def _(X_train, pairwise_kernels):
    pairwise_kernels(X_train, metric='rbf', gamma=1.0/(X_train.shape[1] * X_train.var()))
    return


@app.cell
def _(X_train, gamma, pairwise_kernels, plt):
    # L = X_scaled.dot(X_scaled.T)
    L = pairwise_kernels(X_train, metric='rbf', gamma=gamma)

    plt.figure(figsize=(10, 8))
    plt.imshow(L, cmap="coolwarm")
    plt.colorbar()
    plt.tight_layout()
    plt.show()
    return (L,)


@app.cell
def _(y_train):
    n_train = len(y_train)
    n_train
    return (n_train,)


@app.cell
def _(data_test, fp):
    X_test = fp(data_test["SMILES"].values)
    X_test
    return (X_test,)


@app.cell
def _(data_train):
    y_train = data_train["label"].values 
    y_train
    return (y_train,)


@app.cell
def _(data_test):
    y_test = data_test["label"].values 
    y_test
    return (y_test,)


@app.cell
def _(y_test):
    n_test = len(y_test)
    n_test
    return


@app.cell
def _(np):
    rng = np.random.RandomState(123)
    return (rng,)


@app.cell
def _(X_train):
    X_train.shape[0]
    return


@app.cell
def _(FiniteDPP, L, X_train, grab_ids_train, k, rng):
    DPP = FiniteDPP('likelihood', **{'L': L})
    DPP.sample_mcmc_k_dpp(
        size=k, 
        random_state=rng, 
        s_init=grab_ids_train(X_train, k, L, "uniform"),
        nb_iter=100
    )
    len(DPP.list_of_samples[-1])
    return


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
def _(X_train):
    X_train
    return


@app.cell
def _(L, X_train, grab_ids_train):
    _ids_train = grab_ids_train(X_train, 10, L, "DPP")
    X_train[_ids_train].shape
    return


@app.cell
def _(ExtraTreesClassifier, f1_score, grab_ids_train, y_train):
    def train_test(X_train, X_test, y_test, k, L, sample_method):
        ids_train = grab_ids_train(X_train, k, L, sample_method)

        erts = ExtraTreesClassifier()
        erts.fit(X_train[ids_train], y_train[ids_train])

        y_test_pred = erts.predict(X_test)
        return f1_score(y_test, y_test_pred)

    return (train_test,)


@app.cell
def _(n_train):
    n_train
    return


@app.cell
def _(L, X_test, X_train, pd, train_test, y_test):
    ks = [25, 50, 75, 100, 125, 150]
    n_runs = 25

    rows = []
    ml_data = pd.DataFrame({"sample method": []})
    for sample_method in ["DPP", "uniform"]:
        for k in ks:
            for r in range(n_runs):
                f1 = train_test(X_train, X_test, y_test, k, L, sample_method)
                rows.append(
                    {"sample_method": sample_method, "k": k, "run": r, "f1": f1}
                )
    ml_data = pd.DataFrame(rows)
    ml_data
    return k, ml_data


@app.cell
def _(ml_data, sns):
    sns.boxplot(data=ml_data, x="k", y="f1", hue="sample_method", dodge=True)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
