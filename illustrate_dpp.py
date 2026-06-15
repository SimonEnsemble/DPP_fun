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
    import matplotlib
    from itertools import combinations
    from sklearn.decomposition import PCA
    from matplotlib.transforms import Bbox, TransformedBbox
    from matplotlib.legend_handler import HandlerBase
    from matplotlib.image import BboxImage
    from sklearn.metrics.pairwise import pairwise_kernels
    from rdkit import Chem
    from rdkit.Chem import Draw
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
        Bbox,
        BboxImage,
        Chem,
        Draw,
        FiniteDPP,
        HandlerBase,
        PCA,
        StandardScaler,
        TransformedBbox,
        combinations,
        dc,
        matplotlib,
        mo,
        np,
        pairwise_kernels,
        pd,
        plt,
        sns,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # specify molecules

    cannabanis componds [here](https://cannabisdatabase.ca/downloads)
    """)
    return


@app.cell
def _():
    molecule_to_smiles = {
        "THC": "CCCCCC1=CC(=C2[C@@H]3C=C(CC[C@H]3C(OC2=C1)(C)C)C)O",
        "CBD": "CCCCCC1=CC(=C(C(=C1)O)[C@@H]2C=C(CC[C@H]2C(=C)C)C)O",
        "CBG": "CCCCCC1=CC(=C(C(=C1)O)C/C=C(\C)/CCC=C(C)C)O",
        "CBN": "Oc2cc(cc1OC(c3c(c12)cc(cc3)C)(C)C)CCCCC",
        "CBC": "CCCCCC1=CC(=C2C=CC(OC2=C1)(C)CCC=C(C)C)O",
        # "terephthalic acid": "O=C(O)c1ccc(C(O)=O)cc1",
        # "malic acid": "O=C(O)CC(O)C(=O)O"
    }
    molecules = list(molecule_to_smiles.keys())
    smiles = list(molecule_to_smiles.values())
    molecule_to_smiles
    return molecules, smiles


@app.cell
def _(molecules):
    nb_pairs = len(molecules) * (len(molecules) - 1)
    nb_pairs
    return


@app.cell
def _(Chem, Draw):
    def draw_molecules(molecules, smiles, row=4):
        mols = [Chem.MolFromSmiles(smiles) for smiles in smiles]
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=row,
            subImgSize=(200, 200),
            legends=molecules
        )
        return img

    return (draw_molecules,)


@app.cell
def _(draw_molecules, molecules, smiles):
    draw_molecules(molecules, smiles)
    return


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
def _(featurizer, smiles):
    X = featurizer.featurize(smiles)
    X
    return (X,)


@app.cell
def _(StandardScaler, X):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    return (X_scaled,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # PCA
    """)
    return


@app.cell
def _(PCA, X_scaled):
    pca = PCA(n_components=2)
    X_reduced = pca.fit_transform(X_scaled)
    X_reduced
    return X_reduced, pca


@app.cell
def _(pca):
    pca.explained_variance_ratio_
    return


@app.cell
def _(
    Bbox,
    BboxImage,
    HandlerBase,
    TransformedBbox,
    X_reduced,
    draw_molecules,
    matplotlib,
    molecules,
    np,
    pca,
    plt,
    smiles,
):
    class HandlerLineImage(HandlerBase):
        def __init__(self, image_data, space=15, offset=10):
            self.space = space
            self.offset = offset
            self.image_data = image_data
            super().__init__()

        def create_artists(
            self,
            legend,
            orig_handle,
            xdescent,
            ydescent,
            width,
            height,
            fontsize,
            trans,
        ):
            # plot marker, adjust position of marker
            l = matplotlib.lines.Line2D(
                [
                    xdescent + self.offset * 2 - self.space,
                    xdescent + self.offset * 2 - self.space,
                ],
                [ydescent + height / 6.0, ydescent + height / 6.0],
                marker="o",
                markersize=40,
            )
            l.update_from(orig_handle)
            l.set_clip_on(False)
            l.set_transform(trans)

            # adjust position of molecule
            bb = Bbox.from_bounds(
                xdescent + (width + self.space) / 3.0 + self.offset / 2,
                ydescent - self.space,
                height * self.image_data.shape[1] / self.image_data.shape[0]
                - self.space / 2,
                height - self.space / 2,
            )
            tbb = TransformedBbox(bb, trans)
            image = BboxImage(tbb)
            image.set_data(self.image_data)

            self.update_prop(image, orig_handle, legend)
            return [l, image]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.axhline(0, color="black", linewidth=1)
    ax.axvline(0, color="black", linewidth=1)

    lines = [[] for _ in range(len(smiles))]
    imgs = [[] for _ in range(len(smiles))]
    markers = ["o", "s", "^", "v", "<", ">", "D", "d", "p", "h", "X", "+", "*", "8", "|"]
    for _i in range(len(molecules)):
        smi = smiles[_i]
        mol = molecules[_i]
        imgs[_i] = np.asarray(draw_molecules([mol], [smi], row=1))
    
        (line,) = ax.plot(
            X_reduced[_i, 0],
            X_reduced[_i, 1],
            linestyle="None",
            marker=markers[_i],
            markersize=8,
            # markerfacecolor=smiles_to_color[smi],
            # markeredgecolor=smiles_to_color[smi],
        )
        lines[_i] = line

    handler_map = {line: HandlerLineImage(img) for line, img in zip(lines, imgs)}
    ax.legend(
        lines,
        [""] * len(lines),
        handler_map=handler_map,
        ncol=2,
        loc="upper left",
        bbox_to_anchor=(0.9, 1.4),
        columnspacing=0.1,
        handlelength=0.8,
        labelspacing=0.1,
        fontsize=80,
        handleheight=1.2,
    )
    ax.set_xlabel("PC1 (" + str(np.round(pca.explained_variance_ratio_[0], 2)) + ")")
    ax.set_ylabel("PC2 (" + str(np.round(pca.explained_variance_ratio_[1], 2)) + ")")
    # plt.savefig("pca.png", format="pdf", bbox_inches="tight")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # L-matrix for DPP
    """)
    return


@app.cell
def _(X):
    gamma = 1.0 / (X.shape[1] * X.var())
    return


@app.cell
def _(X_scaled, pairwise_kernels):
    L = pairwise_kernels(
        X_scaled, metric='rbf'#, gamma=gamma
    )
    L
    return (L,)


@app.cell
def _(L, molecules, pd, sns):
    sns.heatmap(
        pd.DataFrame(L, index=molecules, columns=molecules),
        annot=True, fmt=".2f"
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # the DPP
    """)
    return


@app.cell
def _(np):
    def grab_minor(L, ids):
        return np.array([[L[i, j] for j in ids] for i in ids])

    return (grab_minor,)


@app.cell
def _(L, grab_minor):
    grab_minor(L, [1, 3])
    return


@app.cell
def _(combinations, grab_minor, np):
    def norm_factor(k, L):
        n = L.shape[0]
    
        f = 0.0
        for ids in combinations(range(n), k):
            L_minor = grab_minor(L, ids)
            f += np.linalg.det(L_minor)
        return f

    return (norm_factor,)


@app.cell
def _(grab_minor, norm_factor, np):
    def prob(ids, L):
        k = len(ids)
        L_minor = grab_minor(L, ids)
        return np.linalg.det(L_minor) / norm_factor(k, L)

    return (prob,)


@app.cell
def _(L, norm_factor):
    norm_factor(3, L)
    return


@app.cell
def _(L, prob):
    prob([1, 2], L)
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
    def sample_molecules(L, k):
        DPP = FiniteDPP('likelihood', **{'L': L})
        DPP.sample_mcmc_k_dpp(
            size=k, 
            random_state=rng,
            nb_iter=1000
        )
        return np.sort(DPP.list_of_samples[-1][-1])

    return (sample_molecules,)


@app.cell
def _(L, sample_molecules):
    sample_molecules(L, 2)
    return


@app.cell
def _(L, np, sample_molecules):
    n_sim = 1000
    np.sum(
        [sample_molecules(L, 2) for s in range(n_sim)] == np.array([1, 2])
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
