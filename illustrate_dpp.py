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
    from rdkit.Chem import PandasTools
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
        PandasTools,
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
def _(PandasTools):
    cannabinoids = PandasTools.LoadSDF(
        "compounds.sdf", smilesName="SMILES", molColName="Molecule"
    )
    cannabinoid_ids = (
        ["CDB00000" + str(i) for i in range(1, 10)] + 
        ["CDB0000" + str(i) for i in range(10, 44)] # + 
        # ["CDB0063" + str(i) for i in range(46, 53)] # terpenes
    )
    cannabinoids = cannabinoids[
        cannabinoids["Cannabis Database ID"].isin(cannabinoid_ids)
    ]

    # cory's names
    cannabinoids["GENERIC_NAME"] = cannabinoids["GENERIC_NAME"].replace(
        {
            # see Cannabis compound database
            "Cannabidiol": "CBD",
            "Delta-9-tetrahydrocannabinol": "THC",
            "Cannabigerol": "CBG",
            "Cannabigerovarinic acid": "CBGVA",
            "Cannabidivarin": "CBDV",
            "Cannabichromene": "CBC",
            "Cannabichromevarinic acid": "CBCVA",
            "Cannabidiolic acid": "CBDA",
            "Cannabigerolic acid monomethylether": "CBGAM",
            "Cannabidivarinic acid": "CBDVA",
            "Cannabigerolic acid": "CBGA",
            "Delta-8-tetrahydrocannabinol": "delta-8-THC",
            "Delta-9-tetrahydrocannabinol-C4": "THC-C4",
            "Cannabicyclolic acid": "CBLA",
            "Cannabicyclol": "CBL",
            "Delta-9-cis-tetrahydrocannabinol": "cis-THC",
            "Cannabielsoin": "CBE",
            "Cannabinolic acid": "CBNA",
            "Cannabinol": "CBN",
            "Cannabinol methylether": "CBNM", 
            "Cannabinol-C4": "CBN-C4", 
            "Cannabivarin": "CBV"
        }
    )

    # drop indistinguishable
    drop_indistinguishable = True
    if drop_indistinguishable:
        cannabinoids = cannabinoids[cannabinoids["GENERIC_NAME"] != "cis-THC"]

    cannabinoids = cannabinoids.reset_index()
    cannabinoids
    return cannabinoids, drop_indistinguishable


@app.cell
def _(cannabinoids):
    assert cannabinoids["GENERIC_NAME"].nunique() == cannabinoids.shape[0]
    return


@app.cell
def _(cannabinoids):
    cannabinoids["SMILES"].nunique()
    return


@app.cell
def _(Chem, Draw):
    def draw_molecules(cannabinoids, row=4):
        mols = [Chem.MolFromSmiles(smiles) for smiles in cannabinoids["SMILES"]]
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=row,
            subImgSize=(200, 200),
            legends=cannabinoids["GENERIC_NAME"].values
        )
        return img

    return (draw_molecules,)


@app.cell
def _(cannabinoids):
    cannabinoid_subset = [
        "CBD", "THC", "CBG", "CBGVA", "CBDV", "CBC", "CBCVA",
        "CBDA", "CBGAM", "CBDVA", "CBGA", "delta-8-THC", "THC-C4", "CBL", "CBLA", 
        "CBE", "CBNA", "CBV", "CBN-C4"
    ]
    for _c in cannabinoid_subset:
        assert _c in cannabinoids["GENERIC_NAME"].values, f"{_c} not found"
    cannabinoid_subset = ["THC", "CBGA", "CBGAM", "CBLA", "CBDV", "CBD", "CBGVA"]
    return (cannabinoid_subset,)


@app.cell
def _(cannabinoid_subset, np):
    assert len(np.unique(cannabinoid_subset)) == len(cannabinoid_subset)
    return


@app.cell
def _(cannabinoids, draw_molecules):
    draw_molecules(cannabinoids)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # featurize the molecules
    """)
    return


@app.cell
def _(dc):
    # featurizer = dc.feat.MordredDescriptors(ignore_3D=True)
    featurizer = dc.feat.RDKitDescriptors()
    return (featurizer,)


@app.cell
def _(cannabinoids, featurizer):
    X = featurizer.featurize(cannabinoids["SMILES"].values)
    X
    return (X,)


@app.cell
def _(X):
    X[0].shape # number of features
    return


@app.cell
def _(combinations, np):
    def look_for_duplicates(X):
        n = X.shape[0]
        same_rep = []
        for i, j in combinations(range(n), 2):
            if np.allclose(X[i], X[j]):
                same_rep.append((i, j))
        return same_rep

    return (look_for_duplicates,)


@app.cell
def _(X, look_for_duplicates):
    dups = look_for_duplicates(X)

    # if drop_indistinguishable:
    #     assert dups == []
    
    dups
    return


@app.cell
def _(cannabinoids, draw_molecules, drop_indistinguishable):
    if not drop_indistinguishable:
        draw_molecules(cannabinoids.iloc[[17, 37]])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # dim reduction
    """)
    return


@app.cell
def _(StandardScaler, X):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    return (X_scaled,)


@app.cell
def _(np):
    def diffusion_map(L):
        n = L.shape[0]
    
        # make row-stochastic matrix
        L = L / L.sum(axis=1, keepdims=True)

        evals, evecs = np.linalg.eig(L)
        evals = evals.real
        evecs = evecs.real

        # sort by descending eigenvalue
        order = np.argsort(evals)[::-1]
        evals = evals[order]
        evecs = evecs[:, order]

        assert np.isclose(evals[0], 1.0), f"largest eigenvalue is {evals[0]}, not 1"

        X_reduced = np.zeros((n, 2))
        X_reduced[:, 0] = evecs[:, 1] * evals[1]  # 2nd eigenvector
        X_reduced[:, 1] = evecs[:, 2] * evals[2] # 3rd eigenvector

        return X_reduced

    return (diffusion_map,)


@app.cell
def _(L, diffusion_map):
    X_reduced_dm = diffusion_map(L)
    X_reduced_dm
    return


@app.cell
def _(X_reduced):
    X_reduced.shape
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
    cannabinoid_subset,
    cannabinoids,
    draw_molecules,
    matplotlib,
    np,
    pca,
    plt,
    sns,
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

    def viz(X_reduced, cannabinoids, molecule_subset):
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.axhline(0, color="black", linewidth=1)
        ax.axvline(0, color="black", linewidth=1)



        lines = [[] for _ in range(len(molecule_subset))]
        imgs = [[] for _ in range(len(molecule_subset))]
        markers = [
            "s", "o", "^", "v", "<", ">", "D", "d", "p", "h", "X", "+", "*", "8", "|"
        ]
        markers = markers + markers
        colors = sns.color_palette("husl", len(molecule_subset))

        id = -1
        for i, row in cannabinoids.iterrows():
            if not row["GENERIC_NAME"] in molecule_subset:
                ax.plot(
                    X_reduced[i, 0],
                    X_reduced[i, 1],
                    linestyle="None",
                    marker="o",
                    markersize=8,
                    markerfacecolor="black",
                    markeredgecolor="black", zorder=0
                )
            else:
                id += 1
                imgs[id] = np.asarray(draw_molecules(cannabinoids.iloc[[i]], row=1))
                (line,) = ax.plot(
                    X_reduced[i, 0],
                    X_reduced[i, 1],
                    linestyle="None",
                    marker=markers[id],
                    markersize=10,
                    markerfacecolor=colors[id],
                    markeredgecolor=colors[id]
                )
                lines[id] = line
    
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
        ax.set_xlabel("PC1 (" + str(np.round(pca.explained_variance_ratio_[0]*100, 0)) + "%)")
        ax.set_ylabel("PC2 (" + str(np.round(pca.explained_variance_ratio_[1]*100, 0)) + "%)")
        ax.set_aspect('equal', 'box')
        # plt.savefig("pca.png", format="pdf", bbox_inches="tight")
        plt.show()

    viz(X_reduced, cannabinoids, cannabinoid_subset)
    return (viz,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # L-matrix for DPP
    """)
    return


@app.cell
def _(X_scaled, pairwise_kernels):
    L = pairwise_kernels(
        X_scaled, metric='rbf'
    )
    L
    return (L,)


@app.cell
def _(np):
    def grab_L_minor(L, cannabinoids, cannabinoid_subset):
        ids = [
            cannabinoids[cannabinoids["GENERIC_NAME"] == m].index[0] 
            for m in cannabinoid_subset
        ]

        return L[np.ix_(ids, ids)]

    return (grab_L_minor,)


@app.cell
def _(grab_L_minor, pd, plt, sns):
    def viz_L_subset(L, cannabinoids, cannabinoid_subset):
        L_sub = grab_L_minor(L, cannabinoids, cannabinoid_subset)
    
        df = pd.DataFrame(
            L_sub, index=cannabinoid_subset, columns=cannabinoid_subset
        )
    
        sns.heatmap(df, annot=True, fmt=".2f", vmin=0.0, vmax=1.0, cmap="viridis")
        plt.xlabel("cannabinoid")
        plt.ylabel("cannabinoid")
        plt.show()

    return (viz_L_subset,)


@app.cell
def _(L, cannabinoid_subset, cannabinoids, viz_L_subset):
    viz_L_subset(L, cannabinoids, cannabinoid_subset)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # the DPP
    """)
    return


@app.cell
def _(L, cannabinoid_subset, cannabinoids, grab_L_minor):
    L_sub = grab_L_minor(L, cannabinoids, cannabinoid_subset)
    L_sub
    return (L_sub,)


@app.cell
def _(combinations, grab_L_minor, np):
    def norm_factor(k, L, cannabinoids):
        n = L.shape[0]

        f = 0.0
        for cannabinoid_subset in combinations(cannabinoids["GENERIC_NAME"].values, k):
            L_minor = grab_L_minor(L, cannabinoids, cannabinoid_subset)
            f += np.linalg.det(L_minor)
        return f

    return (norm_factor,)


@app.cell
def _(L, cannabinoids, norm_factor):
    p_norm = norm_factor(3, L, cannabinoids)
    return


@app.cell
def _(grab_L_minor, np):
    def prob_prop(L, cannabinoids, cannabinoid_subset):
        L_minor = grab_L_minor(L, cannabinoids, cannabinoid_subset)
        return np.linalg.det(L_minor)

    return (prob_prop,)


@app.cell
def _(L, cannabinoids, prob_prop):
    prob_prop(L, cannabinoids, ["THC", "CBD", "CBLA"])
    return


@app.cell
def _(L, cannabinoids, prob_prop):
    prob_prop(L, cannabinoids, ["CBD", "CBDV", "CBGVA"])
    return


@app.cell
def _(L, cannabinoids, prob_prop):
    prob_prop(L, cannabinoids, ["CBGA", "CBGAM", "CBGVA"])
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
    def sample_molecules(L, k, cannabinoid_list):
        DPP = FiniteDPP('likelihood', **{'L': L})
        DPP.sample_mcmc_k_dpp(
            size=k,
            random_state=rng,
            nb_iter=1000
        )
        return cannabinoid_list[np.sort(DPP.list_of_samples[-1][-1])]

    return (sample_molecules,)


@app.cell
def _(L_sub, cannabinoid_subset, np, sample_molecules):
    sample_molecules(L_sub, 3, np.array(cannabinoid_subset))
    return


@app.cell
def _(L, cannabinoids, sample_molecules):
    dpp_cannabinoid_subset = sample_molecules(L, 3, cannabinoids["GENERIC_NAME"].values)
    dpp_cannabinoid_subset
    return (dpp_cannabinoid_subset,)


@app.cell
def _(X_reduced, cannabinoids, dpp_cannabinoid_subset, viz):
    viz(X_reduced, cannabinoids, dpp_cannabinoid_subset)
    return


@app.cell
def _(L, cannabinoids, np, sample_molecules):
    n_sim = 1000
    np.sum(
        [
            set(sample_molecules(L, 3, cannabinoids["GENERIC_NAME"].values))
            for s in range(n_sim)
        ] == set(["THC", "CBD", "CBLA"])
    ) / n_sim
    return


if __name__ == "__main__":
    app.run()
