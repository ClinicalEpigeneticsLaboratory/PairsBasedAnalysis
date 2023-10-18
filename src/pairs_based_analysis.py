import subprocess
from os.path import join
from os.path import exists
from os import makedirs
import typing as t

from datetime import datetime
from statsmodels.stats.multitest import multipletests

from venn import pseudovenn
from venn import venn

from pathlib import Path
import pandas as pd
import numpy as np
from glob import glob
from tqdm import tqdm
import scipy.stats as sts
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from matplotlib.patches import Patch

np.seed = 101


class PairsBasedAnalysis:
    def __init__(
        self,
        project_dir_path: str,
        mask: str,
        mynorm: str,
        sample_sheet: str,
        manifest: str,
        non_co_methylation_threshold: float,
        alpha: float = 0.05,
    ):
        self.paths = {
            "base": project_dir_path,
            "models": join(project_dir_path, "models"),
            "pairs": join(project_dir_path, "pairs"),
            "plots": join(project_dir_path, "plots"),
            "files": join(project_dir_path, "files"),
        }

        self.mask = self.__load(mask)
        self.mynorm = self.__load(mynorm)
        self.sample_sheet = self.__load(sample_sheet)
        assert (
            "Sample_Group" in self.sample_sheet.columns
        ), "No Sample_Group column in sample sheet."

        self.mynorm = self.mynorm[self.sample_sheet.index]
        assert not self.mynorm.empty, "Mynorm is empty."

        print(f"Mynorm shape: {self.mynorm.shape}")
        print(f"Sample sheet shape: {self.sample_sheet.shape}")
        self.delta_frame = None

        self.manifest = self.__load(
            manifest,
            [
                "CHR",
                "MAPINFO",
                "UCSC_RefGene_Name",
                "Relation_to_UCSC_CpG_Island",
                "UCSC_RefGene_Group",
            ],
        )
        self.manifest["Relation_to_UCSC_CpG_Island"] = self.manifest[
            "Relation_to_UCSC_CpG_Island"
        ].fillna("OpenSea")

        self.non_co_methylation_threshold = non_co_methylation_threshold
        self.alpha = alpha
        self.metadata = {
            "Project": project_dir_path,
            "Start": datetime.now(),
            "Mynorm": mynorm,
            "Sample_sheet": sample_sheet,
            "Manifest": manifest,
            "Non co-methylation threshold": self.non_co_methylation_threshold,
            "Mask": mask,
            "Mynorm shape": self.mynorm.shape,
            "Samples": self.sample_sheet.Sample_Group.value_counts(),
        }
        self.models = {}
        self.ea_manifest = {}
        self.type_to_type_comp = None
        self.comparison_frame = None
        self.pairs_per_marker_type = {}
        self.non_co_methylated_pairs_per_group = {}

    @staticmethod
    def __load(path: str, columns_to_select: t.Optional[list] = False) -> pd.DataFrame:
        if path.endswith(".csv"):
            if columns_to_select:
                return pd.read_csv(path, index_col=0)[columns_to_select]

            return pd.read_csv(path, index_col=0)

        if path.endswith(".parquet"):
            if columns_to_select:
                return pd.read_parquet(path, columns=columns_to_select)

            return pd.read_parquet(path)

        raise Exception("Wrong format of input file.")

    def __build_output_dir(self) -> None:
        for path in self.paths.values():
            makedirs(path, exist_ok=True)

    def __update_mask(self) -> None:
        cpgs = set(self.mynorm.index)
        pairs = []
        self.metadata["Initial_number_of_pairs"] = self.mask.Pair.nunique()

        for pair in tqdm(self.mask.Pair):
            base, successive = pair.split("-")

            if base in cpgs and successive in cpgs:
                pairs.append(pair)

        self.metadata["Number_of_pairs_in_current_dataset"] = len(pairs)
        self.mask = self.mask[self.mask.Pair.isin(pairs)]

    def __fit_model_per_pair(self, samples: list, group: str) -> pd.DataFrame:
        models = []

        for pair in tqdm(self.mask.Pair, desc=group):
            base, successive = pair.split("-")

            base_met_level = self.mynorm.loc[base, samples]
            successive_met_level = self.mynorm.loc[successive, samples]

            y = pd.concat([base_met_level, successive_met_level]).values
            x = [*[0] * len(samples), *[1] * len(samples)]

            model = sts.linregress(x=x, y=y)
            slope, intercept, rvalue, pvalue = (
                model.slope,
                model.intercept,
                model.rvalue,
                model.pvalue,
            )

            models.append(
                {
                    "Pair": pair,
                    "Intercept": intercept,
                    "Slope": slope,
                    "R": rvalue,
                    "p-value": pvalue,
                    "Base CpG mean methylation level": base_met_level.mean(),
                    "Successive CpG mean methylation level": successive_met_level.mean(),
                }
            )

        models = pd.DataFrame(models).set_index("Pair")
        _, models["FDR"], _, _ = multipletests(models["p-value"], method="fdr_bh")
        models.to_csv(join(self.paths["models"], f"{group}.csv"))

        return models

    def _add_status(self, _pvalue, _slope):
        if _pvalue <= self.alpha and abs(_slope) >= self.non_co_methylation_threshold:
            return "discordantly methylated"
        elif _pvalue > self.alpha or abs(_slope) < 0.05:
            return "co-methylated"
        else:
            return "between cutoffs"

    def fit_model_per_group(self) -> None:
        existing_files = glob(join(self.paths["models"], "*.csv"))
        existing_files = [file for file in existing_files if not "stats.csv" in file]

        if existing_files:
            print(f"In dir: {self.paths['models']} found {existing_files}.")
            for file in existing_files:
                stats_frame = pd.read_csv(file, index_col=0)
                group = Path(file).name.replace(".csv", "")
                stats_frame["Status"] = [
                    self._add_status(p, s)
                    for p, s in zip(stats_frame.FDR, stats_frame.Slope)
                ]
                self.models[group] = stats_frame

        else:
            print(
                f"Dir: {self.paths['models']} is empty, fitting models per each pair in each group of samples."
            )
            for group, samples in self.sample_sheet.groupby(
                "Sample_Group"
            ).groups.items():
                samples = samples.tolist()
                stats_frame = self.__fit_model_per_pair(samples, group)
                stats_frame["Status"] = [
                    self._add_status(p, s)
                    for p, s in zip(stats_frame.FDR, stats_frame.Slope)
                ]
                self.models[group] = stats_frame

    def build_comparison_frame(self) -> None:
        comparison_frame = []

        for group, model in self.models.items():
            model = model[["Slope", "Intercept", "FDR", "Status"]]
            model = model.add_prefix(f"{group}_")
            comparison_frame.append(model)

        comparison_frame = pd.concat(comparison_frame, axis=1)
        self.comparison_frame = comparison_frame

    def mark_types_of_markers(self) -> None:
        comparison_frame = self.comparison_frame.copy()
        status_columns = [
            column for column in comparison_frame.columns if "Status" in column
        ]

        for pair, record in tqdm(comparison_frame.iterrows()):
            record = record[status_columns]
            record = record.values

            if all(record == "discordantly methylated"):
                status = "common discordantly methylated"

            elif all(record == "co-methylated"):
                status = "common co-methylated"

            else:
                status = "between cutoffs"

            comparison_frame.loc[pair, "marker_status"] = status

        self.comparison_frame = comparison_frame
        self.metadata["Status_per_pair"] = comparison_frame.marker_status.value_counts()

    def mark_co_methylated_DMRS(self):
        dmr = []
        frame = self.comparison_frame
        frame = frame[frame["marker_status"] == "common co-methylated"]
        frame = frame[[name for name in frame.columns if "Intercept" in name]]

        for pair, record in tqdm(frame.iterrows()):
            for c1_intercept in record.values:
                for c2_intercept in record.values:
                    diff = abs(c2_intercept - c1_intercept)
                    if diff >= self.non_co_methylation_threshold:
                        dmr.append(pair)

        self.comparison_frame.loc[dmr, "marker_status"] = "co-methylated cell specific"

    def create_delta_norm(self):
        base_cpgs = [pair.split("-")[0] for pair in self.mask.Pair]
        successive_cpgs = [pair.split("-")[1] for pair in self.mask.Pair]

        mynorm_delta = (
            self.mynorm.loc[base_cpgs, :].values
            - self.mynorm.loc[successive_cpgs, :].values
        )
        new_index = [
            f"{base}-{nearby}" for base, nearby in zip(base_cpgs, successive_cpgs)
        ]

        delta_frame = pd.DataFrame(
            mynorm_delta, columns=self.mynorm.columns, index=new_index
        )
        self.delta_frame = delta_frame

    def calculate_type_to_type_stats(self):
        if exists(join(self.paths["models"], f"stats.csv")):
            self.type_to_type_comp = pd.read_csv(
                join(self.paths["models"], f"stats.csv"), index_col=0
            )
        else:
            delta_frame = pd.concat((self.delta_frame.T, self.sample_sheet), axis=1)
            stats = []

            for pair in tqdm(self.comparison_frame.index):
                if (
                    not "discordantly methylated"
                    in self.comparison_frame.loc[pair].values
                ):
                    continue

                temp_frame = delta_frame[[pair, "Sample_Group"]]
                to_iterate_over = temp_frame.groupby("Sample_Group").groups.items()

                for group_a, samples_a in to_iterate_over:
                    for group_b, samples_b in to_iterate_over:
                        if group_a != group_b:
                            a = temp_frame.loc[samples_a, pair]
                            b = temp_frame.loc[samples_b, pair]

                            mean_a, mean_b = a.mean(), b.mean()
                            delta = mean_b - mean_a

                            _, pval = sts.kruskal(a, b)
                            stats.append(
                                {
                                    "Pair": pair,
                                    "A": group_a,
                                    "B": group_b,
                                    "Delta": delta,
                                    "p-value": pval,
                                }
                            )

            stats = pd.DataFrame(stats).set_index("Pair")
            _, stats["FDR"], _, _ = multipletests(stats["p-value"], method="fdr_bh")
            stats.to_csv(join(self.paths["models"], f"stats.csv"))

            self.type_to_type_comp = stats

    def mark_markers(self):
        frame = self.type_to_type_comp
        frame = frame[
            (frame["FDR"] <= self.alpha)
            & (frame.Delta.abs() > self.non_co_methylation_threshold)
        ]

        self.type_to_type_comp = frame

        frame = frame.loc[~frame.index.duplicated(keep="first")]
        self.comparison_frame.loc[
            frame.index, "marker_status"
        ] = "discordantly methylated cell specific"
        self.pairs_per_marker_type["discordantly methylated cell specific"] = frame

    def mark_specific_markers(self):
        frame = self.type_to_type_comp
        specific_markers = []

        for pair in frame.index.unique():
            temp = frame.loc[pair]
            temp = temp[(temp.Delta.abs() >= 0.3) & (temp.FDR <= 0.05)]
            if temp.shape[0] == self.sample_sheet.Sample_Group.nunique():
                specific_markers.append(pair)

        self.comparison_frame.loc[specific_markers, "Specific"] = True

    def export_frames_per_type(self) -> None:
        for mtype, pairs in self.comparison_frame.groupby(
            "marker_status"
        ).groups.items():
            pairs_per_mtype = self.comparison_frame.loc[pairs]

            self.pairs_per_marker_type[mtype] = pairs_per_mtype
            pairs_per_mtype.to_csv(join(self.paths["pairs"], f"{mtype}.csv"))

    def plots(
        self,
        collections: t.Dict[str, pd.DataFrame],
        show: bool = False,
        limit: int = 1000,
    ) -> None:
        mynorm = self.mynorm.T
        border = self.non_co_methylation_threshold / 2

        for mtype, pairs in collections.items():
            to_plot = pairs.index[:limit]

            for pair in tqdm(to_plot, desc=f"generating 2D plots for {mtype}"):
                base, successive = pair.split("-")
                fig = px.scatter(
                    mynorm,
                    x=base,
                    y=successive,
                    color=self.sample_sheet.Sample_Group,
                    title=pair,
                )

                fig.update_layout(
                    width=1200,
                    height=1000,
                    font=dict(size=32),
                    xaxis=dict(range=(0, 1), title="Base CpG"),
                    yaxis=dict(range=(0, 1), title="Consecutive CpG"),
                    legend=dict(title=""),
                )

                fig.add_trace(
                    go.Scatter(
                        x=[0, 1],
                        y=[-border, 1 - border],
                        mode="lines",
                        name=f"non co-methylation border +/- {border}",
                        showlegend=False,
                        line=dict(color="red", width=2, dash="dash"),
                    )
                )

                fig.add_trace(
                    go.Scatter(
                        x=[0, 1],
                        y=[0, 1],
                        mode="lines",
                        name="co-methylation",
                        line=dict(color="red", width=2, dash="dot"),
                    )
                )

                fig.add_trace(
                    go.Scatter(
                        x=[0, 1],
                        y=[border, 1 + border],
                        mode="lines",
                        name=f"non co-methylation border +/- {border}",
                        line=dict(color="red", width=2, dash="dash"),
                    )
                )
                fig.update_layout(showlegend=False)
                if show:
                    fig.show()

                out_path = join(self.paths["plots"], mtype)
                makedirs(out_path, exist_ok=True)
                fig.write_image(join(out_path, f"{pair}.png"), engine="kaleido")

    def cluster_map(
        self,
        collections: t.Dict[str, pd.DataFrame],
        features_to_skip: t.Optional[list] = None,
        show: bool = True,
        font_size: float = 1.2,
        dpi: int = 450,
        fig_size: tuple = (13, 13),
        cmap: str = "husl",
        method: str = "ward",
        metric: str = "seuclidean",
        bbox_to_anchor: tuple = (1.05, 1.1),
        legend_loc: str = "upper right",
        legend_font_size: int = 14,
    ) -> None:

        for mtype, pairs in collections.items():
            base_cpgs = [pair.split("-")[0] for pair in pairs.index]
            successive_cpgs = [pair.split("-")[1] for pair in pairs.index]

            mynorm_delta = (
                self.mynorm.loc[base_cpgs, :].values
                - self.mynorm.loc[successive_cpgs, :].values
            )
            new_index = [
                f"{base}-{nearby}" for base, nearby in zip(base_cpgs, successive_cpgs)
            ]

            delta_frame = pd.DataFrame(
                mynorm_delta, columns=self.mynorm.columns, index=new_index
            )

            # Add color bar
            features = self.sample_sheet.copy()
            if features_to_skip:
                features = features.drop(features_to_skip, axis=1)

            global_lut = {}

            for column in features.columns:
                single_feature = features[column]
                pal = sns.color_palette(cmap, single_feature.nunique())
                lut = dict(zip(map(str, single_feature.unique()), pal))
                global_lut = global_lut | lut

                features[column] = features[column].map(lut)

            features.index.name = ""

            sns.set(font_scale=font_size)
            plt.figure(figsize=fig_size)

            fig = sns.clustermap(
                delta_frame,
                method=method,
                cmap="vlag",
                metric=metric,
                col_colors=features,
                tree_kws=dict(linewidths=1.2),
                yticklabels=False,
                xticklabels=False,
            )

            handles = [Patch(facecolor=global_lut[name]) for name in global_lut.keys()]
            plt.legend(
                handles,
                global_lut,
                title="",
                bbox_to_anchor=bbox_to_anchor,
                bbox_transform=plt.gcf().transFigure,
                loc=legend_loc,
                fontsize=legend_font_size,
                edgecolor="white",
                facecolor="white",
            )

            fig.savefig(
                join(self.paths["plots"], f"{mtype}.png"),
                dpi=dpi,
                format="png",
                transparent=True,
            )

            if show:
                plt.show()

    def intersection_map(
        self,
        alpha: float = 0.05,
        dpi=450,
        show: bool = True,
        method="ward",
        metric="seuclidean",
        fig_size: tuple = (13, 13),
        font_size=1.5,
    ) -> None:
        pairwise_matrix = pd.DataFrame(
            columns=self.models.keys(), index=self.models.keys()
        )

        for group, models in self.models.items():
            discordant = models[
                (models["FDR"] <= alpha)
                & (models.Slope.abs() > self.non_co_methylation_threshold)
            ].index

            discordant = set(discordant)
            self.metadata[f"Number_of_discordant_pairs_in_{group}"] = len(discordant)
            self.non_co_methylated_pairs_per_group[group] = discordant

            for paired_group, paired_model in self.models.items():

                paired_discordant = paired_model[
                    (paired_model["FDR"] <= alpha)
                    & (paired_model.Slope.abs() >= self.non_co_methylation_threshold)
                ].index
                paired_discordant = set(paired_discordant)

                intersection = len(set.intersection(discordant, paired_discordant))
                union = len(set.union(discordant, paired_discordant))

                j_index = intersection / union
                pairwise_matrix.loc[group, paired_group] = j_index

        pairwise_matrix = pairwise_matrix.astype(float)
        sns.set(font_scale=font_size)
        plt.figure(figsize=fig_size)

        fig = sns.clustermap(
            pairwise_matrix,
            cmap="vlag",
            annot=True,
            method=method,
            metric=metric,
            linewidth=1.1,
            tree_kws=dict(linewidths=1.2),
            vmin=0,
            vmax=1,
        )
        if show:
            plt.show()

        fig.savefig(
            join(self.paths["plots"], "intersection_map.png"), dpi=dpi, format="png"
        )

    def venn(
        self,
        show: bool = True,
        figsize: tuple = (20, 20),
        fontsize: int = 16,
    ) -> t.Optional[str]:

        # REAL DISCORDANT ONLY
        if len(self.non_co_methylated_pairs_per_group.items()) == 6:
            fig = pseudovenn(
                self.non_co_methylated_pairs_per_group,
                figsize=figsize,
                fontsize=fontsize,
            )
        elif 6 > len(self.non_co_methylated_pairs_per_group) > 1:
            fig = venn(
                self.non_co_methylated_pairs_per_group,
                figsize=figsize,
                fontsize=fontsize,
            )
        else:
            return "Not Applicable in case of < 2 OR > 6 samples groups."

        leg = fig.get_legend()
        leg.get_frame().set_edgecolor("b")
        leg.get_frame().set_linewidth(0.0)

        fig = fig.get_figure()

        if show:
            plt.show()

        fig.savefig(join(self.paths["plots"], "venn.png"))

    def export_bed(self, collections: t.Dict[str, pd.DataFrame]) -> None:
        for mtype, pairs in collections.items():

            base_cpgs = [pair.split("-")[0] for pair in pairs.index]
            successive_cpgs = [pair.split("-")[1] for pair in pairs.index]

            base_cpgs = self.manifest.loc[base_cpgs, ["CHR", "MAPINFO"]].astype(int)
            base_cpgs.MAPINFO = base_cpgs.MAPINFO - 1
            base_cpgs.columns = ["chr", "start"]

            successive_cpgs = self.manifest.loc[successive_cpgs, ["MAPINFO"]].astype(
                int
            )
            successive_cpgs.columns = ["end"]

            bed = pd.concat(
                (
                    base_cpgs.reset_index(drop=True),
                    successive_cpgs.reset_index(drop=True),
                ),
                axis=1,
            )

            bed["pair"] = pairs.index
            bed = bed.sort_values(by="chr", ascending=True)

            bed["chr"] = bed["chr"].map(lambda n_chr: f"chr{n_chr}")
            path = join(self.paths["files"], f"{mtype}.bed")

            self.paths[f"BED_{mtype}"] = path
            bed.to_csv(path, sep="\t", index=False, header=None)

    def export_genes_list(self, collections: t.Dict[str, pd.DataFrame]) -> None:
        for mtype, pairs in collections.items():

            base_cpgs = [pair.split("-")[0] for pair in pairs.index]
            successive_cpgs = [pair.split("-")[1] for pair in pairs.index]

            cpgs = list(set.union(set(base_cpgs), set(successive_cpgs)))
            genes = self.manifest.loc[cpgs, "UCSC_RefGene_Name"]
            genes = genes.str.split(";").explode().drop_duplicates().dropna()

            name = f"{mtype}_genes.csv"
            genes.to_csv(join(self.paths["files"], name))

    def homer(self, threads: int = 10) -> None:
        bg_bed = self.paths["BED_BG"]

        for mtype in self.pairs_per_marker_type.keys():
            input_bed = self.paths[f"BED_{mtype}"]
            output = join(self.paths["files"], f"HOMER_{mtype}")
            self.paths[f"HOMER_{mtype}"] = output

            if exists(output):
                print(f"Homer results for {mtype} already exists. Skipping.")
                continue

            command = [
                "findMotifsGenome.pl",
                f"{input_bed}",
                "hg19",
                f"{output}",
                "-mask",
                "-cpg",
                "-bg",
                f"{bg_bed}",
                "-h",
                "-nomotif",
                "-p",
                f"{threads}",
            ]
            subprocess.call(command, shell=False)

    def make_homer_plots(self, show: bool = True) -> t.Optional[str]:
        for mtype in self.pairs_per_marker_type.keys():
            homer_path = join(self.paths[f"HOMER_{mtype}"], "knownResults.txt")
            data = pd.read_table(homer_path, index_col=0)
            data = data[data["q-value (Benjamini)"] <= 0.05]

            if data.empty:
                print("Not applicable if less than 1 significant TF found")
                continue

            data["% of Target Sequences with Motif"] = (
                data["% of Target Sequences with Motif"]
                .str.replace("%", "")
                .astype(float)
            )
            data["% of Background Sequences with Motif"] = (
                data["% of Background Sequences with Motif"]
                .str.replace("%", "")
                .astype(float)
            )

            data["FC"] = (
                data["% of Target Sequences with Motif"]
                / data["% of Background Sequences with Motif"]
            )

            data.index = data.index.map(lambda name: name.split("/")[0])

            fig = px.bar(
                data_frame=data,
                y=data.index,
                x="FC",
                orientation="h",
                title=mtype,
            )
            fig.update_layout(
                width=500,
                height=800,
                font=dict(size=16),
                yaxis={"categoryorder": "total descending"},
            )

            fig.write_image(join(self.paths["plots"], f"homer_{mtype}.png"))
            if show:
                fig.show()

    def lola(
        self,
        lola_script: str,
        lola_db: str,
    ) -> None:
        bg_bed = self.paths["BED_BG"]

        for mtype in self.pairs_per_marker_type.keys():
            input_bed = self.paths[f"BED_{mtype}"]
            output = join(self.paths["files"], f"LOLA_{mtype}.csv")
            self.paths[f"LOLA_{mtype}"] = output

            command = [
                "Rscript",
                f"{lola_script}",
                f"{lola_db}",
                f"{input_bed}",
                f"{bg_bed}",
                f"{output}",
            ]
            subprocess.call(command, shell=False)

    def update_lola_output(self) -> None:
        for mtype in self.pairs_per_marker_type.keys():
            output = pd.read_csv(self.paths[f"LOLA_{mtype}"], index_col=0)
            output["raw p-value"] = output["pValueLog"].map(
                lambda x: 10**-x
            )  # calculate p-val form -log10(p-value)
            _, output["FDR"], _, _ = multipletests(
                output["raw p-value"], method="fdr_bh"
            )

            output["log2(oddsRatio)"] = output["oddsRatio"].map(
                lambda odds: np.log2(odds)
            )
            output["significant"] = (output["log2(oddsRatio)"].abs() > 1) & (
                output["FDR"] <= 0.05
            )

            output.to_csv(self.paths[f"LOLA_{mtype}"])

    def lola_plot(self, collection: str, show: bool = True) -> None:
        for mtype in self.pairs_per_marker_type.keys():

            output = pd.read_csv(self.paths[f"LOLA_{mtype}"], index_col=0)
            output = output[output.collection == collection]

            if output["cellType"].nunique() > 1:
                fig = px.bar(
                    data_frame=output,
                    x="description",
                    y="oddsRatio",
                    color="cellType",
                    title=mtype,
                )

            else:
                fig = px.bar(
                    data_frame=output,
                    x="description",
                    y="oddsRatio",
                    title=mtype,
                )

            fig.update_layout(
                width=1200,
                height=800,
                font=dict(size=28),
                showlegend=False,
                barmode="group",
                legend={"title": ""},
                xaxis={"categoryorder": "total descending", "title": ""},
                yaxis={"title": "Odds Ratio"},
            )

            if collection == "EPIC_Relation_To_Island":
                fig.update_xaxes(
                    categoryorder="array",
                    categoryarray=[
                        "N_Shelf",
                        "N_Shore",
                        "Island",
                        "S_Shore",
                        "S_Shelf",
                        "Opensea",
                    ],
                )

            elif collection == "EPIC_UCSC_RefGene_Group":
                fig.update_xaxes(
                    categoryorder="array",
                    categoryarray=[
                        "TSS1500",
                        "TSS200",
                        "1stExon",
                        "5_UTR",
                        "Body",
                        "ExonBnd",
                        "3_UTR",
                    ],
                )

            else:
                fig.update_yaxes(range=(0, 13))
                fig.update_xaxes(
                    categoryorder="array",
                    categoryarray=[
                        "TssA",
                        "TssAFlnk",
                        "TxFlnk",
                        "Tx",
                        "TxWk",
                        "EnhG",
                        "Enh",
                        "ZNF Rpts",
                        "Het",
                        "TssBiv",
                        "BivFlnk",
                        "EnhBiv",
                        "ReprPC",
                        "ReprPCWk",
                        "Quies",
                    ],
                )

            fig.write_image(join(self.paths["plots"], f"lola_{mtype}_{collection}.png"))
            if show:
                fig.show()

    def lola_grouped_plot(
        self,
        collections: t.Collection[str] = (
            "EPIC_Relation_To_Island",
            "EPIC_UCSC_RefGene_Group",
        ),
        mtypes_to_skip: t.Collection[str] = ("between cutoffs"),
        show: bool = True,
    ) -> None:

        for collection in collections:
            to_plot = []

            for mtype in self.pairs_per_marker_type.keys():
                if mtype == mtypes_to_skip:
                    continue

                output = pd.read_csv(self.paths[f"LOLA_{mtype}"], index_col=0)
                output = output[output.collection == collection]
                output["mType"] = mtype
                to_plot.append(output)

            to_plot = pd.concat(to_plot)
            fig = px.bar(
                to_plot, x="description", y="oddsRatio", color="mType", title=collection
            )

            fig.update_layout(
                width=1200,
                height=800,
                font=dict(size=28),
                barmode="group",
                legend={"title": ""},
                xaxis={"title": ""},
                yaxis={"title": "Odds Ratio"},
            )

            if collection == "EPIC_Relation_To_Island":
                fig.update_xaxes(
                    categoryorder="array",
                    categoryarray=[
                        "N_Shelf",
                        "N_Shore",
                        "Island",
                        "S_Shore",
                        "S_Shelf",
                        "Opensea",
                    ],
                )

            elif collection == "EPIC_UCSC_RefGene_Group":
                fig.update_xaxes(
                    categoryorder="array",
                    categoryarray=[
                        "TSS1500",
                        "TSS200",
                        "1stExon",
                        "5_UTR",
                        "Body",
                        "ExonBnd",
                        "3_UTR",
                    ],
                )
            else:
                fig.update_yaxes(range=(0, 13))
                fig.update_xaxes(categoryorder="total descending")

            fig.write_image(join(self.paths["plots"], f"lola_grouped_{collection}.png"))
            if show:
                fig.show()

    def export_metadata(self) -> None:
        with open(join(self.paths["files"], "README.txt"), "w") as file:
            for key, value in self.metadata.items():
                file.write(f"{key} --> {value}")
                file.write("\n")

    def run(self) -> None:
        self.__build_output_dir()
        self.__update_mask()

        self.fit_model_per_group()
        self.build_comparison_frame()
        self.mark_types_of_markers()
        self.mark_co_methylated_DMRS()

        self.create_delta_norm()
        self.calculate_type_to_type_stats()

        self.mark_markers()
        self.mark_specific_markers()
        self.export_frames_per_type()
        self.comparison_frame.to_csv(join(self.paths["files"], f"comparison_frame.csv"))

        self.plots(self.pairs_per_marker_type, limit=100)
        self.cluster_map(self.pairs_per_marker_type)
        self.intersection_map()

        self.export_bed(self.pairs_per_marker_type)
        self.export_bed({"BG": self.mask.set_index("Pair")})  # Export BG

        self.export_genes_list(self.pairs_per_marker_type)
        self.export_genes_list({"BG": self.mask.set_index("Pair")})  # Export BG

        self.homer()
        self.make_homer_plots()

        self.lola()
        self.update_lola_output()
        self.lola_plot(collection="WBC_15_state_core_model")
        self.lola_plot(collection="EPIC_Relation_To_Island")
        self.lola_grouped_plot()

        self.export_metadata()
        print("THE END")
