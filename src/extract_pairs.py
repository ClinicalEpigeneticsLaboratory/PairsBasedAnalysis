import pandas as pd
from tqdm import tqdm


class ExtractPairs:
    def __init__(self, manifest: str, hg19: str):
        manifest = pd.read_parquet(manifest, columns=["CHR", "MAPINFO"]).dropna()
        manifest = manifest[~manifest["CHR"].isin(["X", "Y"])]
        manifest = manifest.loc[[cg for cg in manifest.index if cg.startswith("cg")]]
        manifest["MAPINFO"] = manifest["MAPINFO"].astype(int)
        manifest["CHR"] = manifest["CHR"].map(lambda n: f"chr{n}")

        self.hg19 = pd.read_parquet(hg19)
        self.manifest = manifest
        self.raw_pairs = {}
        self.selected_pairs = {}

    def find_pairs(self, threshold: int = 50) -> None:
        for chr_ in tqdm(self.manifest["CHR"].unique()):

            partial_frame = self.manifest[self.manifest["CHR"] == chr_][
                "MAPINFO"
            ].astype(int)

            partial_frame = partial_frame.sort_values(ascending=True)
            base = partial_frame.to_frame().reset_index()
            base.columns = ["BaseCpg", "BasePos"]

            successive = (
                partial_frame.drop(partial_frame.index[0])
                .to_frame()
                .dropna()
                .reset_index()
            )

            successive.columns = ["SuccessiveCpg", "SuccessivePos"]
            distance_frame = pd.concat((base, successive), axis=1).dropna()

            distance_frame["Pair"] = (
                distance_frame.BaseCpg + "-" + distance_frame.SuccessiveCpg
            )
            distance_frame["Distance"] = (
                distance_frame.SuccessivePos - distance_frame.BasePos
            )

            selected_pairs = distance_frame[distance_frame.Distance <= threshold]
            self.raw_pairs[chr_] = selected_pairs

    def remove_non_consecutive(self) -> None:
        for chr_, selected_pairs in self.raw_pairs.items():
            partial_hg19 = self.hg19[self.hg19.chr == chr_]
            selected_pairs_per_chr = []

            for _, record in tqdm(selected_pairs.iterrows(), desc=chr_):
                pair, start, end = (
                    record["Pair"],
                    record["BasePos"],
                    record["SuccessivePos"],
                )
                cg_in_range = partial_hg19[
                    (partial_hg19.start > start) & (partial_hg19.start < end)
                ]

                if cg_in_range.empty:
                    selected_pairs_per_chr.append(pair)

            self.selected_pairs[chr_] = selected_pairs[
                selected_pairs.Pair.isin(selected_pairs_per_chr)
            ]

    @property
    def get_raw_pairs(self) -> pd.DataFrame:
        frame = []
        for chr_, partial_frame in tqdm(self.raw_pairs.items()):
            partial_frame["CHR"] = chr_
            frame.append(partial_frame)

        frame = pd.concat(frame)[
            ["Pair", "BasePos", "SuccessivePos", "Distance", "CHR"]
        ]
        return frame

    @property
    def get_selected_pairs(self) -> pd.DataFrame:
        frame = []
        for chr_, partial_frame in tqdm(self.selected_pairs.items()):
            partial_frame["CHR"] = chr_
            frame.append(partial_frame)

        frame = pd.concat(frame)[
            ["Pair", "BasePos", "SuccessivePos", "Distance", "CHR"]
        ]
        return frame
