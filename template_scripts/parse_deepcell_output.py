import json, tifffile
from skimage.measure import regionprops_table
import pandas as pd

# Parse division file
with open("divisions.json", "r") as f:
    division_data = json.load(f)
division_data = pd.DataFrame(division_data)
division_data[["daughter_1", "daughter_2"]] = pd.DataFrame(
    division_data["daughters"].tolist(), index=division_data.index
)
division_data = division_data.drop(columns="daughters")
# Note that deeplabel time starts from 0, so add 1 to t
division_data["t"] = division_data["t"] + 1
# Parse cell annotation file
with open("cellTypes.json", "r") as f:
    cell_annotation = json.load(f)
# Parse mask image and compute position
mask = tifffile.imread("y.ome.tiff")
tracking_table = None
for i in range(mask.shape[0]):
    mask_slice = mask[i, ...]
    props = pd.DataFrame(
        regionprops_table(
            mask_slice,
            properties=("label", "centroid"),
        )
    )
    props["t"] = i + 1
    if tracking_table is None:
        tracking_table = props
    else:
        tracking_table = pd.concat([tracking_table, props], axis=0)
# merge division and annotation
# join daughters first
tracking_table = (
    tracking_table.merge(
        division_data.drop(columns="t"), left_on="label", right_on="parent", how="left"
    )
    .drop(columns="parent")
    .rename(columns={"daughter_1": "daughter_1_f", "daughter_2": "daughter_2_f"})
)
# join parent
tracking_table = (
    tracking_table.merge(
        division_data.drop(columns="t"),
        left_on="label",
        right_on="daughter_1",
        how="left",
    )
    .rename(columns={"parent": "parent_1"})
    .drop(columns=["daughter_1", "daughter_2"])
)
# Assemble result table
tracking_table = (
    tracking_table.merge(
        division_data.drop(columns="t"),
        left_on="label",
        right_on="daughter_2",
        how="left",
    )
    .rename(columns={"parent": "parent_2"})
    .drop(columns=["daughter_1", "daughter_2"])
)
# Resolve parent and daughter columns
tracking_table["parent"] = (
    tracking_table["parent_1"]
    .combine_first(tracking_table["parent_2"])
    .fillna(0)
    .astype(int)
)
tracking_table = tracking_table.drop(columns=["parent_1", "parent_2"]).rename(
    columns={"daughter_1_f": "daughter_1", "daughter_2_f": "daughter_2"}
)
# Annotation
tracking_table["annotation"] = "Unknown"
for i in range(len(cell_annotation)):
    lab_name = cell_annotation[i]["name"]
    tracking_table.loc[
        pd.Series.isin(tracking_table["label"], cell_annotation[i]["cells"]),
        "annotation",
    ] = lab_name
# Export to csv
tracking_table.to_csv("tracking_result.csv")
