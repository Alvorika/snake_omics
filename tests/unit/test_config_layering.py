import copy
import unittest
from pathlib import Path

import pandas as pd
from snakemake.common.configfile import load_configfile
from snakemake.utils import update_config, validate


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULTS_PATH = REPOSITORY_ROOT / "config" / "defaults.yaml"
ACTIVE_CONFIG_PATH = REPOSITORY_ROOT / "config" / "config.yaml"
TEMPLATE_PATH = REPOSITORY_ROOT / "config" / "config.template.yaml"
SCHEMA_PATH = REPOSITORY_ROOT / "workflow" / "schemas" / "config.schema.yaml"
SAMPLE_TEMPLATE_PATH = REPOSITORY_ROOT / "config" / "samples.template.tsv"
SAMPLE_SCHEMA_PATH = REPOSITORY_ROOT / "workflow" / "schemas" / "samples.schema.yaml"
PATHWAY_TEMPLATE_PATH = REPOSITORY_ROOT / "config" / "pathway_gene_sets.template.tsv"
QC_PROFILE_PATH = (
    REPOSITORY_ROOT / "config" / "qc_profiles" / "unconfigured_v1.yaml"
)

RECOMMENDED_SAMPLE_FIELDS = {
    "animal_id",
    "biological_replicate",
    "technical_batch",
    "slide_id",
    "capture_area",
    "library_id",
    "assay",
    "species",
    "genome_reference",
    "probe_set",
    "probe_set_checksum",
    "section_level",
    "orientation",
}


class ConfigLayeringTests(unittest.TestCase):
    def test_active_project_config_is_not_shipped(self) -> None:
        self.assertFalse(ACTIVE_CONFIG_PATH.exists())

    def test_nested_override_does_not_remove_sibling_defaults(self) -> None:
        effective = copy.deepcopy(load_configfile(DEFAULTS_PATH))
        update_config(
            effective,
            {"qc": {"numeric_metrics": {"detected_genes": False}}},
        )

        self.assertFalse(effective["qc"]["numeric_metrics"]["detected_genes"])
        self.assertTrue(effective["qc"]["numeric_metrics"]["total_counts"])
        self.assertTrue(effective["qc"]["numeric_metrics"]["in_tissue"])

    def test_template_is_a_valid_override_of_defaults(self) -> None:
        effective = load_configfile(DEFAULTS_PATH)
        update_config(effective, load_configfile(TEMPLATE_PATH))
        validate(effective, SCHEMA_PATH)

        self.assertEqual(
            effective["project"]["name"],
            "my-spatial-transcriptomics-project",
        )
        self.assertEqual(effective["execution"]["python"], "python")
        self.assertIsNone(effective["analysis"]["condition"]["genotype_reference"])
        self.assertIsNone(effective["analysis"]["pathway"]["expected_rankings"])
        self.assertEqual(effective["modules"]["enabled"], ["qc"])

    def test_sample_template_matches_schema_and_recommended_contract(self) -> None:
        samples = pd.read_csv(
            SAMPLE_TEMPLATE_PATH,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )
        validate(samples, SAMPLE_SCHEMA_PATH)

        self.assertTrue(RECOMMENDED_SAMPLE_FIELDS.issubset(samples.columns))
        self.assertEqual(set(samples["input_type"]), {"spaceranger"})
        self.assertTrue(samples["input_path"].str.startswith("../").all())
        self.assertTrue(samples["roi_path"].str.startswith("../").all())

        optional_roi = samples.copy()
        optional_roi.loc[0, [
            "roi_path",
            "roi_barcode_column",
            "roi_label_column",
        ]] = ""
        validate(optional_roi, SAMPLE_SCHEMA_PATH)

    def test_pathway_manifest_template_is_disabled_and_portable(self) -> None:
        manifest = pd.read_csv(
            PATHWAY_TEMPLATE_PATH,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )

        self.assertFalse(manifest.empty)
        self.assertEqual(set(manifest["enabled"].str.lower()), {"no"})
        self.assertTrue(
            manifest["gmt_path"]
            .map(lambda value: not Path(value).is_absolute())
            .all()
        )
        self.assertTrue(manifest["sha256"].str.startswith("REPLACE_").all())

    def test_default_qc_profile_is_explicitly_uncalibrated(self) -> None:
        profile = load_configfile(QC_PROFILE_PATH)
        self.assertEqual(profile["profile_id"], "unconfigured_v1")
        for metric in (
            "total_counts",
            "detected_genes",
            "mitochondrial_fraction",
        ):
            self.assertTrue(
                all(
                    value is None
                    for value in profile["thresholds"][metric].values()
                )
            )


if __name__ == "__main__":
    unittest.main()
