# Unit tests

Run the current standard-library tests with:

```bash
python -m unittest discover -s tests/unit -v
```

The tests cover config layering and portable templates, complete and partial Space Ranger layouts, parent `outs/` resolution, canonical spatial coordinates, H5/MTX ingestion, numeric, complexity and spatial QC figures, exact image/scalefactor matching, H&E alignment-review placeholders, coordinate fallback and validation, unavailable mitochondrial metrics, H5AD round trips, the split-rule target layout, and the 20-item manual-review record contract.
