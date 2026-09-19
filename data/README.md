# data/

Nothing in this directory is tracked by git. The eighteen deposited Q Exactive
runs of PXD001584 total 47.3 GiB before conversion, measured from the PRIDE
API rather than estimated.

Layout after a full run:

```
data/
  fasta/    search_database.fasta, the UniProt target proteome, contaminants,
            and the entrapment database built by 06_entrapment.sh
  raw/      Thermo RAW files, transient. 04_convert_mzml.sh deletes each one
            as soon as its mzML has been written and verified
  mzml/     gzipped centroided mzML, retained. Sage needs every run present
            at once for retention time alignment and label free quantification,
            so these cannot be streamed away one at a time
```

Peak disk is the retained mzML set plus the single RAW being converted at that
moment. For all eighteen runs that projects to roughly 14 GiB. `00_configure.sh`
computes the projection for your machine and refuses to start a run that cannot
finish, which is the failure this pipeline is built to avoid.

To put the data on another filesystem, make `data/` a symlink. Nothing else in
the pipeline assumes it lives inside the repository.
