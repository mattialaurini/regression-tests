# NIST data files

`BoxBOD.dat`, `MGH17.dat`, and `Lanczos2.dat` include the official data blocks so that smoke tests and the first pilot cases can run offline.

Run the following MATLAB command to obtain all files directly from the official
NIST StRD site:

```matlab
download_nist_data();
```

Downloaded files are placed in `data/raw/`. The parser always reads the final
NIST data block and checks the expected number of observations and predictors.
