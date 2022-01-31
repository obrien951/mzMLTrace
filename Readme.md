To install mzMLTrace, first obain the code by running

```
cd TOP_DIR && git clone https://github.com/obrien951/mzMLTrace
```

To compile the code, while in `TOP_DIR/mzMLTrace`, run

```
cmake .
make
```

Then to make mzMLTrace available by import, run:

```
python3 -m pip install -e .
```


Data can be extracted by the code as follows:

```
    k = mzMLTrace.specToChrom()
    mt_lists = []
    for mzml in mzmlFiles:
        chrommzml = mzml.removesuffix(".mzML")
        if k.is_set():
            k.reset()
        k.set_filename(mzml)
        k.set_minimum_intensity(PARAMETERS["min_intensity_threshold"])
        k.readSpectra()
        k.findChromatograms()
        mt_lists.append([])
        ind = len(mt_lists) - 1
        for i in range(k.get_nchrom()):
            count = k.get_chrom_len(i)
            mz = [0.0 for i in range(count) ]
            RT = [0.0 for i in range(count) ]
            INTS = [0.0 for i in range(count) ]
            k.get_chrom(i, mz, INTS, RT)
            mt_lists[ind].append(ext_MassTrace())

```
