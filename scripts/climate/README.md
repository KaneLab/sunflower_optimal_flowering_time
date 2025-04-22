# Climate data compilation
To download and compile the *daymet_timeseries_cleaned.csv* dataset, use the
`make` command from this directory:

For example, first ensure you are in this directory:
```bash
pwd
#USER/sunflower_optimal_flowering_time/scripts/climate
```

Then run, `make`:
```bash
make
```

The script *0_01_extract_daymet_climate_data.R* will then install any missing
R packages (see **Requirements** below), download the daymet data into an
intermediate data one site at a time, and finally compile all the data into a
single file: *daymet_timeseries_cleaned.csv* which will be output into
`../../derived_data/`

## Requirements
> Ensure you have GNU Make installed. Run `make -v`. If you do not have
> make, you may install it for free [here](https://www.gnu.org/software/make/).

### Command line tools
- GNU bash, version ~ 5.1.16(1)-release (tested on x86_64-pc-linux-gnu)
- GNU Make ~ 4.3
- R ~ 4.3.2

### R Packages:
- `dplyr` >=1.1.4
- `lubridate` >=1.9.3
- `stringr` >= 1.5.1
- `daymetr` 1.7.1
