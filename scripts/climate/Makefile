SHELL := /bin/bash

# This makefile is designed to allow reextracting of all climate data by simply
# running the command: "make all" or "make" in the directory where this makefile
# lives. Data are output to $DATA_DIR

# User Options - you can probably leave this alone.
DATA_DIR=../../derived_data

################################################################################
# Recipes
################################################################################

# Default target
all: daymet

daymet: ../../derived_data/sunflower_data_simple_v1.csv
	@echo 'Extracting DayMet data at geocoordinates...'
	@Rscript 0_01_extract_daymet_climate_data.R $(DATA_DIR)
	@echo 'Updated DayMet time series data.'

### ----------------------------------------------------------------------------
.PHONY: all daymet
