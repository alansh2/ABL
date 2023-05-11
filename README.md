# Atmospheric Boundary Layer (ABL)
This repository contains the codes for AE 514 Project 2

## Conventionally Neutral ABL
The conventionally neutral ABL is characterized by zero surface temperature potential flux and a temperature inversion cap. Field data of the atmosphere up to 30,000 m is sourced from https://www1.ncdc.noaa.gov/pub/data/igra/. Only recent soundings conducted at midnight are used. File downloads from `data/data-y2d` will be of the form `Station-ID-data.txt`.

__`cleanData.m`__
The script `cleanData.m` reads from `Station-ID-data.txt` and outputs `Station-ID-data-beg2021.txt` containing only the sounding data collected at midnight.

## Stable Surface Layer
