# Atmospheric Boundary Layer (ABL)
This repository contains the codes for AE 514 Project 2

## Conventionally Neutral ABL
The conventionally neutral ABL is characterized by zero surface temperature potential flux and a temperature inversion cap. Field data of the atmosphere up to 30,000 m is sourced from https://www1.ncdc.noaa.gov/pub/data/igra/. Only recent soundings conducted at midnight are used. File downloads from `data/data-y2d` will be of the form `Station-ID-data.txt`.

### `cleanData.m` ###
The script `cleanData.m` reads from `Station-ID-data.txt` and outputs `Station-ID-data-beg2021.txt` containing only the sounding data collected at midnight.

### `dataToModel.m` ###
The script `dataToModel.m` is the top-level program that calls subroutines to fetch raw data and apply the CNABL model to those data. Since the ABL is fully defined within the Troposphere, all data are truncated to the altitude of Tropopause. Since the IGRA data reports geopotential altitude instead of geometric altitude, this conversion is performed and saved in the output structure array. Both these procedures are done in the subroutine `getData.m`.

The friction velocity is computed via regression of the velocity magnitude in the subroutine `bestFricVel.m`. It finds the friction velocity that minimizes the absolute error between the measurements and the prediction from `getModel.m`, which is the implementation of the CNABL model proposed by Liu and Stevens (2022). Following this, the coordinate axes are oriented such that the cross-flow velocity, *V*, approaches zero at the surface. This is done via call to the subroutine `orientAxes.m`.

## Stable Surface Layer
