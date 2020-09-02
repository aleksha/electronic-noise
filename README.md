
# electronic-noise

This project contains a tool to create a spectrum of the electronic noise
for the Flash ADC used with ACTAF. The code is tested with Python 3, but 
is expected to run with Python 2.7 as well.

## How to use tool

To generate events:
```bash
sudo apt install wget git python3-numpy python3-scipy python3-matplotlib
git clone https://github.com/aleksha/electronic-noise.git
cd electronic-noise
wget http://adzyuba.web.cern.ch/adzyuba/d/dump.txt
python3 anode_noise.py
```

Generated events will be stored into **_noise_events.data_** file.
Each line will contain a ADC value separated by space-character.

## Documentation

For the deatils see **_doc/MC4NOISE.docx_**.

## How to produce dump-file

 1. Ask Alexey Dzyuba or Alexander Inglessi for the root-file with data
 2. Use **_dump/no_gui.C_** script to process it.

Scripts has been tested with CERN ROOT v5.36.

## Resolution

First one have to convert CSV file to the one with space-separator.
After one can use ROOT script `resolution.C`.

```bash
cp dump.txt Dump.txt
sed -i 's/,/ /g' Dump.txt
root -b -q resolution.C
```

1 channel = 27.5 eV (if one use 5.5 MeV alpha calibration => 200000 ch.)

For the real data results are:
 * Nev  = 1000 events
 * GATE = 6 us
 * MEAN = 1.22889e+06 a.u.
 * SDV. = 762.773 a.u.
 * RESO = 20.9763 keV

For the simulated noise with same gate / number of events condition:
 * MEAN = 1.22845e+06 a.u.
 * SDV. = 776.29 a.u.
 * RESO = 21.348 keV

They are in reasonable agreement.
