
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
