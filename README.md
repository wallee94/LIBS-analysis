# LIBS-analysis

This library is NOT completed yet. It can be used to fit lorentzian functions to experimental data, but most of the methods to improve this fit are still not working.

## Instructions

You should start installing the requirements:

    pip install -r requirements.txt

To make an analysis using the configuration file, just run
    
    python script.py
    
Using python3. 

You can edit the configuration changing the params.txt file.

## Retrieving data from NIST

In the configuration file you can add elements to the "elements_from_NIST" array. This elements will be retrieved from NIST online database public here:

https://physics.nist.gov/PhysRefData/ASD/lines_form.html
