
# Overview
Provides an example of the output of the `determineOrbit()` function
Note: The actual output of the function is a string. In python `json.loads()` should be used to format the data for use as a dict


# Example Output
The output is a dictionary of dictionaries each representing a 

Note: Once you have loaded the file into a python dictionary the data will not remain ordered. To get ordered data do the following
```python
outputs = json.loads(inputs, object_pairs_hook=OrderedDict)  
```
This will load the OD output data into an OrderedDict

{{
 'EstimatedCovariance': [[Nx1 Floats]]
 'EstimatedState': [Nx1 Floats],
 'InnovationCovariance': [[MxM Floats]],
 'PostFit': {'Range': [0.703046889891979], 'RangeRate': [37701413.22956683]},
 'PreFit': {'Range': [0.703002641655857], 'RangeRate': [37701412.657939374]},
 'Time': '2018-01-01T23:52:54.000Z'
 }, 
 ...
 }
