
# Overview
Provides an example of the output of the `determineOrbit()` function
Note: The actual output of the function is a string. In python `json.loads()` should be used to format the data for use as a dict


# Example Output
The output is a dictionary of dictionaries. The basic structure is as follows: 
```python
{
  "Filter: <FILTER_NAME>
  // Details what filter was used for generating data
  "Estimation":{
  // Contains all the data on estimated states and covariances
  }, 
  "Propagation":{
  // Gives the basic starting state and Time used for 
  }
}
```

## Filter Example
```python
{
"Filter": "UKF"
}
```

## Propagation Example
```python
{
  "Propagation": { 
     "Time": '2018-01-02T00:02:54Z'
     "State": list[1xN Floats]
  }
}
```

## Estimation 
The Estimation section is a list of dictionaries containing the estimates at each measurement epoch

```python
[{
 'EstimatedCovariance': [[Nx1 Floats]]
 'EstimatedState': [Nx1 Floats],
 'InnovationCovariance': [[MxM Floats]],
 'PostFit': {'Range': [0.703046889891979], 'RangeRate': [37701413.22956683]},
 'PreFit': {'Range': [0.703002641655857], 'RangeRate': [37701412.657939374]},
 'Time': '2018-01-01T23:52:54.000Z'
 }, 
 ...
 ]
```
