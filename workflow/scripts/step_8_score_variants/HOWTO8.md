
# howto run it
### i.e. how i did it
----

## TODO:
- add authors and in snakerule 
- add train test score constraint for amount of chromosomes, before processing data

```
wildcard_constraints:   
     part="[chr][0-9a-zA-z_]+",

# Function to gather all outputs from checkpoint 
CHROMSOME_LIST = config['chromosomes']['score']
```

- TODO: what to be temp. and what to save, more intermediate files that can be removed?
- TODO: move output of problem_out; see why problem and if they can be incorporated later
- TODO: intermediate zipping?


pandas.errors.ParserError: Error tokenizing data. C error: Expected 8 fields in line 3, saw 23

