# FindPeek
Find peak areas in location series  
* The first method is to judge the region as a peak according to the continuous rise and fall.  
* the second method is to give a background value (CutMethod = T). If it is higher than the background
 value continuously, the region is considered as a peak. The background value can be given manually or 
 automatically generated according to the maximum value 
## installation
### Checkout the latest release of FindPeek from GitHub
```https://github.com/yiliu4234/FindPeek.git```
### Install R dependencies (in R)
 ```install.packages("data.table") # version > 1.10.4```

### Install the FindPeek package
from the command line and in the directory where FindPeek github was cloned.
```R CMD INSTALL FindPeek ```

