\name{RadOnc}
\alias{johndoe}
\alias{janedoe}
\alias{RadOnc-data}
\docType{data}
\title{
DVH data for 'John Doe' and 'Jane Doe'
}
\description{
List of 10 DVH objects corresponding to structures (e.g. liver) for patients 'John Doe' and 'Jane Doe' }
\usage{data(package="RadOnc")}
\format{
  The format for \code{johndoe} is:
\preformatted{
Formal class 'DVH.list' [package "RadOnc"] with 1 slots
  ..@ structures:List of 10
  .. ..$ LIVER       :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "LIVER"
  .. .. .. ..@ structure.volume: num 1367
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 92.9
  .. .. .. ..@ dose.min        : num 0
  .. .. .. ..@ dose.mean       : num 4.9
  .. .. .. ..@ dose.median     : num 0
  .. .. .. ..@ dose.mode       : num 0
  .. .. .. ..@ dose.STD        : num 10.8
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1026] 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
  .. .. .. ..@ dose.type       : chr "relative"
  .. .. .. ..@ volumes         : num [1:1026] 1367 414 392 378 368 ...
  .. .. .. ..@ volume.type     : chr "absolute"
  .. ..$ SMALL_BOWEL :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "SMALL_BOWEL"
  .. .. .. ..@ structure.volume: num 206
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 99.8
  .. .. .. ..@ dose.min        : num 0
  .. .. .. ..@ dose.mean       : num 4.3
  .. .. .. ..@ dose.median     : num 0
  .. .. .. ..@ dose.mode       : num 0
  .. .. .. ..@ dose.STD        : num 15.3
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1026] 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
  .. .. .. ..@ dose.type       : chr "relative"
  .. .. .. ..@ volumes         : num [1:1026] 206.2 43.1 39.5 37.4 35.8 ...
  .. .. .. ..@ volume.type     : chr "absolute"
  .. ..$ DUODENUM    :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "DUODENUM"
  .. .. .. ..@ structure.volume: num 93.1
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 102
  .. .. .. ..@ dose.min        : num 0
  .. .. .. ..@ dose.mean       : num 70.7
  .. .. .. ..@ dose.median     : num 81.3
  .. .. .. ..@ dose.mode       : num 100
  .. .. .. ..@ dose.STD        : num 30.3
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1026] 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
  .. .. .. ..@ dose.type       : chr "relative"
  .. .. .. ..@ volumes         : num [1:1026] 93.1 93.1 93 93 92.9 ...
  .. .. .. ..@ volume.type     : chr "absolute"
  .. ..$ STOMACH     :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "STOMACH"
  .. .. .. ..@ structure.volume: num 304
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 101
  .. .. .. ..@ dose.min        : num 0
  .. .. .. ..@ dose.mean       : num 7.9
  .. .. .. ..@ dose.median     : num 0
  .. .. .. ..@ dose.mode       : num 0
  .. .. .. ..@ dose.STD        : num 22.8
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1026] 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
  .. .. .. ..@ dose.type       : chr "relative"
  .. .. .. ..@ volumes         : num [1:1026] 303.6 78 71.5 67.8 65.1 ...
  .. .. .. ..@ volume.type     : chr "absolute"
  .. ..$ CTV         :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "CTV"
  .. .. .. ..@ structure.volume: num 88.4
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 103
  .. .. .. ..@ dose.min        : num 96.8
  .. .. .. ..@ dose.mean       : num 100
  .. .. .. ..@ dose.median     : num 100
  .. .. .. ..@ dose.mode       : num 100
  .. .. .. ..@ dose.STD        : num 0.7
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1026] 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
  .. .. .. ..@ dose.type       : chr "relative"
  .. .. .. ..@ volumes         : num [1:1026] 88.4 88.4 88.4 88.4 88.4 ...
  .. .. .. ..@ volume.type     : chr "absolute"
  .. ..$ PTV         :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "PTV"
  .. .. .. ..@ structure.volume: num 156
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 103
  .. .. .. ..@ dose.min        : num 84.1
  .. .. .. ..@ dose.mean       : num 99.6
  .. .. .. ..@ dose.median     : num 99.8
  .. .. .. ..@ dose.mode       : num 99.9
  .. .. .. ..@ dose.STD        : num 1.2
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1026] 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
  .. .. .. ..@ dose.type       : chr "relative"
  .. .. .. ..@ volumes         : num [1:1026] 156 156 156 156 156 ...
  .. .. .. ..@ volume.type     : chr "absolute"
  .. ..$ BODY        :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "BODY"
  .. .. .. ..@ structure.volume: num 17666
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 103
  .. .. .. ..@ dose.min        : num 0
  .. .. .. ..@ dose.mean       : num 5.6
  .. .. .. ..@ dose.median     : num 0
  .. .. .. ..@ dose.mode       : num 0
  .. .. .. ..@ dose.STD        : num 17.1
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1026] 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
  .. .. .. ..@ dose.type       : chr "relative"
  .. .. .. ..@ volumes         : num [1:1026] 17666 3424 3290 3206 3143 ...
  .. .. .. ..@ volume.type     : chr "absolute"
  .. ..$ LEFT_KIDNEY :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "LEFT_KIDNEY"
  .. .. .. ..@ structure.volume: num 154
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 44.4
  .. .. .. ..@ dose.min        : num 0
  .. .. .. ..@ dose.mean       : num 4.7
  .. .. .. ..@ dose.median     : num 0
  .. .. .. ..@ dose.mode       : num 0
  .. .. .. ..@ dose.STD        : num 9.5
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1026] 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
  .. .. .. ..@ dose.type       : chr "relative"
  .. .. .. ..@ volumes         : num [1:1026] 154.2 48.3 46.4 45.3 44.4 ...
  .. .. .. ..@ volume.type     : chr "absolute"
  .. ..$ RIGHT_KIDNEY:Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "RIGHT_KIDNEY"
  .. .. .. ..@ structure.volume: num 155
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 98.5
  .. .. .. ..@ dose.min        : num 0
  .. .. .. ..@ dose.mean       : num 24.2
  .. .. .. ..@ dose.median     : num 24.3
  .. .. .. ..@ dose.mode       : num 0
  .. .. .. ..@ dose.STD        : num 22.8
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1026] 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
  .. .. .. ..@ dose.type       : chr "relative"
  .. .. .. ..@ volumes         : num [1:1026] 155 134 132 130 129 ...
  .. .. .. ..@ volume.type     : chr "absolute"
  .. ..$ CORD        :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "CORD"
  .. .. .. ..@ structure.volume: num 40.7
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 55.5
  .. .. .. ..@ dose.min        : num 0
  .. .. .. ..@ dose.mean       : num 20.7
  .. .. .. ..@ dose.median     : num 4.1
  .. .. .. ..@ dose.mode       : num 0
  .. .. .. ..@ dose.STD        : num 23.2
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1026] 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
  .. .. .. ..@ dose.type       : chr "relative"
  .. .. .. ..@ volumes         : num [1:1026] 40.7 23 22.7 22.4 22.2 ...
  .. .. .. ..@ volume.type     : chr "absolute"
}


  The format for \code{janedoe} is:
\preformatted{
Formal class 'DVH.list' [package "RadOnc"] with 1 slots
  ..@ structures:List of 10
  .. ..$ LIVER       :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "LIVER"
  .. .. .. ..@ structure.volume: num 1636
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 5634
  .. .. .. ..@ dose.min        : num 42.7
  .. .. .. ..@ dose.mean       : num 707
  .. .. .. ..@ dose.median     : num 276
  .. .. .. ..@ dose.mode       : num 99.5
  .. .. .. ..@ dose.STD        : num 917
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1133] 0 5 10 15 20 25 30 35 40 45 ...
  .. .. .. ..@ dose.type       : chr "absolute"
  .. .. .. ..@ volumes         : num [1:1133] 100 100 100 100 100 ...
  .. .. .. ..@ volume.type     : chr "relative"
  .. ..$ LEFT_KIDNEY :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "LEFT_KIDNEY"
  .. .. .. ..@ structure.volume: num 196
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 3847
  .. .. .. ..@ dose.min        : num 75.8
  .. .. .. ..@ dose.mean       : num 1021
  .. .. .. ..@ dose.median     : num 703
  .. .. .. ..@ dose.mode       : num 119
  .. .. .. ..@ dose.STD        : num 802
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1133] 0 5 10 15 20 25 30 35 40 45 ...
  .. .. .. ..@ dose.type       : chr "absolute"
  .. .. .. ..@ volumes         : num [1:1133] 100 100 100 100 100 100 100 100 100 100 ...
  .. .. .. ..@ volume.type     : chr "relative"
  .. ..$ STOMACH     :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "STOMACH"
  .. .. .. ..@ structure.volume: num 695
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 5353
  .. .. .. ..@ dose.min        : num 59
  .. .. .. ..@ dose.mean       : num 1280
  .. .. .. ..@ dose.median     : num 1302
  .. .. .. ..@ dose.mode       : num 111
  .. .. .. ..@ dose.STD        : num 1062
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1133] 0 5 10 15 20 25 30 35 40 45 ...
  .. .. .. ..@ dose.type       : chr "absolute"
  .. .. .. ..@ volumes         : num [1:1133] 100 100 100 100 100 100 100 100 100 100 ...
  .. .. .. ..@ volume.type     : chr "relative"
  .. ..$ DUODENUM    :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "DUODENUM"
  .. .. .. ..@ structure.volume: num 34.2
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 5620
  .. .. .. ..@ dose.min        : num 2708
  .. .. .. ..@ dose.mean       : num 4755
  .. .. .. ..@ dose.median     : num 4943
  .. .. .. ..@ dose.mode       : num 5365
  .. .. .. ..@ dose.STD        : num 635
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1133] 0 5 10 15 20 25 30 35 40 45 ...
  .. .. .. ..@ dose.type       : chr "absolute"
  .. .. .. ..@ volumes         : num [1:1133] 100 100 100 100 100 100 100 100 100 100 ...
  .. .. .. ..@ volume.type     : chr "relative"
  .. ..$ RIGHT_KIDNEY:Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "RIGHT_KIDNEY"
  .. .. .. ..@ structure.volume: num 224
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 4202
  .. .. .. ..@ dose.min        : num 102
  .. .. .. ..@ dose.mean       : num 1511
  .. .. .. ..@ dose.median     : num 1624
  .. .. .. ..@ dose.mode       : num 1638
  .. .. .. ..@ dose.STD        : num 714
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1133] 0 5 10 15 20 25 30 35 40 45 ...
  .. .. .. ..@ dose.type       : chr "absolute"
  .. .. .. ..@ volumes         : num [1:1133] 100 100 100 100 100 100 100 100 100 100 ...
  .. .. .. ..@ volume.type     : chr "relative"
  .. ..$ CTV         :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "CTV"
  .. .. .. ..@ structure.volume: num 147
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 5647
  .. .. .. ..@ dose.min        : num 5169
  .. .. .. ..@ dose.mean       : num 5500
  .. .. .. ..@ dose.median     : num 5505
  .. .. .. ..@ dose.mode       : num 5500
  .. .. .. ..@ dose.STD        : num 59.8
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1133] 0 5 10 15 20 25 30 35 40 45 ...
  .. .. .. ..@ dose.type       : chr "absolute"
  .. .. .. ..@ volumes         : num [1:1133] 100 100 100 100 100 100 100 100 100 100 ...
  .. .. .. ..@ volume.type     : chr "relative"
  .. ..$ PTV         :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "PTV"
  .. .. .. ..@ structure.volume: num 239
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 5665
  .. .. .. ..@ dose.min        : num 4750
  .. .. .. ..@ dose.mean       : num 5471
  .. .. .. ..@ dose.median     : num 5493
  .. .. .. ..@ dose.mode       : num 5500
  .. .. .. ..@ dose.STD        : num 98.6
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1133] 0 5 10 15 20 25 30 35 40 45 ...
  .. .. .. ..@ dose.type       : chr "absolute"
  .. .. .. ..@ volumes         : num [1:1133] 100 100 100 100 100 100 100 100 100 100 ...
  .. .. .. ..@ volume.type     : chr "relative"
  .. ..$ SMALL_BOWEL :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "SMALL_BOWEL"
  .. .. .. ..@ structure.volume: num 232
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 4934
  .. .. .. ..@ dose.min        : num 59.6
  .. .. .. ..@ dose.mean       : num 294
  .. .. .. ..@ dose.median     : num 171
  .. .. .. ..@ dose.mode       : num 127
  .. .. .. ..@ dose.STD        : num 391
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1133] 0 5 10 15 20 25 30 35 40 45 ...
  .. .. .. ..@ dose.type       : chr "absolute"
  .. .. .. ..@ volumes         : num [1:1133] 100 100 100 100 100 100 100 100 100 100 ...
  .. .. .. ..@ volume.type     : chr "relative"
  .. ..$ CORD        :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "CORD"
  .. .. .. ..@ structure.volume: num 64.9
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 3443
  .. .. .. ..@ dose.min        : num 0
  .. .. .. ..@ dose.mean       : num 732
  .. .. .. ..@ dose.median     : num 136
  .. .. .. ..@ dose.mode       : num 6.9
  .. .. .. ..@ dose.STD        : num 953
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1133] 0 5 10 15 20 25 30 35 40 45 ...
  .. .. .. ..@ dose.type       : chr "absolute"
  .. .. .. ..@ volumes         : num [1:1133] 100 98.5 95.9 93.8 90.8 ...
  .. .. .. ..@ volume.type     : chr "relative"
  .. ..$ BODY        :Formal class 'DVH' [package "RadOnc"] with 18 slots
  .. .. .. ..@ structure.name  : chr "BODY"
  .. .. .. ..@ structure.volume: num 25508
  .. .. .. ..@ type            : chr "cumulative"
  .. .. .. ..@ dose.max        : num 5665
  .. .. .. ..@ dose.min        : num 0
  .. .. .. ..@ dose.mean       : num 478
  .. .. .. ..@ dose.median     : num 66.4
  .. .. .. ..@ dose.mode       : num 0.4
  .. .. .. ..@ dose.STD        : num 946
  .. .. .. ..@ conf.index      : num 0
  .. .. .. ..@ equiv.sphere    : num 0
  .. .. .. ..@ gradient        : num 0
  .. .. .. ..@ dose.rx         : num 5500
  .. .. .. ..@ dose.fx         : num 0
  .. .. .. ..@ doses           : num [1:1133] 0 5 10 15 20 25 30 35 40 45 ...
  .. .. .. ..@ dose.type       : chr "absolute"
  .. .. .. ..@ volumes         : num [1:1133] 100 77.7 75.1 73 71.3 ...
  .. .. .. ..@ volume.type     : chr "relative"
}
}

\examples{
data(list="RadOnc",package="RadOnc")
}
\keyword{datasets}
