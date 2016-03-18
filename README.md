# AreaCon
A C++ library for area-constrained partitioning, written by Jeffrey R. Peters

##About
AreaCon is a light-weight, C++ library for carrying out area-constrained partitioning operations in floating-point precision. The library is primarily a numerical implementation of the area-constrained partitioning algorithm by Patel, Frasca, and Bullo, 2014 [[Reference](http://www.areacon.org/References)].

Information about AreaCon including documentation, examples, and other tips can be found on this page, the AreaCon [main page](http://www.areacon.org), and the [Wiki](https://github.com/jrpeters/AreaCon/wiki). 

##Components/Dependencies

AreaCon consists of a single compilation unit (one .hpp and one .cpp file), which has dependence on a single free and open source third party library (the Polygon Clipper Library). For convenience, source code for the version of the Clipper library that is required by AreaCon is included with the source code in this repository. More information and the most recent version of the Clipper Library is available [here](http://www.angusj.com/delphi/clipper.php).

##Usage 

AreaCon is free to download and use, subject to the terms and conditions laid out in the license. Please cite AreaCon whenever possible. A bibtex citation is provided below for convenience:
<pre>
@Misc{AreaCon:16,
author =        {J.R. Peters and Contributors},
title =         {The {AreaCon} Library},
howpublished =  {\texttt{http://www.areacon.org}},
year =          {2016},
note =          {v. 1},
abstract =      {A C++ library for area-constrained partitioning.},
}
</pre>

##Contact

Contributions are most easily made through the GitHub site. Emails with questions, etc. can also be sent directly to [the author](www.jeffreyrpeters.com). However, before sending an email, please check the material on the AreaCon [main page](http://www.areacon.org) and the [Wiki](https://github.com/jrpeters/AreaCon/wiki) to be sure that the question/concern/information is not already posted there.
