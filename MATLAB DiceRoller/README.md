# Overview
This project provides classes and functions for rolling virtual dice:
- Make numerical dice with NumericDie.m
- Make text dice with Categoric.m

Being non-physical, dice can easily be defined with any number of sides. A weight parameter can also be set to make some sides favored over others!

Some basic usage examples are included in Examples.m.

# Object-Oriented Design...
Demonstrates OOP design in MATLAB including:
- App design
- Input validation with arguments blocks
- Properties and attributes
- Inheritance
- Packages
- Static methods
- Unit tests and test suites

# GMRoller...

# Ideas
- PercentileDice: two d10s bundled together to make a pseudo d100.
- ContinuousDie: a die with an infinite number of sides which have values within a certain range. Nonphysical but might be interesting...Could also consider allowing multiple ranges, e.g. a die which can roll 1 to 5 or 10 to 15. Useful or not, weighting would be very interesting.
- SymbolicDie: for images instead of text. Uses could include a mood die for emojis or a weather die with various weather icons. It seems a categoric die accomplishs all the same things, though, and is much easier to use for any form of data analysis.
- 3D die simulator using Simulink...compare random and chaotic...single die builder.

# References
- https://www.mathworks.com/help/matlab/math/create-arrays-of-random-numbers.html
- https://www.mathworks.com/help/matlab/class-based-unit-tests.html?s_tid=CRUX_lftnav
- https://stackoverflow.com/questions/2977497/weighted-random-numbers-in-matlab
- https://www.mathworks.com/help/matlab/ref/histcounts.html
- [Odd-sided dice for sale](https://mathartfun.com/d357.html)