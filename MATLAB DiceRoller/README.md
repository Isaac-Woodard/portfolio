# Overview
This project provides classes and functions for rolling virtual dice:
- Make numerical dice with NumericDie.m
- Make text dice with Categoric.m

Being non-physical, dice can easily be defined with any number of sides. A weight parameter can also be set to make some sides favored over others!

Some basic usage examples are included in Examples.m.

# Contents
- +probability - A package of custom probability-related functions.
- tests/ - Directory for tests. Test are organized into files by class and package.
- CategoricDie.m - A class for a text based die. Stores data in categoric arrays.
- Die.m - An abstract base class for a die. Inherited by CategoricDie.m and NumericDie.m.
- Examples.m - A MATLAB script with some usage examples.
- GMRoller.mlapp - An App Designer app for assisting game masters with their dice rolling.
- NumericDie.m - A class for a numeral die. 
- PercentileDice.m - A class for a pair of d10s. Not implemented.
- test_suite.m - A MATLAB script for creating and running a test suite from the files in tests/.

# Object-Oriented Design
Object-oriented programming is a powerful way to keep a project organized and scalable. This project demonstrates many of MATLABs OOP features including:
- Classes
- Properties and attributes
- Static methods
- Inheritance

Also demonstrated are some of MATLABs features for maintaining code quality:
- Input validation with arguments blocks
- Packages
- Unit testing

# GMRoller
GMRoller is an app for helping roleplaying game masters with their dice rolling. Unfortunately, the MATLAB Compiler isn't available with a Home license, so to see the app in action you'll need a MATLAB license yourself.

The app features the five most common numeric dice along with a d2 and a d100.
- d2, d4, d6, d10, d12, d20, d100

With physical dice a pair of d10s are typically used to provide a range from 1-100, but with virtual dice it's easier to just make a d100! 

In addition to numeric dice, there are also a few special dice with text options.
- Weather: Clear, Clouds, Rain, Snow
- Direction: North, East, South, West
- Mood: Happy, Sad, Angry, Excited, Anxious

# References
- [MATLABs core random number functions](https://www.mathworks.com/help/matlab/math/create-arrays-of-random-numbers.html)
- [Class based unit tests](https://www.mathworks.com/help/matlab/class-based-unit-tests.html?s_tid=CRUX_lftnav)
- [Some very tidy solutions to generating weighted random numbers](https://stackoverflow.com/questions/2977497/weighted-random-numbers-in-matlab)
- [Documentation for the function I actually used](https://www.mathworks.com/help/matlab/ref/histcounts.html)
- [Odd-sided dice for sale](https://mathartfun.com/d357.html)

# Ideas
- Add icons for app title and numeric dice.
- PercentileDice: two d10s bundled together to make a pseudo d100.
- ContinuousDie: a die with an infinite number of sides which have values within a certain range. Nonphysical but might be interesting...Could also consider allowing multiple ranges, e.g. a die which can roll 1 to 5 or 10 to 15. Useful or not, weighting would be very interesting.
- SymbolicDie: for images instead of text. Uses could include a mood die for emojis or a weather die with various weather icons. It seems a categoric die accomplishs all the same things, though, and is much easier to use for any form of data analysis.
- 3D die simulator using Simulink...compare random and chaotic...single die builder.

# TODO
- Fix bug where table clears on right clicking before selecting context menu option.