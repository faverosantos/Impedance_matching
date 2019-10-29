# RF library
A pythonic and very simplistic way of doing what Iowa Hills' Smith Chart does. \
\
It surely has no same interface nor functionality as Iowa Hills' software does, but this python program allows the user to better visualize constant reactance and susceptance movements by drawing real time curves as the mouse moves.

## Why am I doing this?
This is being used for my Ph.D. studies as I needed to better understand the movements of lumped components on the Smith Chart and how to realize an impedance match.
Final objective is to make this a standard RF library for RF designers.

## How to use:
I. Mouse right button: set target impedance ('x' marker) \
II. Mouse left button: set de initial/intermediate impedances (circle markers). By cascating left mouse button clicks, impedance movements are "real time" traced, allowing user add series or shunt capacitors/inductors and to see how the movement on the Smith Chart happens.\
III. Mouse center button: toggle show mode on and off. Show mode allows to draw all movements until the last lef click. \
IV. By pressing "1", "2", "3" or "4" on the keyboard, user can select which element type and connection to add: series capacitance, shunt capacitance, series inductance or shunt inducatance, respectivelly.

Note it's important to call "set_operation_frequency(my_chosen_frequency)" before "plot_smith_chart()" as this sets the frequency which all impedances will be calculated.

## Observations

## To do
I. Add non-ideal passive models
II. Optimize movement plot
III. Add a GUI to show how the circuit is being created (similar to the Iowa Hills SW)
IV. Add a clear all command (preferible by pressing c in the keyboard)
