# exposureTimeCalculator
Links
=====

`Full Documentation <https://github.com/KeckObservatory/exposureTimeCalculator/blob/master/docs/Exposure%20Time%20Calculator%20Documentation.pdf>`

`WMKO Exposure Time Calculator <https://www2.keck.hawaii.edu/inst/PILogin/etcgui/>`

`exposureTimeCalculator GitHub Repository <https://github.com/KeckObservatory/exposureTimeCalculator>`

Requirements
============

We recommend that this project should be run with Anaconda Python 3.5 or above.

The current Python library requirements for this project are:

	astropy
	bokeh [1.0.4]
	flask
	getpass
	importlib
	matplotlib
	numpy
	psutil
	pysynphot
	scipy

How to Download and Run the ETC
===============================

Navigate to the desired location on your computer and clone this repository by running the command:

	git clone https://github.com/KeckObservatory/exposureTimeCalculator.git


Edit line 27 of `etc.py` to change the port to an unoccupied value (if the default of 50008 is unoccupied, there's no need to change this unless you want to). Run the following command to start the ETC:

	python3 etc.py

This will start a webpage on your localhost at the specified port number. In your web browser, enter the following line in the URL bar (substitute new port number if you changed it earlier):

	localhost:50008/etcgui/
