# exposureTimeCalculator
Links
=====

`Full Documentation <https://github.com/KeckObservatory/exposureTimeCalculator/blob/master/docs/Exposure%20Time%20Calculator%20Documentation.pdf>`

`WMKO Exposure Time Calculator <https://www2.keck.hawaii.edu/inst/PILogin/etcgui/>`

`exposureTimeCalculator GitHub Repository <https://github.com/KeckObservatory/exposureTimeCalculator>`

Requirements
============

This project is written in Python 3 and should be run with Python 3.5 or above.

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

The file structure will be as follows:

	./exposureTimeCalculator/*.py
				/static
				/templates
				/datafiles

Edit line 27 of `etc.py` to change the port to an unoccupied value (if the default of 50008 is unoccupied, there's no need to change this unless you want to). Then edit line 52 of :code:`manager.py` to point to your local install of python3. The next commands are dependent on aliasing the local python3 executable to "python3". It can be substituted with the full path to that executable instead if desired.

Afterwards, run the following command in the exposureTimeCalculator directory:

	 python3 manager.py etc start

You can also run the equivalent command:

	python3 etc.py

(More details on the `manager.py` script are provided on the `Using the ETC` section.)

This will start a webpage on your localhost at the specified port number.

In your web browser, enter the following line in the URL bar (substitute new port number if you changed it earlier):

	localhost:50008/etcgui/
