#############
Instructions:

- The source code is in the folder 'src', with the header files in the folder 'header'. The main program is 'Application.cpp'. It needs to be run in 32bit (x86) mode.

- The .dll files 'ANN.dll' and 'glut32.dll' need to be in the folder from where the program is run.

- The additional libraries are in the folder 'Dependencies'.

- The normalized sample data base is in the folder 'DB'. This folder needs to be in the source directory of the program.

- The dimensionality reduction function is in the Python file 'dr_inderactive.py'. It needs the folder 'DB' to be in the same directory.

- For an interactive plot (the hover functionality) the environment needs to allow for an interactive window. Exampel when using with Spyder:
Tools > Preferences > IPython Console > Graphics > Backend and change it from "Inline" to "Automatic

- The program will use the console to ask the user to either perform a query, extract features or perform an evaluation.
Queries need a path to a model file and a 'outputStand.csv' file in the source directory. Feature extraction needs a path to the sample data base 'DB', same with evaluation. These two functions will also output their respective csv files.

- The csv file 'outputStand.csv' is the extracted standardized features, necessary to perform a query, is already generated, as well as the evaluation csv files.
#############