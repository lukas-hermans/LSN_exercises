# LSN Exercise Delivery

This is my (= Lukas Hermans') repository for the exercise delivery for the "Numerical Simulation Laboratory" course at the University of Milan. For detailed information on the course's content visit the dedicated [wep page](https://www.unimi.it/en/education/degree-programme-courses/2021/numerical-simulation-laboratory).

## Structure of the Repository
For each of the 12 exercises I created a folder with the name "exi" where i is the number of the corresponding exercise. Each folder contains a Jupyter-Notebook (e.g. "LSN_Exercises_01.ipynb") in which I present my solution to the exercise. The C++ code that was used to generate the data in the subfolder "exi/data" can be found in the subfolder "exi/cpp_code". There is a "makefile" from which the structure of my C++ files should become obvious.

So, the folder structure for the first exercise is as follows:

* ex1/
    * **<span style="color:red">LSN_Exercises_01.ipynb</span>**
    * data/
        * contains data (usually .txt files)
    * cpp_code/
        * makefile
        * other C++ files