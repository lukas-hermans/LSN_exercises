# LSN Exercise Delivery

This is my (= Lukas Hermans') repository for the exercise delivery for the "Numerical Simulation Laboratory" course at the University of Milan. For detailed information on the course's content visit the dedicated [wep page](https://www.unimi.it/en/education/degree-programme-courses/2021/numerical-simulation-laboratory).

## Structure of the Repository
The 12 original exercise sheets (distributed by the professor) can be found in the folder "ex_descriptions". For each of the 12 exercise sheets, I created a folder with the name "exi" where i is the number of the corresponding exercise sheet that contain my solution to the exercise sheets. 

Each folder contains a Jupyter-Notebook (e.g. "LSN_Exercises_01.ipynb"), a sub-folder "cpp_code", and a sub-folder "data". In the Jupyter-Notebook, I present my solution to the exercise sheet with figures and detailed explanations. The sub-folder "cpp_code" includes all the C++ code used for the calculations in each exercise sheet. It also contains a "makefile", from which the hierarchy of the C++ and header files should become obvious. The generated data is saved in the sub-folder "data", which is accessed by the Jupyter-Notebook. The last two exercise sheets 11 and 12 do not contain the subfolders because all the computations (with Python and Keras) are included in the Jupyter-Notebook.

In addition to the folders of the 12 exercise sheets, I created a folder called "tools" that includes libraries (e.g., the random number generated or routines for the Metropolis algorithm) that are used in several exercise sheets.


