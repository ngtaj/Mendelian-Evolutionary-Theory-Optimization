# Mendelian-Evolutionary-Theory-Optimization
Evolutionary Computation Algorithm: "Mendelian Evolutionary Theory Optimization Algorithm"
To test the Algorithm please run "METOSimulationMain.m" file.
Matlab code of the METO is given in "METO.m" file.
List of Testfunctions are given in "TestFunctionsMETOPaper.m", which is called in "METOSimulationMain.m" by 
"[CostFunction,  VarMin, VarMax] = TestFunctionsMETOPaper(func_num, nVar);" in line 27.

This algorithm is based on the Paper: https://www.techrxiv.org/articles/Mendelian_Evolutionary_Theory_Optimization_Algorithm/12095802

This Optimizers is best suited to solve Mixed variables such as problems having Real, integer and Binary varibales. In the algorithm variables are assigned with the following code:

Real variable = 1
Integer variable = 2
Binary variable = 3 To assign the type of variables please modify the variable "vType" accordingly, where each entry is corresponding to a particular variable. By default it is set to 1, for all real varibales.
All benchmark test functions are from the links: PowerSystemsandEvolutionaryAlgorithms[On-line]http://alroomi.org/benchmarks/unconstrained/n-dimensions.
Benchmarks for Evaluation of Evolutionary Algorithms [On-line]http://www.ntu.edu.sg/home/epnsugan/
